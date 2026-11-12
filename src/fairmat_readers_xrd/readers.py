#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import collections
import xml.etree.ElementTree as ET
from typing import TYPE_CHECKING, Any, Dict

import numpy as np
import pint

from fairmat_readers_xrd.bruker_raw_parser import BrukerRAW4Parser
from fairmat_readers_xrd.ikz import BRMLfile, RASXfile
from fairmat_readers_xrd.models import validate_output
from fairmat_readers_xrd.metadata import (
    brml_metadata,
    rasx_metadata,
    raw_metadata,
    xrdml_metadata,
)
from fairmat_readers_xrd.timing import (
    complete_scan_times,
    counting_time,
    rasx_counting_time,
)

# from pynxtools.dataconverter.convert import transfer_data_into_template
from fairmat_readers_xrd.utils import (
    detect_scan_type,
    intensity_quantity,
    modify_scan_data,
)

if TYPE_CHECKING:
    from structlog.stdlib import (
        BoundLogger,
    )

ureg = pint.get_application_registry()


def transfer_data_into_template(**kwargs):
    raise NotImplementedError


def read_file(file_path: str) -> dict:
    """
    Main function to figure out which specific file format to read.

    Args:
        file_path (str): The path to the file to be read.

    Returns:
        Dict[str, Any]: The X-ray diffraction data in a Python dictionary.
    """
    if file_path.endswith('.rasx'):
        return read_rigaku_rasx(file_path)
    if file_path.endswith('.xrdml'):
        return read_panalytical_xrdml(file_path)
    if file_path.endswith('.brml'):
        return read_bruker_brml(file_path)
    if file_path.endswith('.raw'):
        return read_bruker_raw(file_path)
    raise NotImplementedError('Unknown file type.')


def _common_scan_value(values, field, logger):
    values = set(values)
    if len(values) == 1:
        return values.pop()
    if logger is not None and values:
        logger.warning(f'Mixed scan {field}; omitting common metadata.')
    return None


@validate_output
def read_panalytical_xrdml(
    file_path: str, logger: 'BoundLogger' = None
) -> dict[str, Any]:
    """
    Function for reading the X-ray diffraction data in a Panalytical `.xrdml` file.

    Args:
        file_path (str): The path to the `.xrdml` file.
        logger (BoundLogger): A structlog logger.

    Returns:
        Dict[str, Any]: The X-ray diffraction data in a Python dictionary.
    """
    with open(file_path, encoding='utf-8') as file:
        element_tree = ET.parse(file)
    root = element_tree.getroot()
    ns_version = root.tag.split('}')[0].strip('{')
    ns = {'xrdml': ns_version}

    def find_float(path):
        result = root.find(path, ns)
        if result is not None:
            unit = result.get('unit', '')
            if unit == 'Angstrom':
                unit = 'angstrom'
            return float(result.text) * ureg(unit)
        return None

    def find_string(path):
        result = root.find(path, ns)
        if result is not None:
            return result.text
        return None

    xrd_measurement = root.find('xrdml:xrdMeasurement', ns)
    scans = xrd_measurement.findall('xrdml:scan', ns)

    axes = collections.defaultdict(list)
    counting_time_value = []
    attenuation_values = []
    intensities_array = []
    counts_array = []

    for scan_index, scan in enumerate(scans):
        data_points = scan.find('xrdml:dataPoints', ns)
        intensities = data_points.find('xrdml:intensities', ns)
        time_element = data_points.find('xrdml:commonCountingTime', ns)
        if time_element is None:
            time_element = data_points.find('xrdml:countingTimes', ns)
        attenuation = data_points.find('xrdml:beamAttenuationFactors', ns)

        counting_time_value.append(
            counting_time(time_element.text, time_element.get('unit'), logger)
            if time_element is not None
            else None
        )
        if attenuation is not None:
            attenuation_values.append(
                np.fromstring(attenuation.text, sep=' ') * ureg.dimensionless
            )
            scaling_factor = attenuation_values[-1]
        else:
            attenuation_values.append(None)
            scaling_factor = 1
        if intensities is not None:
            intensities_array.append(
                intensity_quantity(
                    np.fromstring(intensities.text, sep=' '), intensities.get('unit')
                )
            )
            counts_array.append(None)
        else:
            counts = data_points.find('xrdml:counts', ns)
            counts_array.append(
                intensity_quantity(
                    np.fromstring(counts.text, sep=' '), counts.get('unit')
                )
            )
            if not counts_array[-1].is_compatible_with('dimensionless'):
                raise ValueError('XRDML counts must have count units, not rate units.')
            intensities_array.append(
                counts_array[-1] * scaling_factor * ureg.dimensionless
            )
        n_points = len(intensities_array[-1])

        for axis in data_points.findall('xrdml:positions', ns):
            name = axis.get('axis')
            if name not in axes:
                axes[name] = [None] * len(scans)
            unit = axis.get('unit')
            listed = axis.find('xrdml:listPositions', ns)
            start = axis.find('xrdml:startPosition', ns)
            end = axis.find('xrdml:endPosition', ns)
            common = axis.find('xrdml:commonPosition', ns)
            if listed is not None:
                axes[name][scan_index] = np.fromstring(listed.text, sep=' ') * ureg(
                    unit
                )
            elif start is not None and end is not None:
                axes[name][scan_index] = np.linspace(
                    float(start.text), float(end.text), n_points
                ) * ureg(unit)
            elif common is not None:
                axes[name][scan_index] = np.array([float(common.text)]) * ureg(unit)
            elif logger is not None:
                logger.warning(f'Unknown data for {name} axis.')

    scan_data = collections.defaultdict(list)
    scan_data.update(axes)
    if any(
        not signal.is_compatible_with(intensities_array[0])
        for signal in intensities_array[1:]
    ):
        raise ValueError('Mixed counts and count-rate scans are not supported.')
    scan_data['intensity'] = intensities_array
    scan_data['counts'] = (
        counts_array if all(value is not None for value in counts_array) else None
    )
    scan_data['countTime'] = complete_scan_times(counting_time_value, logger)
    scan_data['beamAttenuationFactors'] = (
        attenuation_values
        if all(value is not None for value in attenuation_values)
        else None
    )

    scan_type = detect_scan_type(scan_data)
    modified_scan_data = modify_scan_data(scan_data, scan_type)

    return {
        **modified_scan_data,
        'countTime': modified_scan_data.get('countTime'),
        'scanmotname': _common_scan_value(
            (s.get('scanAxis') for s in scans), 'axis', logger
        ),
        'metadata': {
            'details': xrdml_metadata(root, ns, logger),
            'sample_id': find_string('xrdml:sample/xrdml:id'),
            'measurement_type': xrd_measurement.get('measurementType', None),
            'sample_mode': xrd_measurement.get('sampleMode', None),
            'instrument_id': find_string(
                'xrdml:xrdMeasurement/xrdml:scan/xrdml:header/xrdml:source'
                '/xrdml:instrumentID'
            ),
            'start_time': scans[0].findtext(
                'xrdml:header/xrdml:startTimeStamp', namespaces=ns
            ),
            'end_time': scans[-1].findtext(
                'xrdml:header/xrdml:endTimeStamp', namespaces=ns
            ),
            'source': {
                'current': find_float(
                    'xrdml:xrdMeasurement/xrdml:incidentBeamPath/xrdml:xRayTube'
                    '/xrdml:current'
                ),
                'voltage': find_float(
                    'xrdml:xrdMeasurement/xrdml:incidentBeamPath/xrdml:xRayTube'
                    '/xrdml:tension'
                ),
                'kAlpha1': find_float(
                    'xrdml:xrdMeasurement/xrdml:usedWavelength/xrdml:kAlpha1'
                ),
                'kAlpha2': find_float(
                    'xrdml:xrdMeasurement/xrdml:usedWavelength/xrdml:kAlpha2'
                ),
                'kBeta': find_float(
                    'xrdml:xrdMeasurement/xrdml:usedWavelength/xrdml:kBeta'
                ),
                'ratioKAlpha2KAlpha1': find_float(
                    'xrdml:xrdMeasurement/xrdml:usedWavelength/xrdml:ratioKAlpha2KAlpha1'
                ),
                'anode_material': find_string(
                    'xrdml:xrdMeasurement/xrdml:incidentBeamPath/xrdml:xRayTube'
                    '/xrdml:anodeMaterial'
                ),
            },
            'scan_mode': _common_scan_value(
                (s.get('mode') for s in scans), 'mode', logger
            ),
            'scan_axis': _common_scan_value(
                (s.get('scanAxis') for s in scans), 'axis', logger
            ),
            'scan_type': scan_type,
        },
    }


@validate_output
def read_rigaku_rasx(file_path: str, logger: 'BoundLogger' = None) -> dict[str, Any]:
    """
    Reads .rasx files from Rigaku instruments
        - reader is based on IKZ module

    Args:
        file_path (string): absolute path of the file.
        logger (BoundLogger): A structlog logger for propagating errors and warnings.

    Returns:
        Dict[str, Any]: The X-ray diffraction data in a Python dictionary.
    """
    reader = RASXfile(file_path, verbose=False)
    scan_data = reader.get_scan_data(logger)
    source = reader.get_source_info()

    scan_infos = [
        meta.get('ScanInformation', {}) for meta in reader.meta[: len(reader.data)]
    ]
    scan_data['countTime'] = complete_scan_times(
        [rasx_counting_time(info, logger) for info in scan_infos], logger
    )
    scan_type = detect_scan_type(scan_data)
    modified_scan_data = modify_scan_data(scan_data, scan_type)
    scan_axis = _common_scan_value(
        (s.get('AxisName') for s in scan_infos), 'axis', logger
    )

    return {
        **modified_scan_data,
        'countTime': modified_scan_data.get('countTime'),
        'metadata': {
            'details': rasx_metadata(reader, logger),
            'sample_id': None,
            'scan_axis': scan_axis,
            'scan_type': scan_type,
            'scan_mode': _common_scan_value(
                (s.get('Mode') for s in scan_infos), 'mode', logger
            ),
            'start_time': scan_infos[0].get('StartTime') if scan_infos else None,
            'end_time': scan_infos[-1].get('EndTime') if scan_infos else None,
            'source': {
                'anode_material': source.get('TargetName'),
                'kAlpha1': source.get('WavelengthKalpha1'),
                'kAlpha2': source.get('WavelengthKalpha2'),
                'kBeta': source.get('WavelengthKbeta'),
                'voltage': source.get('Voltage'),
                'current': source.get('Current'),
                'ratioKAlpha2KAlpha1': None,
            },
        },
    }


@validate_output
def read_bruker_brml(file_path: str, logger: 'BoundLogger' = None) -> dict[str, Any]:
    """
    Reads .brml files from Bruker instruments
        - reader is based on IKZ module

    Args:
        file_path (string): absolute path of the file.
        logger (BoundLogger): A structlog logger for propagating errors and warnings.

    Returns:
        Dict[str, Any]: The X-ray diffraction data in a Python dictionary.
    """
    reader = BRMLfile(file_path, verbose=False)
    scan_info = reader.get_scan_info()
    source = reader.get_source_info()
    scan_data = reader.get_scan_data(logger)
    scan_data['countTime'] = complete_scan_times(
        [
            counting_time(value, 's', logger)
            for value in reader.data['TimePerStepEffective']
        ],
        logger,
    )
    scan_type = detect_scan_type(scan_data)
    modified_scan_data = modify_scan_data(scan_data, scan_type)
    modified_scan_data['Omega'] = modified_scan_data.pop(
        'Theta', None
    )  # theta and omega are synonymous in .brml

    return {
        **modified_scan_data,
        'countTime': modified_scan_data.get('countTime'),
        'metadata': {
            'details': brml_metadata(reader, logger),
            'sample_id': None,
            'scan_axis': _common_scan_value(reader.data['ScanName'], 'axis', logger),
            'scan_type': scan_type,
            'scan_mode': _common_scan_value(reader.data['ScanMode'], 'mode', logger),
            'start_time': scan_info.get('start_time'),
            'end_time': scan_info.get('end_time'),
            'source': {
                'anode_material': source.get('TubeMaterial'),
                'kAlpha1': source.get('WaveLengthAlpha1'),
                'kAlpha2': source.get('WaveLengthAlpha2'),
                'kBeta': source.get('WaveLengthBeta'),
                'ratioKAlpha2KAlpha1': source.get('WaveLengthRatio'),
                'voltage': source.get('Voltage'),
                'current': source.get('Current'),
            },
        },
    }


@validate_output
def read_bruker_raw(file_path: str, logger: 'BoundLogger' = None) -> Dict[str, Any]:
    """Read the validated single-range RAW4 powder geometry and native metadata."""
    parser = BrukerRAW4Parser(file_path)
    if logger is not None:
        parser.logger = logger
    raw = parser.parse()
    meta, params = raw['metadata'], raw['scan_params']
    source = {
        'anode_material': meta.get('anode_material'),
        'kAlpha1': meta.get('wavelength_kalpha1'),
        'kAlpha2': meta.get('wavelength_kalpha2'),
        'kBeta': meta.get('wavelength_kbeta'),
        'ratioKAlpha2KAlpha1': meta.get('kalpha2_kalpha1_ratio'),
        'voltage': meta['voltage'] * ureg.kilovolt,
        'current': meta['current'] * ureg.milliampere,
    }
    return {
        '2Theta': np.asarray(parser.angles) * ureg.degree,
        'intensity': np.asarray(parser.intensities) * ureg.dimensionless,
        'counts': None,
        'countTime': counting_time(params.get('counting_time_ms'), 'ms', logger),
        'beamAttenuationFactors': None,
        'scanmotname': params['scan_axis'],
        'metadata': dict(
            details=raw_metadata(meta),
            sample_id=meta.get('sample_id'),
            scan_type='line',
            scan_axis=params['scan_axis'],
            scan_mode=params['scan_mode'],
            source=source,
        ),
    }
