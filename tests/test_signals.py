# SPDX-FileCopyrightText: The NOMAD Authors
# SPDX-License-Identifier: Apache-2.0

"""Focused reader data contracts."""

from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np
import pint
import pytest

from fairmat_readers_xrd import read_bruker_brml, read_panalytical_xrdml
from fairmat_readers_xrd.timing import counting_time, rasx_counting_time
from fairmat_readers_xrd.utils import detect_scan_type, modify_scan_data


def intensity_xrdml_scan(unit, values='2 4', times='2 3', omega=1):
    return f"""<scan scanAxis="2Theta"><dataPoints>
      <countingTimes unit="seconds">{times}</countingTimes>
      <intensities unit="{unit}">{values}</intensities>
      <positions axis="2Theta" unit="deg"><startPosition>10</startPosition>
      <endPosition>20</endPosition></positions>
      <positions axis="Omega" unit="deg"><commonPosition>{omega}</commonPosition>
      </positions></dataPoints></scan>"""


def intensity_write_xrdml(tmp_path, scans):
    path = tmp_path / 'scan.xrdml'
    path.write_text(
        '<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.1"><xrdMeasurement>'
        + ''.join(scans)
        + '</xrdMeasurement></xrdMeasurements>'
    )
    return path


@pytest.mark.parametrize(
    ('unit', 'target', 'expected'),
    [
        ('counts', 'dimensionless', [2, 4]),
        ('counts per second', 'count/s', [2, 4]),
        ('counts/minute', 'count/s', [2 / 60, 4 / 60]),
    ],
)
def test_xrdml_declared_signal_scale_is_decoded(tmp_path, unit, target, expected):
    data = read_panalytical_xrdml(
        str(intensity_write_xrdml(tmp_path, [intensity_xrdml_scan(unit)]))
    )
    np.testing.assert_allclose(data['intensity'].to(target).magnitude, expected)
    np.testing.assert_allclose(data['countTime'].to('s').magnitude, [2, 3])


@pytest.mark.parametrize('unit', ['watts', 'not_a_signal_unit'])
def test_unknown_signal_units_are_not_treated_as_counts(tmp_path, unit):
    with pytest.raises(ValueError, match='XRD signal unit'):
        read_panalytical_xrdml(
            str(intensity_write_xrdml(tmp_path, [intensity_xrdml_scan(unit)]))
        )


def test_mixed_count_and_rate_scans_are_not_combined(tmp_path):
    path = intensity_write_xrdml(
        tmp_path, [intensity_xrdml_scan('counts'), intensity_xrdml_scan('cps', omega=2)]
    )
    with pytest.raises(ValueError, match='Mixed counts and count-rate'):
        read_panalytical_xrdml(str(path))


ureg = pint.get_application_registry()


@pytest.mark.parametrize('channels', [1, 7])
@pytest.mark.parametrize('rows', [1, 3])
def test_constant_signals_keep_all_points_and_repeated_scans(channels, rows):
    scan = {
        'intensity': [np.full(channels, 5) * ureg.dimensionless] * rows,
        'counts': [np.full(channels, 2) * ureg.dimensionless] * rows,
        '2Theta': [np.arange(channels) * ureg.degree] * rows,
        'Omega': [np.full(channels, i % 2) * ureg.degree for i in range(rows)],
        'Chi': None,
    }
    kind = detect_scan_type(scan)
    assert kind == ('line' if rows == 1 else 'rsm')
    result = modify_scan_data(scan, kind)
    shape = (channels,) if rows == 1 else (rows, channels)
    np.testing.assert_array_equal(result['intensity'].magnitude, np.full(shape, 5))
    np.testing.assert_array_equal(result['counts'].magnitude, np.full(shape, 2))
    if rows > 1:
        np.testing.assert_array_equal(result['Omega'].magnitude, [0, 1, 0])


def test_brml_omega_matches_recorded_scan_coordinates(tmp_path):
    path = Path(__file__).parent / 'data/23-012-AG_2thomegascan_long.brml'
    with ZipFile(path) as source:
        root = ET.fromstring(source.read('Experiment0/RawData0.xml'))
    route = root.find('./DataRoutes/DataRoute')
    view = next(
        v
        for v in route.findall('./DataViews/RawDataView')
        if v.find('./Varying[@LogicName="ScanAxes"]') is not None
    )
    fields = view.findall('./Varying/FieldDefinitions')
    index = int(view.get('Start')) + next(
        (i for i, f in enumerate(fields) if f.get('FieldName') == 'Theta')
    )
    recorded = np.array([row.text.split(',') for row in route.findall('Datum')], float)
    result = read_bruker_brml(str(path))
    np.testing.assert_array_equal(
        result['Omega'].to('degree').magnitude, recorded[:, index]
    )
    recorded[10, index] += 0.0002
    route.findall('Datum')[10].text = ','.join(map(str, recorded[10]))
    changed = tmp_path / 'irregular.brml'
    with ZipFile(path) as source, ZipFile(changed, 'w') as target:
        for name in source.namelist():
            target.writestr(
                name,
                ET.tostring(root)
                if name == 'Experiment0/RawData0.xml'
                else source.read(name),
            )
    result = read_bruker_brml(str(changed))
    np.testing.assert_array_equal(
        result['Omega'].to('degree').magnitude, recorded[:, index]
    )


@pytest.mark.parametrize(
    ('mode', 'speed', 'unit', 'step', 'expected'),
    [
        ('STEP', 1.2, 'sec', 0.02, 1.2),
        ('STEP', 0.1, 'minute', 0.02, 6),
        ('CONTINUOUS', 0.04, 'deg/min', 0.0048, 7.2),
        ('CONTINUOUS', 1, 'deg/min', -0.04, 2.4),
    ],
)
def test_rasx_modes_have_distinct_time_meanings(mode, speed, unit, step, expected):
    result = rasx_counting_time(
        dict(Mode=mode, Speed=speed, SpeedUnit=unit, Step=step, PositionUnit='deg')
    )
    assert result.to('s').magnitude == pytest.approx([expected])


@pytest.mark.parametrize('value', [None, 0, np.inf])
def test_invalid_optional_times_are_omitted(value):
    assert counting_time(value, 's') is None
