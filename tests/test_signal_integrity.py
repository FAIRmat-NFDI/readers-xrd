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
"""Measurement multiplicity and optional scan metadata are independent."""

from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np
import pint
import pytest

from fairmat_readers_xrd import read_bruker_brml, read_panalytical_xrdml
from fairmat_readers_xrd.utils import detect_scan_type, modify_scan_data

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


def test_shared_coordinates_are_compared_and_stored_in_common_units():
    scan = {
        'intensity': [np.array([1, 2]) * ureg.dimensionless] * 2,
        '2Theta': [
            np.array([0, 180]) * ureg.degree,
            np.array([0, np.pi]) * ureg.radian,
        ],
        'Omega': [np.array([1]) * ureg.degree, np.array([2]) * ureg.degree],
    }
    assert detect_scan_type(scan) == 'rsm'
    result = modify_scan_data(scan, 'rsm')
    np.testing.assert_allclose(result['2Theta'].to('degree').magnitude, [0, 180])


@pytest.mark.parametrize('reverse', [False, True])
def test_mixed_optional_xrdml_fields_do_not_break_scan_order(tmp_path, reverse):
    scans = [
        '<intensities>5 5 5</intensities>',
        '<counts>2 2 2</counts><beamAttenuationFactors>3 3 3</beamAttenuationFactors>',
    ]
    if reverse:
        scans.reverse()
    body = ''.join(
        f'<scan scanAxis="2Theta"><dataPoints>{signal}'
        '<positions axis="2Theta" unit="deg"><startPosition>10</startPosition>'
        '<endPosition>30</endPosition></positions>'
        f'<positions axis="Omega" unit="deg"><commonPosition>{i}</commonPosition>'
        '</positions></dataPoints></scan>'
        for i, signal in enumerate(scans)
    )
    path = tmp_path / 'mixed.xrdml'
    path.write_text(
        '<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.1">'
        f'<xrdMeasurement>{body}</xrdMeasurement></xrdMeasurements>'
    )
    result = read_panalytical_xrdml(str(path))
    expected = np.array([[5, 5, 5], [6, 6, 6]])
    np.testing.assert_array_equal(
        result['intensity'].magnitude, expected[::-1] if reverse else expected
    )
    assert result.get('counts') is None
    assert result.get('beamAttenuationFactors') is None


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
        i for i, f in enumerate(fields) if f.get('FieldName') == 'Theta'
    )
    recorded = np.array([row.text.split(',') for row in route.findall('Datum')], float)
    result = read_bruker_brml(str(path))
    np.testing.assert_array_equal(
        result['Omega'].to('degree').magnitude, recorded[:, index]
    )
    # Actual point positions need not equal the nominal start/stop raster.
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
