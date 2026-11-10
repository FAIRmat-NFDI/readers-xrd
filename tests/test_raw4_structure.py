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
"""RAW4 axes and samples follow declared metadata, not fixture offsets."""

import struct
from pathlib import Path
from types import SimpleNamespace
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np
import pytest

from fairmat_readers_xrd import (
    read_bruker_raw,
    read_panalytical_xrdml,
    read_rigaku_rasx,
)
from fairmat_readers_xrd.bruker_raw_parser import BrukerRAW4Parser
from fairmat_readers_xrd.ikz import RASXfile
from fairmat_readers_xrd.utils import detect_scan_type, modify_scan_data


def raw4(points=5, record_size=4, extra_size=0, mode='Locked Coupled'):
    preamble = bytearray(61)
    preamble[:8] = b'RAW4.00\x00'
    metadata = bytearray(8 + extra_size)
    struct.pack_into('<II', metadata, 0, 99, len(metadata))
    source = bytearray(136)
    struct.pack_into('<II', source, 0, 30, len(source))
    struct.pack_into('<5d', source, 72, 1.541, 1.541, 0, 1.392, 0)
    source[116:120] = b'Cu\x00\x00'
    scan = bytearray(160)
    scan[32 : 32 + len(mode)] = mode.encode()
    struct.pack_into('<2dIf', scan, 72, 10, 0.01, points, 250)
    struct.pack_into('<2f', scan, 100, 40, 30)
    struct.pack_into('<II', scan, 136, record_size, 64)
    drive = bytearray(64)
    struct.pack_into('<III', drive, 0, 50, 64, 2)
    drive[12:18] = b'2Theta'
    struct.pack_into('<d', drive, 56, 10)
    samples = bytearray(points * record_size)
    for i in range(points):
        struct.pack_into('<f', samples, i * record_size, i + 1)
    return bytes(preamble + metadata + source + scan + drive + samples)


@pytest.mark.parametrize('record_size', [4, 8, 12])
@pytest.mark.parametrize('extra_size', [0, 37])
def test_raw4_uses_declared_sizes(tmp_path, record_size, extra_size):
    path = tmp_path / 'scan.raw'
    path.write_bytes(raw4(record_size=record_size, extra_size=extra_size))
    result = read_bruker_raw(str(path))
    np.testing.assert_array_equal(result['intensity'].magnitude, [1, 2, 3, 4, 5])
    np.testing.assert_allclose(result['2Theta'].magnitude, 10 + np.arange(5) * 0.01)
    np.testing.assert_allclose(result['countTime'].to('s').magnitude, [0.25])
    assert result['metadata']['scan_axis'] == '2Theta'
    assert result['metadata']['source']['kAlpha1'] == pytest.approx(1.541)
    assert result['metadata']['source']['kAlpha2'] is None
    assert result['metadata']['source']['ratioKAlpha2KAlpha1'] == 0
    assert result['metadata']['source']['voltage'].to('V').magnitude == pytest.approx(
        40000
    )


def test_raw4_point_count_is_not_truncated_to_16_bits(tmp_path):
    path = tmp_path / 'large.raw'
    points = 65537
    path.write_bytes(raw4(points=points))
    result = BrukerRAW4Parser(str(path)).parse()
    assert len(result['intensities']) == points
    assert result['intensities'][-1] == points


@pytest.mark.parametrize('bad', ['mode', 'truncated', 'block_length'])
def test_raw4_rejects_unsupported_or_incomplete_data(tmp_path, bad):
    content = raw4(mode='Rocking Curve' if bad == 'mode' else 'Locked Coupled')
    if bad == 'truncated':
        content = content[:-1]
    if bad == 'block_length':
        content = bytearray(content)
        struct.pack_into('<I', content, 65, 100000)
    path = tmp_path / 'bad.raw'
    path.write_bytes(content)
    with pytest.raises(ValueError):
        read_bruker_raw(str(path))


@pytest.mark.parametrize('points', [1, 2, 5])
def test_rasx_profile_preserves_point_and_column_dimensions(tmp_path, points):
    source = Path(__file__).parent / 'data/TwoTheta_scan_powder.rasx'
    with ZipFile(source) as archive:
        profile = next(n for n in archive.namelist() if 'Profile' in n)
        conditions = profile.replace('Profile', 'MesurementConditions')[:-4] + '.xml'
        root = ET.fromstring(archive.read(conditions))
    angles = 15 + np.arange(points) * 0.25
    intensity = 2.5 + np.arange(points)
    info = root.find('ScanInformation')
    for field, value in dict(
        Mode='STEP',
        Start=angles[0],
        Stop=angles[-1],
        Step=0.25,
        Speed=1,
        SpeedUnit='s',
        DataCount=points,
    ).items():
        info.find(field).text = str(value)
    for axis in root.find('Axes'):
        if axis.get('Name') == 'TwoTheta':
            axis.set('Position', str(angles[0]))
            axis.set('EndPosition', str(angles[-1]))
    path = tmp_path / 'profile.rasx'
    with ZipFile(path, 'w') as archive:
        rows = ''.join(
            f'{angle}\t{value}\t1\n' for angle, value in zip(angles, intensity)
        )
        archive.writestr(profile, rows.encode('utf-8-sig'))
        archive.writestr(conditions, ET.tostring(root))
    result = read_rigaku_rasx(str(path))
    assert result['metadata']['scan_type'] == 'line'
    assert result['intensity'].shape == (points,)
    assert result['2Theta'].shape == (points,)
    np.testing.assert_array_equal(result['intensity'].magnitude, intensity)
    np.testing.assert_array_equal(result['2Theta'].to('degree').magnitude, angles)


def test_rasx_omega_profile_is_not_two_theta(tmp_path):
    source = Path(__file__).parent / 'data/TwoTheta_scan_powder.rasx'
    with ZipFile(source) as archive:
        profile = next(n for n in archive.namelist() if 'Profile' in n)
        conditions = profile.replace('Profile', 'MesurementConditions')[:-4] + '.xml'
        root = ET.fromstring(archive.read(conditions))
    root.find('./ScanInformation/AxisName').text = 'Omega'
    for axis in root.find('Axes'):
        if axis.get('Name') == 'TwoTheta':
            axis.set('Position', '44')
    path = tmp_path / 'omega.rasx'
    with ZipFile(path, 'w') as archive:
        archive.writestr(profile, '15\t2\t1\n16\t3\t1\n17\t4\t1\n'.encode('utf-8-sig'))
        archive.writestr(conditions, ET.tostring(root))
    result = read_rigaku_rasx(str(path))
    np.testing.assert_array_equal(result['Omega'].to('degree').magnitude, [15, 16, 17])
    np.testing.assert_array_equal(result['2Theta'].to('degree').magnitude, [44])


def test_partial_angle_does_not_become_shared_metadata(tmp_path):
    body = ''.join(
        '<scan><dataPoints><counts>1 2 3</counts>'
        '<positions axis="2Theta" unit="deg"><startPosition>10</startPosition>'
        '<endPosition>30</endPosition></positions>'
        f'<positions axis="Omega" unit="deg"><commonPosition>{i}</commonPosition>'
        '</positions>'
        + (
            '<positions axis="Chi" unit="deg"><commonPosition>5</commonPosition>'
            '</positions>'
            if i
            else ''
        )
        + '</dataPoints></scan>'
        for i in range(2)
    )
    path = tmp_path / 'partial.xrdml'
    path.write_text(
        '<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.1">'
        f'<xrdMeasurement>{body}</xrdMeasurement></xrdMeasurements>'
    )
    with pytest.raises(ValueError, match='Incomplete XRD axis Chi'):
        read_panalytical_xrdml(str(path))


def test_mixed_scan_metadata_is_not_taken_from_one_xrdml_scan(tmp_path):
    scans = ''.join(
        f'<scan scanAxis="{axis}" mode="{mode}"><dataPoints><counts>1 2 3</counts>'
        '<positions axis="2Theta" unit="deg"><startPosition>10</startPosition>'
        '<endPosition>30</endPosition></positions>'
        f'<positions axis="Omega" unit="deg"><commonPosition>{index}</commonPosition>'
        '</positions></dataPoints></scan>'
        for index, (axis, mode) in enumerate(
            [('2Theta', 'Continuous'), ('Omega-2Theta', 'Step')]
        )
    )
    path = tmp_path / 'mixed.xrdml'
    path.write_text(
        '<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.1">'
        f'<xrdMeasurement>{scans}</xrdMeasurement></xrdMeasurements>'
    )
    result = read_panalytical_xrdml(str(path))
    assert result['metadata']['scan_axis'] is None
    assert result['metadata']['scan_mode'] is None
    assert result['intensity'].shape == (2, 3)


def test_rasx_uses_each_profiles_axis_units():
    reader = RASXfile.__new__(RASXfile)
    angles = np.array([36.544, 37.144, 42.544])
    reader.data = [
        np.column_stack((a, [1, 2, 3], [1, 1, 1])) for a in (angles, np.deg2rad(angles))
    ]
    reader.meta = [
        dict(
            ScanInformation=dict(AxisName='TwoTheta', PositionUnit=unit),
            Axes={
                'Omega': SimpleNamespace(Name='Omega', Position=position, Unit='deg')
            },
        )
        for unit, position in [('deg', 15), ('rad', 16)]
    ]
    scanned = reader.get_scan_data()
    assert detect_scan_type(scanned) == 'rsm'
    result = modify_scan_data(scanned, 'rsm')
    np.testing.assert_allclose(result['2Theta'].to('degree').magnitude, angles)
    assert result['2Theta'].shape == angles.shape
    np.testing.assert_array_equal(result['Omega'].magnitude, [15, 16])
