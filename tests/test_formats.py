# SPDX-FileCopyrightText: The NOMAD Authors
# SPDX-License-Identifier: Apache-2.0

"""Focused reader data contracts."""

import struct
from pathlib import Path
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
from fairmat_readers_xrd.metadata_extraction import raw_metadata
from fairmat_readers_xrd.wavelengths import K_ALPHA_WAVELENGTHS


def raw4_structure_raw4(points=5, record_size=4, extra_size=0, mode='Locked Coupled'):
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


@pytest.mark.parametrize('record_size', [4, 12])
@pytest.mark.parametrize('extra_size', [37])
def test_raw4_uses_declared_sizes(tmp_path, record_size, extra_size):
    path = tmp_path / 'scan.raw'
    path.write_bytes(
        raw4_structure_raw4(record_size=record_size, extra_size=extra_size)
    )
    result = read_bruker_raw(str(path))
    np.testing.assert_array_equal(result['intensity'].magnitude, [1, 2, 3, 4, 5])
    np.testing.assert_allclose(result['2Theta'].magnitude, 10 + np.arange(5) * 0.01)
    np.testing.assert_allclose(result['countTime'].to('s').magnitude, [0.25])
    assert result['metadata']['scan_axis'] == '2Theta'
    assert result['metadata']['source']['kAlpha1'].to(
        'angstrom'
    ).magnitude == pytest.approx(1.541)
    assert result['metadata']['source']['kAlpha2'] is None
    assert result['metadata']['source']['ratioKAlpha2KAlpha1'] == 0
    assert result['metadata']['source']['voltage'].to('V').magnitude == pytest.approx(
        40000
    )


def test_raw4_point_count_is_not_truncated_to_16_bits(tmp_path):
    path = tmp_path / 'large.raw'
    points = 65537
    path.write_bytes(raw4_structure_raw4(points=points))
    result = BrukerRAW4Parser(str(path)).parse()
    assert len(result['intensities']) == points
    assert result['intensities'][-1] == points


@pytest.mark.parametrize('bad', ['mode', 'truncated', 'block_length'])
def test_raw4_rejects_unsupported_or_incomplete_data(tmp_path, bad):
    content = raw4_structure_raw4(
        mode='Rocking Curve' if bad == 'mode' else 'Locked Coupled'
    )
    if bad == 'truncated':
        content = content[:-1]
    if bad == 'block_length':
        content = bytearray(content)
        struct.pack_into('<I', content, 65, 100000)
    path = tmp_path / 'bad.raw'
    path.write_bytes(content)
    with pytest.raises(ValueError):
        read_bruker_raw(str(path))


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


DATA = Path(__file__).parent / 'data'


def test_xrdml_preserves_explicit_optics_and_avoids_product_inference():
    details = read_panalytical_xrdml(str(DATA / 'm54313_om2th_10.xrdml'))['metadata'][
        'details'
    ]
    mono = details['monochromators'][0]
    assert mono['values']['chemical_formula'] == 'Ge'
    np.testing.assert_array_equal(mono['values']['reflection'], [2, 2, 0])
    assert len(mono['sources']['reflection']) == len(('h', 'k', 'l'))
    assert 'order_no' not in mono['values']
    assert 'wavelength' not in mono['values']
    collimator = details['collimators'][0]
    assert collimator['values']['soller_angle'].to('rad').magnitude == pytest.approx(
        0.04
    )
    assert any(
        a['values']['size'].to('mm').magnitude == 1 for a in details['apertures']
    )
    assert all(
        'beamAttenuationFactors' not in x['values'] for x in details['attenuators']
    )
    assert 'acquisition_mode' not in details['detector'][0]['values']
    assert details['instrument_notes']


def test_raw4_comment_and_creator_version_remain_separate():
    details = raw_metadata(
        {
            'user': 'Operator',
            'comment': 'Comment',
            'creator': 'Recorder',
            'creator_version': '0',
        }
    )
    assert details['notes'][0]['values']['text'] == 'Comment'
    assert details['program'][0]['values']['version'] == '0'
    assert raw_metadata({}) == {'users': [], 'notes': [], 'program': []}


@pytest.mark.parametrize(
    ('material', 'expected'),
    [('Cu', (1.5405929, 1.5444274)), ('Mo', (0.70931715, 0.713607))],
)
def test_kalpha_matches_iucr_reference_wavelengths(material, expected):
    assert K_ALPHA_WAVELENGTHS[material] == expected
