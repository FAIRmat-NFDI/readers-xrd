"""Optional metadata retains source paths, units and scan identities."""

from pathlib import Path
from types import SimpleNamespace
from xml.etree import ElementTree as ET

import numpy as np
import pytest
import structlog
from structlog.testing import capture_logs

from fairmat_readers_xrd import (
    read_bruker_brml,
    read_bruker_raw,
    read_panalytical_xrdml,
    read_rigaku_rasx,
)
from fairmat_readers_xrd.metadata import rasx_metadata, raw_metadata, xrdml_metadata

DATA = Path(__file__).parent / 'data'
NS = {'xrdml': 'http://www.xrdml.com/XRDMeasurement/2.1'}


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


def test_many_headers_create_one_person_and_one_program():
    data = read_panalytical_xrdml(
        str(DATA / 'm82762_rc1mm_1_16dg_src_slit_phi-101_3dg_-420_mesh_long.xrdml')
    )
    details = data['metadata']['details']
    assert len(details['users']) == 1
    assert len(details['program']) == 1
    assert len(details['users'][0]['sources']['name']) == data['intensity'].shape[0]
    assert len(details['program'][0]['sources']['name']) == data['intensity'].shape[0]


@pytest.mark.parametrize('missing', [False, True])
def test_conflicting_or_missing_program_versions_are_not_picked_from_first_scan(
    missing,
):
    other = (
        '<applicationSoftware version="2">Recorder</applicationSoftware>'
        if not missing
        else '<applicationSoftware>Recorder</applicationSoftware>'
    )
    root = ET.fromstring(
        f'<xrdMeasurements xmlns="{NS["xrdml"]}"><xrdMeasurement>'
        '<scan><header><source><applicationSoftware version="1">'
        'Recorder</applicationSoftware>'
        '</source></header></scan><scan><header><source>'
        + other
        + '</source></header></scan></xrdMeasurement></xrdMeasurements>'
    )
    with capture_logs() as logs:
        details = xrdml_metadata(root, NS, structlog.get_logger())
    assert details['program'][0]['values'] == {'name': 'Recorder'}
    assert any('Mixed scan metadata' in row['event'] for row in logs)


def test_incomplete_program_name_omits_common_program():
    root = ET.fromstring(
        f'<xrdMeasurements xmlns="{NS["xrdml"]}"><xrdMeasurement>'
        '<scan><header><source><applicationSoftware>Recorder</applicationSoftware>'
        '</source></header></scan><scan/></xrdMeasurement></xrdMeasurements>'
    )
    assert not xrdml_metadata(root, NS)['program']


def test_invalid_optional_aperture_unit_does_not_guess_or_abort():
    root = ET.fromstring(
        f'<xrdMeasurements xmlns="{NS["xrdml"]}"><xrdMeasurement>'
        '<incidentBeamPath><mask id="a" name="Mask"><width unit="seconds">1</width>'
        '</mask></incidentBeamPath></xrdMeasurement></xrdMeasurements>'
    )
    with capture_logs() as logs:
        details = xrdml_metadata(root, NS, structlog.get_logger())
    assert details['apertures'][0]['values']['name'] == 'Mask'
    assert 'size' not in details['apertures'][0]['values']
    assert any('invalid optional metadata' in row['event'] for row in logs)


def test_rasx_known_positions_and_general_information():
    details = read_rigaku_rasx(str(DATA / 'Omega-2Theta_scan_high_temperature.rasx'))[
        'metadata'
    ]['details']
    assert details['users'][0]['values']['name']
    assert details['detector'][0]['values']['description']
    assert 'Temp' not in {p['values']['name'] for p in details['positioners']}
    assert any(p['values']['value'].magnitude[0] == 0 for p in details['positioners'])
    assert all(p['values']['value'].shape == (1,) for p in details['positioners'])
    assert all('target_value' not in p['values'] for p in details['positioners'])
    assert all(
        'MesurementConditions0.xml' in p for p in details['users'][0]['sources']['name']
    )
    assert details['instrument_notes']


def test_rasx_mixed_detectors_are_not_reported_as_one_common_detector():
    reader = SimpleNamespace(
        data=[None, None],
        metadata_paths=['first.xml', 'second.xml'],
        meta=[
            {'HardwareConfig': {'Detector': 'First'}},
            {'HardwareConfig': {'Detector': 'Second'}},
        ],
    )
    with capture_logs() as logs:
        details = rasx_metadata(reader, structlog.get_logger())
    assert not details['detector']
    assert any('Mixed scan metadata' in row['event'] for row in logs)


def test_brml_uses_both_tracks_and_preserves_fixed_motor_arrays():
    details = read_bruker_brml(str(DATA / '23-012-AG_2thomegascan_long.brml'))[
        'metadata'
    ]['details']
    assert details['detector'][0]['values']['description'] == 'Scinti'
    assert details['mirrors'] and details['apertures']
    assert details['source_fabrication'][0]['values']['serial_number'] == '541635'
    assert any(
        'PrimaryTracks' in source
        for a in details['apertures']
        for source in a['sources']['size']
    )
    assert any(
        'SecondaryTracks' in source
        for a in details['apertures']
        for source in a['sources']['size']
    )
    assert all(p['values']['value'].shape == (1,) for p in details['positioners'])
    assert not details.get('monochromators')


def test_raw4_creator_and_user_have_tag_provenance():
    details = read_bruker_raw(str(DATA / 'TwoTheta_scan_scrambled.raw'))['metadata'][
        'details'
    ]
    assert details['users'][0]['sources']['name'] == ['RAW4/USER']
    assert details['program'][0]['sources']['name'] == ['RAW4/CREATOR']
    assert 'version' not in details['program'][0]['values']


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
