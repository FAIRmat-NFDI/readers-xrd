"""Declared count/rate units and scan-wise acquisition times."""

from pathlib import Path
from types import SimpleNamespace
from xml.etree import ElementTree
from zipfile import ZipFile

import numpy as np
import pytest

from fairmat_readers_xrd import read_bruker_brml, read_panalytical_xrdml
from fairmat_readers_xrd.ikz import RASXfile
from fairmat_readers_xrd.utils import intensity_quantity


def xrdml_scan(unit, values='2 4', times='2 3', omega=1):
    return f"""<scan scanAxis="2Theta"><dataPoints>
      <countingTimes unit="seconds">{times}</countingTimes>
      <intensities unit="{unit}">{values}</intensities>
      <positions axis="2Theta" unit="deg"><startPosition>10</startPosition>
      <endPosition>20</endPosition></positions>
      <positions axis="Omega" unit="deg"><commonPosition>{omega}</commonPosition>
      </positions></dataPoints></scan>"""


def write_xrdml(tmp_path, scans):
    path = tmp_path / 'scan.xrdml'
    path.write_text(
        '<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.1">'
        '<xrdMeasurement>' + ''.join(scans) + '</xrdMeasurement></xrdMeasurements>'
    )
    return path


@pytest.mark.parametrize(
    ('unit', 'target', 'expected'),
    [
        ('counts', 'dimensionless', [2, 4]),
        ('cps', 'count/s', [2, 4]),
        ('counts per second', 'count/s', [2, 4]),
        ('counts/minute', 'count/s', [2 / 60, 4 / 60]),
    ],
)
def test_xrdml_declared_signal_scale_is_decoded(tmp_path, unit, target, expected):
    data = read_panalytical_xrdml(str(write_xrdml(tmp_path, [xrdml_scan(unit)])))
    np.testing.assert_allclose(data['intensity'].to(target).magnitude, expected)
    np.testing.assert_allclose(data['countTime'].to('s').magnitude, [2, 3])


def test_point_times_remain_a_matrix_for_rsm(tmp_path):
    path = write_xrdml(tmp_path, [xrdml_scan('cps'), xrdml_scan('cps', omega=2)])
    data = read_panalytical_xrdml(str(path))
    assert data['metadata']['scan_type'] == 'rsm'
    assert data['countTime'].shape == (2, 2)
    np.testing.assert_array_equal(data['countTime'].magnitude, [[2, 3], [2, 3]])


@pytest.mark.parametrize('unit', ['watts', 'not_a_signal_unit'])
def test_unknown_signal_units_are_not_treated_as_counts(tmp_path, unit):
    with pytest.raises(ValueError, match='XRD signal unit'):
        read_panalytical_xrdml(str(write_xrdml(tmp_path, [xrdml_scan(unit)])))


def test_mixed_count_and_rate_scans_are_not_combined(tmp_path):
    path = write_xrdml(tmp_path, [xrdml_scan('counts'), xrdml_scan('cps', omega=2)])
    with pytest.raises(ValueError, match='Mixed counts and count-rate'):
        read_panalytical_xrdml(str(path))


def test_rasx_uses_the_reported_intensity_unit():
    reader = SimpleNamespace(
        data=[np.array([[10, 2], [20, 4]])],
        meta=[
            {
                'ScanInformation': {'AxisName': 'TwoTheta', 'IntensityUnit': 'cps'},
                'Axes': {},
            }
        ],
    )
    data = RASXfile.get_scan_data(reader)
    np.testing.assert_allclose(data['intensity'][0].to('count/s').magnitude, [2, 4])


def test_brml_recording_unit_is_used(tmp_path):
    source = Path(__file__).parent / 'data/23-012-AG_2thomegascan_long.brml'
    path = tmp_path / 'rates.brml'
    with ZipFile(source) as original, ZipFile(path, 'w') as modified:
        for name in original.namelist():
            content = original.read(name)
            if 'RawData' in name and name.endswith('.xml'):
                root = ElementTree.fromstring(content)
                for recording in root.iter('Recording'):
                    if 'counter' in recording.get('LogicName', '').lower():
                        recording.find('Unit').set('Base', 'CountsPerSecond')
                content = ElementTree.tostring(root)
            modified.writestr(name, content)
    counts = read_bruker_brml(str(source))['intensity']
    rates = read_bruker_brml(str(path))['intensity']
    np.testing.assert_allclose(rates.to('count/s').magnitude, counts.magnitude)


def test_brml_signal_prefix_is_converted_to_the_common_scale():
    signal = intensity_quantity([2, 4], 'CountsPerSecond', prefix='Kilo')
    np.testing.assert_allclose(signal.to('count/s').magnitude, [2000, 4000])
