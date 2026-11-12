"""Nominal time semantics, scan order, and optional timing metadata."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pint
import pytest

from fairmat_readers_xrd import readers
from fairmat_readers_xrd.timing import counting_time, rasx_counting_time

ureg = pint.get_application_registry()
DATA = Path(__file__).parent / 'data'


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


@pytest.mark.parametrize(
    'value', [None, '', 'invalid', '1 invalid', -1, 0, float('nan'), float('inf')]
)
def test_invalid_optional_times_are_omitted(value):
    assert counting_time(value, 's') is None


def test_incompatible_rasx_time_units_are_omitted():
    assert rasx_counting_time(dict(Mode='STEP', Speed=1, SpeedUnit='deg/min')) is None
    assert (
        rasx_counting_time(
            dict(
                Mode='CONTINUOUS',
                Step=0.1,
                PositionUnit='deg',
                Speed=0,
                SpeedUnit='deg/min',
            )
        )
        is None
    )


def test_rasx_uses_time_of_each_scan_without_channel_expansion(monkeypatch):
    infos = [
        dict(Mode='STEP', AxisName='TwoTheta', Speed=s, SpeedUnit='sec') for s in (1, 2)
    ]
    scan_data = dict(
        intensity=[
            np.array([1, 2, 3]) * ureg.dimensionless,
            np.array([3, 2, 1]) * ureg.dimensionless,
        ],
        **{'2Theta': [np.array([10, 20, 30]) * ureg.degree] * 2},
        Omega=[np.array([1]) * ureg.degree, np.array([2]) * ureg.degree],
    )
    source = dict(
        TargetName='Cu',
        WavelengthKalpha1=None,
        WavelengthKalpha2=None,
        WavelengthKbeta=None,
        Voltage=None,
        Current=None,
    )
    reader = SimpleNamespace(
        data=[None, None],
        meta=[{'ScanInformation': info} for info in infos],
        metadata_paths=['Data0/first.xml', 'Data0/second.xml'],
        get_scan_info=lambda: infos[0],
        get_scan_data=lambda logger: scan_data,
        get_source_info=lambda: source,
    )
    monkeypatch.setattr(readers, 'RASXfile', lambda *args, **kwargs: reader)
    result = readers.read_rigaku_rasx('unused.rasx')
    assert result['intensity'].shape == (2, 3)
    np.testing.assert_allclose(result['countTime'].to('s').magnitude, [1, 2])
    assert result['countTime'].shape == (2,)


@pytest.mark.parametrize(
    ('filename', 'expected'),
    [('23-012-AG_2thomegascan_long.brml', 1), ('EJZ060_13_004_RSM.brml', 9)],
)
def test_brml_effective_times_are_seconds_without_channel_multiplication(
    filename, expected
):
    result = readers.read_bruker_brml(str(DATA / filename))
    np.testing.assert_allclose(result['countTime'].to('s').magnitude, [expected])
    assert result['metadata']['start_time']
    assert result['metadata']['end_time']


def test_rasx_step_fixture_and_measurement_boundaries():
    result = readers.read_rigaku_rasx(str(DATA / 'RSM_111_sdd=350.rasx'))
    np.testing.assert_allclose(result['countTime'].to('s').magnitude, [1.2])
    assert result['metadata']['start_time'] == '2021-05-05T14:27:55Z'
    assert result['metadata']['end_time'] == '2021-05-05T14:35:56Z'


@pytest.mark.parametrize(
    'time_xml', ['', '<commonCountingTime unit="seconds">invalid</commonCountingTime>']
)
def test_xrdml_missing_or_invalid_time_does_not_block_angles(tmp_path, time_xml):
    path = tmp_path / 'scan.xrdml'
    path.write_text(f"""<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/2.1">
      <xrdMeasurement><scan scanAxis="2Theta"><dataPoints>
      {time_xml}<counts>1 2 3</counts>
      <positions axis="2Theta" unit="deg"><startPosition>10</startPosition>
      <endPosition>30</endPosition></positions></dataPoints></scan></xrdMeasurement>
    </xrdMeasurements>""")
    result = readers.read_panalytical_xrdml(str(path))
    assert result['countTime'] is None
    np.testing.assert_allclose(result['2Theta'].magnitude, [10, 20, 30])
