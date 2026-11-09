"""Pin Kα references and their use by the RAW metadata lookup."""

import pytest

from fairmat_readers_xrd.bruker_raw_parser import XRAY_WAVELENGTHS
from fairmat_readers_xrd.wavelengths import K_ALPHA_WAVELENGTHS


@pytest.mark.parametrize(
    ('material', 'expected'),
    [
        ('Cr', (2.2897260, 2.2936510)),
        ('Fe', (1.9360410, 1.9399730)),
        ('Co', (1.7889960, 1.7928350)),
        ('Cu', (1.54059290, 1.54442740)),
        ('Ga', (1.3401270, 1.3440260)),
        ('Mo', (0.70931715, 0.713607)),
        ('Ag', (0.55942178, 0.5638131)),
        ('In', (0.5121251, 0.5165572)),
    ],
)
def test_kalpha_matches_iucr_reference_wavelengths(material, expected):
    # International Tables for Crystallography (2006), Vol. C, Table 4.2.2.1.
    assert K_ALPHA_WAVELENGTHS[material] == expected


@pytest.mark.parametrize('material', ['Cu', 'Mo', 'Co', 'Fe', 'Cr', 'Ag'])
def test_raw_lookup_uses_shared_lines_and_consistent_weighted_mean(material):
    reference = K_ALPHA_WAVELENGTHS[material]
    raw = XRAY_WAVELENGTHS[material]
    assert (raw['K_alpha1'], raw['K_alpha2']) == reference
    ratio = raw['K_alpha2_K_alpha1_ratio']
    assert raw['K_alpha_avg'] == pytest.approx(
        (reference[0] + ratio * reference[1]) / (1 + ratio), rel=0, abs=1e-15
    )
