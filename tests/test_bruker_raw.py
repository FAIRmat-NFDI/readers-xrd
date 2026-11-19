# SPDX-FileCopyrightText: The NOMAD Authors
# SPDX-License-Identifier: Apache-2.0

"""Regression checks for native RAW4 output."""

import os

import pint
import pytest

from fairmat_readers_xrd import read_bruker_raw
from fairmat_readers_xrd.bruker_raw_parser import BrukerRAW4Parser

ureg = pint.get_application_registry()


class TestBrukerRAW4Parser:
    """Test suite for BrukerRAW4Parser class."""

    @pytest.fixture
    def sample_raw_file(self):
        """
        Fixture providing path to a sample RAW file.

        Uses scrambled test data (anonymized from client data).
        """
        test_file = 'tests/data/TwoTheta_scan_scrambled.raw'
        if not os.path.exists(test_file):
            pytest.skip(f'Test RAW file not found: {test_file}')
        return test_file

    def test_parser_initialization(self):
        """Test that parser can be initialized with a file path."""
        parser = BrukerRAW4Parser('dummy.raw')
        assert parser.filepath.name == 'dummy.raw'
        assert parser.metadata == {}
        assert parser.scan_params == {}
        assert parser.intensities is None
        assert parser.angles is None

    def test_header_validation(self, sample_raw_file):
        """Test that parser validates RAW 4.00 header correctly."""
        parser = BrukerRAW4Parser(sample_raw_file)
        try:
            data = parser.parse()
            assert data is not None
        except ValueError as e:
            if 'Invalid RAW file header' in str(e):
                pytest.fail('Valid RAW file rejected by parser')

    def test_metadata_extraction(self, sample_raw_file):
        """Test extraction of metadata from RAW file."""
        parser = BrukerRAW4Parser(sample_raw_file)
        data = parser.parse()
        metadata = data['metadata']
        assert isinstance(metadata, dict)
        expected_fields = ['date', 'time']
        for field in expected_fields:
            assert field in metadata, f'Missing metadata field: {field}'
        assert 'sample_id' in metadata
        assert metadata['sample_id'] == 'HeOx-1001-nsp-sps-900C-10min-01-poliert'

    def test_scan_parameters(self, sample_raw_file):
        """Test extraction of scan parameters."""
        parser = BrukerRAW4Parser(sample_raw_file)
        data = parser.parse()
        scan_params = data['scan_params']
        assert isinstance(scan_params, dict)
        assert 'start_angle' in scan_params
        assert isinstance(scan_params['start_angle'], (int, float))
        assert scan_params['start_angle'] > 0
        assert 'num_points' in scan_params
        assert isinstance(scan_params['num_points'], int)
        assert scan_params['num_points'] > 0
        assert 'scan_axis' in scan_params
        assert scan_params['scan_axis'] == '2Theta'

    def test_intensity_extraction(self, sample_raw_file):
        """Test extraction of intensity data."""
        parser = BrukerRAW4Parser(sample_raw_file)
        data = parser.parse()
        intensities = data['intensities']
        assert intensities is not None
        assert isinstance(intensities, list)
        assert len(intensities) > 0
        assert all(isinstance(x, (int, float)) for x in intensities)
        assert len(intensities) == data['scan_params']['num_points']

    def test_angle_array_generation(self, sample_raw_file):
        """Test that angle array is correctly generated from extracted parameters."""
        parser = BrukerRAW4Parser(sample_raw_file)
        parser.parse()
        angles = parser.angles
        assert len(angles) == parser.scan_params['num_points']
        assert all(angles[i] < angles[i + 1] for i in range(len(angles) - 1))
        assert angles[0] == pytest.approx(parser.scan_params['start_angle'], rel=1e-06)

    def test_xy_export_leaves_out_unknown_start_and_step(self, sample_raw_file):
        """
        The XY header must not invent a start angle or step of 0
        (Bug B34 in https://github.com/FAIRmat-NFDI/pynxtools-plugin-audit/blob/main/pynxtools-unification/FINDINGS.md).
        """
        # With a real file, start and step are known and written as before.
        parser = BrukerRAW4Parser(sample_raw_file)
        parser.parse()
        header = parser.to_xy_format().splitlines()
        assert '# Start: 10.000°' in header
        assert any(line.startswith('# Step: ') for line in header)

        # Without them, the export used to write "Start: 0.000°" and
        # "Step: 0.0000°", which looks like real data. Now the lines are left out.
        parser.scan_params = {}
        header = parser.to_xy_format().splitlines()
        assert not any(line.startswith(('# Start:', '# Step:')) for line in header)

    def test_truncated_range_is_rejected(self, sample_raw_file, tmp_path):
        """A file too short to hold the start angle must not crash the parser."""
        # The start angle sits at byte 0x04FB. A file cut before it used to
        # crash with UnboundLocalError, because start_angle = None was set in
        # the wrong except block. Now the angle is just reported as missing.
        with open(sample_raw_file, 'rb') as f:
            truncated = f.read(0x04FB)
        truncated_file = tmp_path / 'truncated.raw'
        truncated_file.write_bytes(truncated)

        parser = BrukerRAW4Parser(str(truncated_file))
        with pytest.raises(ValueError, match='Truncated'):
            parser.parse()


class TestReadBrukerRaw:
    """Test suite for read_bruker_raw() integration function."""

    @pytest.fixture
    def sample_raw_file(self):
        """Fixture providing path to a sample RAW file."""
        test_file = 'tests/data/TwoTheta_scan_scrambled.raw'
        if not os.path.exists(test_file):
            pytest.skip(f'Test RAW file not found: {test_file}')
        return test_file

    def test_read_function_returns_dict(self, sample_raw_file):
        """Test that read_bruker_raw returns a properly structured dictionary."""
        output = read_bruker_raw(sample_raw_file)
        assert output is not None
        assert isinstance(output, dict)
        required_keys = ['2Theta', 'intensity', 'metadata', 'scanmotname']
        for key in required_keys:
            assert key in output, f'Missing required key: {key}'

    def test_data_has_units(self, sample_raw_file):
        """Test that numerical data has proper pint units."""
        output = read_bruker_raw(sample_raw_file)
        assert hasattr(output['2Theta'], 'magnitude'), (
            '2Theta should be a pint Quantity'
        )
        assert hasattr(output['2Theta'], 'units'), '2Theta should have units'
        assert str(output['2Theta'].units) == 'degree'
        assert len(output['2Theta'].magnitude) > 0
        assert hasattr(output['intensity'], 'magnitude'), (
            'intensity should be a pint Quantity'
        )
        assert hasattr(output['intensity'], 'units'), 'intensity should have units'
        assert str(output['intensity'].units) == 'dimensionless'
        assert len(output['intensity'].magnitude) > 0

    def test_metadata_structure(self, sample_raw_file):
        """Test that metadata has the expected structure."""
        output = read_bruker_raw(sample_raw_file)
        metadata = output['metadata']
        assert isinstance(metadata, dict)
        assert 'scan_type' in metadata
        assert 'scan_axis' in metadata
        assert 'sample_id' in metadata
        assert metadata['scan_type'] in ['line', 'rsm']
        assert metadata['sample_id'] == 'HeOx-1001-nsp-sps-900C-10min-01-poliert'

    def test_missing_start_angle_raises(self, sample_raw_file, tmp_path):
        """read_bruker_raw must raise for a file without a 2Theta axis (B34)."""
        # A file cut before the start angle cannot give a 2Theta axis. The
        # function used to swallow the error and return None, so callers got
        # no output and, without a logger, no message. Now it raises.
        with open(sample_raw_file, 'rb') as f:
            truncated = f.read(0x04FB)
        truncated_file = tmp_path / 'truncated.raw'
        truncated_file.write_bytes(truncated)

        with pytest.raises(ValueError, match='Truncated'):
            read_bruker_raw(str(truncated_file))


class TestBrukerRawIntegration:
    """Integration tests for RAW parser with NOMAD-measurements."""

    def test_data_array_consistency(self, sample_raw_file=None):
        """Test that data arrays have consistent lengths."""
        if sample_raw_file is None:
            sample_raw_file = 'tests/data/TwoTheta_scan_scrambled.raw'
        if not os.path.exists(sample_raw_file):
            pytest.skip(f'Test RAW file not found: {sample_raw_file}')
        output = read_bruker_raw(sample_raw_file)
        assert len(output['2Theta'].magnitude) > 0, '2Theta should contain data'
        assert len(output['intensity'].magnitude) > 0, 'intensity should contain data'
        assert len(output['2Theta'].magnitude) == len(output['intensity'].magnitude), (
            '2Theta and intensity arrays should have the same length'
        )


def test_raw_wavelengths_are_quantities():
    """The RAW reader returns wavelengths in angstrom, like the other readers."""
    source = read_bruker_raw('tests/data/TwoTheta_scan_scrambled.raw')['metadata'][
        'source'
    ]
    for key in ['kAlpha1', 'kBeta']:
        assert source[key].units == ureg.angstrom
    assert source['kAlpha1'].magnitude == pytest.approx(1.5405999422073364)
    assert source['kAlpha2'] is None
    assert source['ratioKAlpha2KAlpha1'].magnitude == 0
    assert source['ratioKAlpha2KAlpha1'].units == ureg.dimensionless


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
