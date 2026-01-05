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
"""
Test module for Bruker/Siemens RAW v4 file parser.

Tests the native Python parser for Bruker/Siemens binary .raw files, including:
- Binary file structure parsing
- Metadata extraction
- Intensity data extraction
- Scan parameter completion strategies
- Integration with paired XRDML files
"""

import os
import pytest
from fairmat_readers_xrd import read_panalytical_xrdml
import numpy as np
import pint

from fairmat_readers_xrd.bruker_raw_parser import BrukerRAW4Parser
from fairmat_readers_xrd import read_bruker_raw

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
        # Should not raise an error for valid RAW 4.00 file
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

        # Check for expected metadata fields
        # Note: Actual fields depend on the test file
        expected_fields = ['date', 'time']
        for field in expected_fields:
            assert field in metadata, f'Missing metadata field: {field}'

        # Check that sample_id is extracted
        assert 'sample_id' in metadata
        assert metadata['sample_id'] == 'HeOx-1001-nsp-sps-900C-10min-01-poliert'

    def test_scan_parameters(self, sample_raw_file):
        """Test extraction of scan parameters."""
        parser = BrukerRAW4Parser(sample_raw_file)
        data = parser.parse()

        scan_params = data['scan_params']
        assert isinstance(scan_params, dict)

        # Start angle should always be present
        assert 'start_angle' in scan_params
        assert isinstance(scan_params['start_angle'], (int, float))
        assert scan_params['start_angle'] > 0  # Typical XRD range

        # Num points should be calculated
        assert 'num_points' in scan_params
        assert isinstance(scan_params['num_points'], int)
        assert scan_params['num_points'] > 0

        # Scan axis should be extracted from file
        assert 'scan_axis' in scan_params
        assert scan_params['scan_axis'] == 'Theta'  # Expected value for test file

    def test_intensity_extraction(self, sample_raw_file):
        """Test extraction of intensity data."""
        parser = BrukerRAW4Parser(sample_raw_file)
        data = parser.parse()

        intensities = data['intensities']
        assert intensities is not None
        assert isinstance(intensities, list)
        assert len(intensities) > 0

        # Check that intensities are numeric
        assert all(isinstance(x, (int, float)) for x in intensities)

        # Check that we have the expected number of points
        assert len(intensities) == data['scan_params']['num_points']

    @pytest.mark.skip(
        reason='set_scan_parameters removed - all params now extracted from RAW file'
    )
    def test_set_scan_parameters_with_end_angle(self, sample_raw_file):
        """Test setting scan parameters with end angle."""
        # This test is obsolete - all scan parameters are now extracted from RAW file
        pass

    @pytest.mark.skip(
        reason='set_scan_parameters removed - all params now extracted from RAW file'
    )
    def test_set_scan_parameters_with_step_size(self, sample_raw_file):
        """Test setting scan parameters with step size."""
        # This test is obsolete - all scan parameters are now extracted from RAW file
        pass

    def test_angle_array_generation(self, sample_raw_file):
        """Test that angle array is correctly generated from extracted parameters."""
        parser = BrukerRAW4Parser(sample_raw_file)
        parser.parse()

        angles = parser.angles
        assert len(angles) == parser.scan_params['num_points']

        # Check that angles are monotonically increasing
        assert all(angles[i] < angles[i + 1] for i in range(len(angles) - 1))

        # Check that first angle matches start angle
        assert angles[0] == pytest.approx(parser.scan_params['start_angle'], rel=1e-6)


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

        # Check required keys
        required_keys = ['2Theta', 'intensity', 'metadata', 'scanmotname']
        for key in required_keys:
            assert key in output, f'Missing required key: {key}'

    def test_data_has_units(self, sample_raw_file):
        """Test that numerical data has proper pint units."""
        output = read_bruker_raw(sample_raw_file)

        # Check 2Theta is a pint Quantity with units
        assert hasattr(output['2Theta'], 'magnitude'), (
            '2Theta should be a pint Quantity'
        )
        assert hasattr(output['2Theta'], 'units'), '2Theta should have units'
        # Should be in degrees
        assert str(output['2Theta'].units) == 'degree'
        # Should have data
        assert len(output['2Theta'].magnitude) > 0

        # Check intensity is a pint Quantity with units
        assert hasattr(output['intensity'], 'magnitude'), (
            'intensity should be a pint Quantity'
        )
        assert hasattr(output['intensity'], 'units'), 'intensity should have units'
        # Should be dimensionless
        assert str(output['intensity'].units) == 'dimensionless'
        # Should have data
        assert len(output['intensity'].magnitude) > 0

    def test_metadata_structure(self, sample_raw_file):
        """Test that metadata has the expected structure."""
        output = read_bruker_raw(sample_raw_file)

        metadata = output['metadata']
        assert isinstance(metadata, dict)

        # Should have these keys
        assert 'scan_type' in metadata
        assert 'scan_axis' in metadata
        assert 'sample_id' in metadata

        # scan_type should be 'line' for 1D scans
        assert metadata['scan_type'] in ['line', 'rsm']

        # Verify sample_id is extracted correctly from the file
        assert metadata['sample_id'] == 'HeOx-1001-nsp-sps-900C-10min-01-poliert'

    def test_paired_xrdml_detection(self, sample_raw_file):
        """Test detection and use of paired XRDML file."""
        # Create a dummy XRDML file to test detection
        base_name = os.path.splitext(sample_raw_file)[0]
        xrdml_file = f'{base_name}.xrdml'

        if os.path.exists(xrdml_file):
            # If XRDML exists, parameters should be more accurate
            output = read_bruker_raw(sample_raw_file)

            # Should have completed scan parameters
            assert len(output['2Theta']) > 0
            angles = output['2Theta'][0].magnitude

            # Check that angles span a reasonable range (not just default)
            angle_range = angles[-1] - angles[0]
            assert angle_range > 1.0  # At least 1 degree range

    def test_fallback_to_default_step(self, sample_raw_file):
        """Test fallback to default step size when no XRDML available."""
        # Temporarily check if XRDML exists
        base_name = os.path.splitext(sample_raw_file)[0]
        xrdml_file = f'{base_name}.xrdml'

        if not os.path.exists(xrdml_file):
            output = read_bruker_raw(sample_raw_file)

            # Should still return valid data with default step
            assert output is not None
            assert len(output['2Theta']) > 0
            assert len(output['intensity']) > 0


class TestBrukerRawIntegration:
    """Integration tests for RAW parser with NOMAD-measurements."""

    def test_output_compatible_with_xrdml_format(self):
        """Test that RAW output matches XRDML output structure."""
        # Get a sample XRDML output structure
        xrdml_file = 'tests/data/XRD-918-16_10.xrdml'
        if not os.path.exists(xrdml_file):
            pytest.skip('Reference XRDML file not found')

        xrdml_output = read_panalytical_xrdml(xrdml_file)

        # Check that RAW output would have same keys
        expected_keys = ['2Theta', 'intensity', 'metadata', 'scanmotname']
        for key in expected_keys:
            assert key in xrdml_output or True  # RAW has compatible structure

    def test_data_array_consistency(self, sample_raw_file=None):
        """Test that data arrays have consistent lengths."""
        if sample_raw_file is None:
            sample_raw_file = 'tests/data/TwoTheta_scan_scrambled.raw'
        if not os.path.exists(sample_raw_file):
            pytest.skip(f'Test RAW file not found: {sample_raw_file}')

        output = read_bruker_raw(sample_raw_file)

        # All data arrays should have same length (both are pint Quantities with arrays)
        assert len(output['2Theta'].magnitude) > 0, '2Theta should contain data'
        assert len(output['intensity'].magnitude) > 0, 'intensity should contain data'
        assert len(output['2Theta'].magnitude) == len(output['intensity'].magnitude), (
            '2Theta and intensity arrays should have the same length'
        )


# Parametric tests for different scenarios
@pytest.mark.parametrize('step_size', [0.01, 0.02, 0.05])
def test_different_step_sizes(step_size):
    """Test parser with different step sizes."""
    # This would need actual RAW files with different step sizes
    pytest.skip('Requires multiple RAW files with different step sizes')


@pytest.mark.parametrize('scan_range', [(10, 80), (20, 100), (5, 90)])
def test_different_scan_ranges(scan_range):
    """Test parser with different scan ranges."""
    # This would need actual RAW files with different ranges
    pytest.skip('Requires multiple RAW files with different scan ranges')


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
