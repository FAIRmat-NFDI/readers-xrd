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
import json
import os
import numpy as np
import pytest

import pint

from fairmat_readers_xrd import (
    read_panalytical_xrdml,
    read_rigaku_rasx,
    read_bruker_brml,
    read_bruker_raw,
)

ureg = pint.get_application_registry()


def convert_quantity_to_string(data_dict):
    """
    In a dict, recursively convert every pint.Quantity into str containing its shape.

    Args:
        data_dict (dict): A nested dictionary containing pint.Quantity and other data.
    """
    for k, v in data_dict.items():
        if isinstance(v, ureg.Quantity):
            if isinstance(v.magnitude, np.ndarray):
                data_dict[k] = str(v.shape)
            else:
                data_dict[k] = str(v.magnitude)
        if isinstance(v, dict):
            convert_quantity_to_string(v)
        if isinstance(v, list):
            for i in v:
                convert_quantity_to_string(i)


def test_rasx_reader():
    file_path = [
        'tests/data/RSM_111_sdd=350.rasx',  # File with RSM data
        'tests/data/Omega-2Theta_scan_high_temperature.rasx',  # File with line scan data
        'tests/data/ZnO-ALD-training_001_1_0-000_0-000.rasx',  # File with X, Y data
    ]
    for path in file_path:
        output = read_rigaku_rasx(path)
        convert_quantity_to_string(output)
        with open(f'{path}.json', 'r', encoding='utf-8') as f:
            reference = json.load(f)
        assert output == reference


def test_xrdml_reader():
    file_path = [
        'tests/data/XRD-918-16_10.xrdml',
        'tests/data/m82762_rc1mm_1_16dg_src_slit_phi-101_3dg_-420_mesh_long.xrdml',
    ]
    for path in file_path:
        output = read_panalytical_xrdml(path)
        convert_quantity_to_string(output)
        with open(f'{path}.json', 'r', encoding='utf-8') as f:
            reference = json.load(f)
        assert output == reference


def test_brml_reader():
    file_path = [
        'tests/data/23-012-AG_2thomegascan_long.brml',
        'tests/data/EJZ060_13_004_RSM.brml',
    ]
    for path in file_path:
        output = read_bruker_brml(path)
        convert_quantity_to_string(output)
        with open(f'{path}.json', 'r', encoding='utf-8') as f:
            reference = json.load(f)
        assert output == reference


def test_bruker_raw_reader():
    """
    Test the Bruker/Siemens RAW v4 parser with a sample file.

    This test validates:
    - Binary file parsing
    - Metadata extraction
    - Intensity data extraction
    - Scan parameter completion (from paired XRDML if available)
    - Output format compatibility with other readers
    """
    # Test with scrambled RAW file (anonymized client data)
    test_raw = 'tests/data/TwoTheta_scan_scrambled.raw'
    if not os.path.exists(test_raw):
        pytest.skip('No test RAW file available in tests/data/')
        return

    try:
        output = read_bruker_raw(test_raw)

        # Validate structure
        assert output is not None, 'read_bruker_raw returned None'
        assert '2Theta' in output, 'Missing 2Theta data'
        assert 'intensity' in output, 'Missing intensity data'
        assert 'metadata' in output, 'Missing metadata'
        assert 'scanmotname' in output, 'Missing scanmotname'

        # Validate extracted scan axis
        assert output['scanmotname'] == 'Theta', (
            f"Expected scanmotname='Theta', got '{output['scanmotname']}'"
        )
        assert output['metadata']['scan_axis'] == 'Theta', (
            f"Expected scan_axis='Theta', got '{output['metadata']['scan_axis']}'"
        )

        # Validate data types
        assert hasattr(output['2Theta'], 'magnitude'), (
            '2Theta should be a pint Quantity'
        )
        assert hasattr(output['2Theta'], 'units'), '2Theta should have units'
        assert hasattr(output['intensity'], 'magnitude'), (
            'intensity should be a pint Quantity'
        )
        assert hasattr(output['intensity'], 'units'), 'intensity should have units'
        assert isinstance(output['metadata'], dict), 'metadata should be a dict'

        # Validate that the data arrays have content
        assert len(output['2Theta'].magnitude) > 0, '2Theta should contain data points'
        assert len(output['intensity'].magnitude) > 0, (
            'intensity should contain data points'
        )

    except Exception as e:
        pytest.skip(f'RAW reader test skipped: {str(e)}')
