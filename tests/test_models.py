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
import glob

import numpy as np
import pint
import pytest
from pydantic import ValidationError

from fairmat_readers_xrd import read_file, readers
from fairmat_readers_xrd.models import XRDOutput

ureg = pint.get_application_registry()

EXAMPLE_FILES = sorted(
    path for path in glob.glob('tests/data/*') if not path.endswith('.json')
)


def _minimal_output():
    return {
        '2Theta': np.linspace(10, 80, 5) * ureg.degree,
        'intensity': np.ones(5) * ureg.dimensionless,
    }


@pytest.mark.parametrize('path', EXAMPLE_FILES)
def test_example_files_fit_the_model(path):
    output = read_file(path)
    XRDOutput.model_validate(output)


def test_minimal_output_is_valid():
    XRDOutput.model_validate(_minimal_output())


@pytest.mark.parametrize('missing', ['2Theta', 'intensity'])
def test_missing_required_field_raises(missing):
    output = _minimal_output()
    del output[missing]
    with pytest.raises(ValidationError, match=missing):
        XRDOutput.model_validate(output)


def test_unknown_field_raises():
    output = _minimal_output()
    output['2theta'] = output['2Theta']
    with pytest.raises(ValidationError, match='2theta'):
        XRDOutput.model_validate(output)


def test_intensity_and_two_theta_lengths_must_match():
    output = _minimal_output()
    output['intensity'] = np.ones(4) * ureg.dimensionless
    with pytest.raises(ValidationError, match='4 points along 2Theta'):
        XRDOutput.model_validate(output)


@pytest.mark.parametrize(
    'reader, path',
    [
        (readers.read_rigaku_rasx, 'tests/data/TwoTheta_scan_powder.rasx'),
        (readers.read_panalytical_xrdml, 'tests/data/XRD-918-16_10.xrdml'),
        (readers.read_bruker_brml, 'tests/data/23-012-AG_2thomegascan_long.brml'),
        (readers.read_bruker_raw, 'tests/data/TwoTheta_scan_scrambled.raw'),
    ],
)
def test_readers_validate_their_output(reader, path, monkeypatch):
    original = readers.modify_scan_data

    def drop_two_theta(scan_data, scan_type):
        modified = original(scan_data, scan_type)
        modified.pop('2Theta', None)
        return modified

    monkeypatch.setattr(readers, 'modify_scan_data', drop_two_theta)
    with pytest.raises(ValidationError, match='2Theta'):
        reader(path)


def test_unknown_nested_field_raises():
    output = _minimal_output()
    output['metadata'] = {'source': {'kalpha1': 1.5}}
    with pytest.raises(ValidationError, match='kalpha1'):
        XRDOutput.model_validate(output)


@pytest.mark.parametrize(
    'field', ['kAlpha1', 'kAlpha2', 'kBeta', 'ratioKAlpha2KAlpha1']
)
def test_source_values_must_be_quantities(field):
    output = _minimal_output()
    output['metadata'] = {'source': {field: 1.5}}
    with pytest.raises(ValidationError, match=field):
        XRDOutput.model_validate(output)
