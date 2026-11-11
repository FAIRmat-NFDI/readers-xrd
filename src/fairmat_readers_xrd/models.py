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
"""Pydantic models describing the output of every `read_*` function."""

import functools
from collections.abc import Callable
from typing import Any

import pint
from pydantic import BaseModel, ConfigDict, Field, model_validator

Quantity = pint.Quantity


class XRDSource(BaseModel):
    """X-ray source information."""

    model_config = ConfigDict(arbitrary_types_allowed=True, extra='forbid')

    anode_material: str | None = None
    kAlpha1: Quantity | float | None = None
    kAlpha2: Quantity | float | None = None
    kBeta: Quantity | float | None = None
    ratioKAlpha2KAlpha1: Quantity | float | None = None
    voltage: Quantity | None = None
    current: Quantity | None = None


class XRDMetadata(BaseModel):
    """Metadata about the sample, the instrument and the scan."""

    model_config = ConfigDict(arbitrary_types_allowed=True, extra='forbid')

    sample_id: str | None = None
    measurement_type: str | None = None
    sample_mode: str | None = None
    instrument_id: str | None = None
    start_time: str | None = None
    scan_mode: str | None = None
    scan_axis: str | None = None
    scan_type: str | None = None
    source: XRDSource | None = None


class XRDOutput(BaseModel):
    """Output of a `read_*` function. Only 2Theta and intensity are required."""

    model_config = ConfigDict(arbitrary_types_allowed=True, extra='forbid')

    two_theta: Quantity = Field(alias='2Theta')
    intensity: Quantity
    Omega: Quantity | None = None
    Chi: Quantity | None = None
    Phi: Quantity | None = None
    X: Quantity | None = None
    Y: Quantity | None = None
    Z: Quantity | None = None
    counts: Quantity | None = None
    countTime: Quantity | None = None
    beamAttenuationFactors: Quantity | None = None
    scanmotname: str | None = None
    metadata: XRDMetadata | None = None

    @model_validator(mode='after')
    def intensity_matches_two_theta(self):
        if self.intensity.shape[-1] != self.two_theta.shape[-1]:
            raise ValueError(
                f'intensity has {self.intensity.shape[-1]} points along 2Theta, '
                f'but 2Theta has {self.two_theta.shape[-1]}'
            )
        return self


def validate_output(reader: Callable) -> Callable:
    """
    Decorator for `read_*` functions: checks their output against `XRDOutput`.

    The output is returned unchanged.

    Raises:
        pydantic.ValidationError: If the output does not fit the model.
    """

    @functools.wraps(reader)
    def wrapper(*args, **kwargs) -> dict[str, Any]:
        output = reader(*args, **kwargs)
        XRDOutput.model_validate(output)
        return output

    return wrapper
