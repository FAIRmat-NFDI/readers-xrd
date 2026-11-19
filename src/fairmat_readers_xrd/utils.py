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
import collections
from typing import Any

import numpy as np
import pint

ureg = pint.get_application_registry()
MAX_SUPPORTED_INTENSITY_DIMENSIONS = 2


def intensity_quantity(values, unit=None, *, prefix=None):
    """Decode counts or an explicit count rate without normalizing by time."""
    aliases = {
        'counts': 'count',
        'cts': 'count',
        'cps': 'count/s',
        'countspersecond': 'count/s',
        'countpersecond': 'count/s',
        'counts/s': 'count/s',
        'cts/s': 'count/s',
    }
    unit = aliases.get(str(unit).casefold().replace(' ', ''), unit)
    try:
        signal = np.asarray(values, dtype=float) * ureg(unit or 'dimensionless')
        if prefix and prefix.casefold() != 'none':
            signal *= ureg(prefix.lower() + 'count').to('count').magnitude
        if signal.is_compatible_with('dimensionless'):
            return signal.to('dimensionless')
        if signal.is_compatible_with('count/s'):
            return signal.to('count/s')
    except pint.errors.PintError as error:
        raise ValueError(f'Invalid XRD signal unit: {unit!r}.') from error
    raise ValueError(f'Unsupported XRD signal unit: {unit!r}.')


def to_pint_quantity(value: Any = None, unit: str = None) -> Any:
    """
    Attempts to generate a pint quantity.
    In case the value is a string, it is returned as is.
    If the value is a pint quantity, it is converted to the given unit.

    Args:
        value (Any): Value of the quantity.
        unit (str): Unit of the quantity.

    Returns:
        Any: Processed quantity with datatype depending on the value.
    """
    if isinstance(value, str) or value is None:
        return value
    if isinstance(value, ureg.Quantity):
        if unit is None:
            return value
        return value.to(unit)
    return value * ureg(unit)


def are_all_identical(arr_list):
    """
    Check if all the arrays in the list are identical. Also works if the arrays are
    pint.Quantity.

    Args:
        arr_list (list): A list of numpy arrays.

    Returns:
        bool: True if all the arrays are identical, False otherwise.
    """
    first_arr = arr_list[0]
    unit = getattr(first_arr, 'units', None)
    if isinstance(first_arr, ureg.Quantity):
        first_arr = first_arr.magnitude

    for arr in arr_list[1:]:
        arr_magnitude = arr.to(unit).magnitude if unit is not None else arr
        # Permit only floating-point unit-conversion roundoff.
        if np.shape(first_arr) != np.shape(arr_magnitude) or not np.allclose(
            first_arr, arr_magnitude, rtol=8 * np.finfo(float).eps, atol=0
        ):
            return False
    return True


def detect_scan_type(scan_data):
    """
    Based on the shape of data vectors, decide whether the scan_type is `line` (single
    line scan), `multiline` (multiple line scans), or `rsm` (reciprocal space mapping).
    For a 2D scan, if the conditions for `rsm` are not met, it is considered
    a `multiline` scan.

    Args:
        scan_data (dict): The X-ray diffraction data in a Python dictionary. Each key is
            a list of scan data as pint.Quantity arrays.

    Returns:
        str: The type of scan.
    """
    count = len(scan_data['intensity'])
    if not count:
        raise ValueError('XRD scan contains no intensity data.')
    # if intensity data is not a regular 2D array, it is not `rsm`
    for scan_intensity in scan_data['intensity'][1:]:
        if scan_intensity.shape != scan_data['intensity'][0].shape:
            return 'multiline'

    for key in ('2Theta', 'Omega', 'Chi', 'Phi', 'Theta', 'X', 'Y', 'Z'):
        values = scan_data.get(key)
        if values is not None and (
            len(values) != count or any(value is None for value in values)
        ):
            raise ValueError(
                f'Incomplete XRD axis {key}; positions cannot be assigned to scans.'
            )
    if len(scan_data['intensity']) == 1:
        return 'line'

    intensity_shape = (
        len(scan_data['intensity']),
        *scan_data['intensity'][0].shape,
    )
    if len(intensity_shape) > MAX_SUPPORTED_INTENSITY_DIMENSIONS:
        raise AssertionError(
            'Scan type not detected. `intensity.ndim` must be 1 or 2. '
            f'Found: {len(intensity_shape)}'
        )

    if scan_data.get('2Theta') is None or not are_all_identical(scan_data['2Theta']):
        return 'multiline'
    # find axis that updates from one scan to other
    var_axis = []
    for key in ['Omega', 'Chi', 'Phi', 'Theta']:
        if scan_data.get(key) is None:
            continue
        data = scan_data[key]
        if not are_all_identical(data):
            var_axis.append(key)
    # if only one var_axis
    # and dimensions of 2theta, var_axis, and intensity are consistent, it is a rsm
    if len(var_axis) == 1:
        two_theta_shape = (len(scan_data['2Theta']), *scan_data['2Theta'][0].shape)
        var_axis_data = scan_data[var_axis[0]]
        if (
            intensity_shape == two_theta_shape
            and len(var_axis_data) == intensity_shape[0]
            and all(np.all(row == row[0]) for row in var_axis_data)
        ):
            return 'rsm'
    return 'multiline'


def modify_scan_data(scan_data: dict, scan_type: str):
    """
    Modifies the scan data based on the scan type:

    If the scan type is `line`, the data is converted to 1D arrays.

    Intensity and counts retain every scan and point, including constant spectra.
    For RSM, signal arrays are two-dimensional. Shared coordinates and optional
    metadata may be reduced to vectors or a single value, in consistent units.

    If the scan type is `multiline`, the data is converted into a list of 1D arrays.
    Currently not implemented.

    Args:
        scan_data (dict): The X-ray diffraction data in a Python dictionary. Each key is
            a list of scan data as pint.Quantity arrays.
        scan_type (str): The type of scan.

    Returns:
        dict: scan_data containing same keys but modified values.
    """
    output = collections.defaultdict(lambda: None)

    if scan_type not in ['line', 'rsm', 'multiline']:
        raise ValueError(f'Invalid scan type: {scan_type}')

    if scan_type == 'line':
        for key, value in scan_data.items():
            if value is None:
                continue
            data = value[0].magnitude
            if key not in ('intensity', 'counts') and np.all(
                np.diff(data, axis=0) == 0
            ):
                # if elements are identical, pick the first one
                data = np.array([data[0]])
            output[key] = data * value[0].units
        return output

    elif scan_type == 'multiline':
        raise NotImplementedError(f'Scan type {scan_type} is not supported.')

    elif scan_type == 'rsm':
        for key, value in scan_data.items():
            if value is None:
                continue
            data = np.array([item.to(value[0].units).magnitude for item in value])
            if key in ('intensity', 'counts') or (
                key == 'countTime' and data.shape[1] > 1
            ):
                output[key] = data * value[0].units
                continue
            # if it is column vector, make it a row vector
            if data.shape[1] == 1 or (
                key in ('Omega', 'Chi', 'Phi', 'Theta') and np.all(data == data[:, :1])
            ):
                data = data[:, 0]
            # if rows (or elements of a row) are identical, pick the first one
            if are_all_identical(data):
                data = data[0].reshape(-1)
            output[key] = data * value[0].units
        return output
