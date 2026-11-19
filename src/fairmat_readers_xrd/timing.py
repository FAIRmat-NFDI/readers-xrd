"""Optional nominal counting times, independent of detector channel count."""

import numpy as np
import pint

ureg = pint.get_application_registry()


def counting_time(value, unit, logger=None):
    """Return positive finite times in seconds, or omit unusable optional metadata."""
    if value is None or unit is None:
        return None
    try:
        values = value.split() if isinstance(value, str) else value
        result = (np.atleast_1d(np.asarray(values, dtype=float)) * ureg(unit)).to('s')
        if (
            not result.size
            or not np.all(np.isfinite(result.magnitude))
            or np.any(result.magnitude <= 0)
        ):
            raise ValueError('Counting time must be positive and finite.')
        return result
    except (ValueError, TypeError, pint.errors.PintError) as error:
        if logger is not None:
            logger.warning(f'Ignoring invalid nominal counting time: {error}')
        return None


def rasx_counting_time(info, logger=None):
    """STEP Speed is a duration; CONTINUOUS Step/Speed is a nominal point time."""
    try:
        mode = str(info.get('Mode', '')).upper()
        if mode == 'STEP':
            return counting_time(info.get('Speed'), info.get('SpeedUnit'), logger)
        if mode == 'CONTINUOUS':
            step = float(info['Step']) * ureg(info['PositionUnit'])
            speed = float(info['Speed']) * ureg(info['SpeedUnit'])
            if not np.isfinite(speed.magnitude) or speed.magnitude == 0:
                raise ValueError('Scan speed must be finite and nonzero.')
            seconds = abs(step / speed).to('s').magnitude
            return counting_time(seconds, 's', logger)
    except (KeyError, TypeError, ValueError, pint.errors.PintError) as error:
        if logger is not None:
            logger.warning(f'Ignoring invalid nominal RASX point time: {error}')
    return None


def complete_scan_times(values, logger=None):
    """Keep scan order; omit incomplete metadata instead of inventing missing times."""
    if not values or any(value is None for value in values):
        if any(value is not None for value in values) and logger is not None:
            logger.warning(
                'Incomplete scan counting times; omitting counting-time metadata.'
            )
        return None
    return values
