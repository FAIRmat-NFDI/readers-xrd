"""
Pure Python parser for Siemens/Bruker RAW v4 X-ray diffraction files.

This module provides native Python parsing of the Siemens/Bruker proprietary
binary .raw format (magic header "RAW4.00") without requiring external tools,
Wine, or .NET libraries.

Based on reverse engineering of RAW v4 file structure from example files.
Validated on Bruker DIFFRAC.EVA generated single-axis powder diffraction scans.

Author: Generated for NOMAD
License: MIT
"""

import logging
import struct
from pathlib import Path
from typing import Any

import numpy as np

from fairmat_readers_xrd.wavelengths import K_ALPHA_WAVELENGTHS

# Kα references are shared with NM; retain the existing Kβ values and 2:1 ratio.
XRAY_WAVELENGTHS = {
    material: {
        'K_alpha1': K_ALPHA_WAVELENGTHS[material][0],
        'K_alpha2': K_ALPHA_WAVELENGTHS[material][1],
        'K_beta': k_beta,
        'K_alpha_avg': (
            2 * K_ALPHA_WAVELENGTHS[material][0] + K_ALPHA_WAVELENGTHS[material][1]
        )
        / 3,
        'K_alpha2_K_alpha1_ratio': 0.5,
    }
    for material, k_beta in {
        'Cu': 1.392250,
        'Mo': 0.632305,
        'Co': 1.620830,
        'Fe': 1.757462,
        'Cr': 2.084920,
        'Ag': 0.497082,
    }.items()
}


# RAW4 block identifiers and fixed header sizes (bytes).
_BLOCK_HEADER = 8
_PREAMBLE = 61
_TEXT_BLOCK = 10
_TEXT_HEADER = 36
_SOURCE_BLOCK = 30
_SOURCE_HEADER = 120
_DRIVE_BLOCK = 50
_DRIVE_HEADER = 64
_INTENSITY_BYTES = 4


class BrukerRAW4Parser:
    """Read a single RAW4 powder range using its declared block and record sizes.

    The range/drive layout follows GSAS-II's G2pwd_BrukerRAW reader. Only the
    locally validated Locked Coupled geometry is accepted; other modes need
    independent reference files before their axis semantics can be supported.
    """

    HEADER_MAGIC = b'RAW4.00\x00'
    HEADER_SIZE = 8
    METADATA_TAGS = {
        'USER': 'user',
        'SITE': 'site',
        'SAMPLEID': 'sample_id',
        'COMMENT': 'comment',
        'CREATOR': 'creator',
        'CREATOR_VERSION': 'creator_version',
    }

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)
        self.metadata = {}
        self.scan_params = {}
        self.intensities = None
        self.angles = None
        self.logger = logging.getLogger(__name__)

    @staticmethod
    def _text(data):
        return data.split(b'\x00', 1)[0].decode('latin-1').strip()

    @staticmethod
    def _block(data, offset, end):
        if offset + _BLOCK_HEADER > end:
            raise ValueError('Truncated RAW4 block header.')
        kind, size = struct.unpack_from('<II', data, offset)
        if size < _BLOCK_HEADER or offset + size > end:
            raise ValueError('Invalid RAW4 block length.')
        return kind, size

    def parse(self) -> dict[str, Any]:
        data = self.filepath.read_bytes()
        if not data.startswith(self.HEADER_MAGIC) or len(data) < _PREAMBLE:
            raise ValueError('Invalid or truncated RAW4 header.')
        # RAW4 timestamps carry no timezone; retain them without inventing one.
        self.metadata['date'] = self._text(data[12:24])
        self.metadata['time'] = self._text(data[24:34])
        offset = _PREAMBLE
        while offset + 4 <= len(data):
            kind = struct.unpack_from('<I', data, offset)[0]
            if kind in (0, 160):
                self._read_range(data, offset)
                break
            kind, size = self._block(data, offset, len(data))
            block = data[offset : offset + size]
            if kind == _TEXT_BLOCK and size >= _TEXT_HEADER:
                key = self.METADATA_TAGS.get(self._text(block[12:36]))
                if key:
                    self.metadata[key] = self._text(block[36:])
            elif kind == _SOURCE_BLOCK and size >= _SOURCE_HEADER:
                self._read_source(block)
            offset += size
        else:
            raise ValueError('RAW4 measurement range is missing.')
        return dict(
            metadata=self.metadata,
            scan_params=self.scan_params,
            intensities=self.intensities,
            angles=self.angles,
        )

    def _read_source(self, block):
        fields = (
            ('wavelength_kalpha_avg', 'K_alpha_avg'),
            ('wavelength_kalpha1', 'K_alpha1'),
            ('wavelength_kalpha2', 'K_alpha2'),
            ('wavelength_kbeta', 'K_beta'),
            ('kalpha2_kalpha1_ratio', 'K_alpha2_K_alpha1_ratio'),
        )
        material = self._text(block[116:120])
        self.metadata['anode_material'] = material or None
        fallback = XRAY_WAVELENGTHS.get(material, {})
        for (key, fallback_key), raw_value in zip(
            fields, struct.unpack_from('<5d', block, 72)
        ):
            value = raw_value
            # Zero wavelengths are unreported; a zero intensity ratio is valid.
            if not np.isfinite(value) or value < 0:
                self.logger.warning(f'Ignoring invalid RAW4 source value {key}.')
                value = fallback.get(fallback_key)
            elif value == 0 and key != 'kalpha2_kalpha1_ratio':
                value = None
            self.metadata[key] = value

    def _read_range(self, data, offset):
        if offset + 160 > len(data):
            raise ValueError('Truncated RAW4 range header.')
        mode = self._text(data[offset + 32 : offset + 56])
        if mode != 'Locked Coupled':
            raise ValueError(f'Unsupported RAW4 scan mode: {mode!r}.')
        start, step = struct.unpack_from('<2d', data, offset + 72)
        count = struct.unpack_from('<I', data, offset + 88)[0]
        time_ms = struct.unpack_from('<f', data, offset + 92)[0]
        voltage, current = struct.unpack_from('<2f', data, offset + 100)
        record_size, drive_size = struct.unpack_from('<2I', data, offset + 136)
        if not count or not np.isfinite(start) or not np.isfinite(step) or step == 0:
            raise ValueError('Invalid RAW4 scan coordinates or point count.')
        if record_size < _INTENSITY_BYTES:
            raise ValueError('RAW4 records cannot contain an intensity value.')
        position = offset + 160
        end = position + drive_size
        if end > len(data):
            raise ValueError('Truncated RAW4 drive records.')
        drives = {}
        while position < end:
            kind, size = self._block(data, position, end)
            if kind == _DRIVE_BLOCK and size >= _DRIVE_HEADER:
                name = self._text(data[position + 12 : position + 36])
                drives[name] = struct.unpack_from('<d', data, position + 56)[0]
            position += size
        if '2Theta' in drives and not np.isclose(
            drives['2Theta'], start, atol=1e-6, rtol=0
        ):
            raise ValueError('RAW4 2Theta drive contradicts the range start.')
        data_end = end + count * record_size
        if data_end != len(data):
            raise ValueError(
                'RAW4 data length does not match one complete range; '
                'multiple ranges are not supported.'
            )
        self.intensities = np.ndarray(
            (count,), dtype='<f4', buffer=data, offset=end, strides=(record_size,)
        ).tolist()
        self.angles = (start + np.arange(count) * step).tolist()
        self.scan_params.update(
            scan_axis='2Theta',
            scan_mode=mode,
            start_angle=start,
            step_size=step,
            num_points=count,
            end_angle=self.angles[-1],
            counting_time_ms=time_ms,
        )
        self.metadata.update(voltage=voltage, current=current)

    def to_xy_format(self) -> str:
        """
        Convert parsed data to simple XY format (angle, intensity).

        Returns:
            String in XY format
        """
        if not self.intensities or not self.angles:
            raise ValueError('No data parsed. Call parse() first.')

        lines = ['# Bruker/Siemens RAW v4 converted to XY format']
        lines.append(f'# Sample: {self.metadata.get("sample_id", "Unknown")}')
        lines.append(
            f'# Date: {self.metadata.get("date", "")} {self.metadata.get("time", "")}'
        )
        if 'start_angle' in self.scan_params:
            lines.append(f'# Start: {self.scan_params["start_angle"]:.3f}°')
        if 'step_size' in self.scan_params:
            lines.append(f'# Step: {self.scan_params["step_size"]:.4f}°')
        lines.append(f'# Points: {len(self.intensities)}')
        lines.append('#')
        lines.append('# 2Theta(°)  Intensity(counts)')

        for angle, intensity in zip(self.angles, self.intensities):
            lines.append(f'{angle:.4f}  {intensity:.2f}')

        return '\n'.join(lines)


def read_bruker_raw4(filepath: str) -> dict[str, Any]:
    """
    Read a Bruker/Siemens RAW v4 file.

    Args:
        filepath: Path to .raw file

    Returns:
        Dictionary with 'metadata', 'scan_params', 'intensities', 'angles'

    Example:
        >>> data = read_bruker_raw4('sample.raw')
        >>> print(data['metadata']['sample_id'])
        >>> print(len(data['intensities']))
    """
    parser = BrukerRAW4Parser(filepath)
    return parser.parse()


def convert_raw_to_xy(raw_filepath: str, xy_filepath: str | None = None) -> str:
    """
    Convert Bruker RAW v4 .raw file to simple XY format.

    Args:
        raw_filepath: Path to input .raw file
        xy_filepath: Optional output path. If None, replaces .raw with .xy

    Returns:
        Path to created XY file
    """
    parser = BrukerRAW4Parser(raw_filepath)
    parser.parse()

    xy_content = parser.to_xy_format()

    if xy_filepath is None:
        xy_filepath = str(Path(raw_filepath).with_suffix('.xy'))

    with open(xy_filepath, 'w') as f:
        f.write(xy_content)

    return xy_filepath
