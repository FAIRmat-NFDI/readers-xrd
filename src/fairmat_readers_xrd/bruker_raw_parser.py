"""
Pure Python parser for Siemens/Bruker RAW v4 X-ray diffraction files.

This module provides native Python parsing of the Siemens/Bruker proprietary binary
.raw format (magic header "RAW4.00") without requiring external tools, Wine, or .NET libraries.

Based on reverse engineering of RAW v4 file structure from example files.
Validated on Bruker DIFFRAC.EVA generated single-axis powder diffraction scans.

Author: Generated for NOMAD
License: MIT
"""

import struct
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
import logging


# X-ray wavelength reference table (in Angstroms)
# Source: International Tables for Crystallography, Volume C
# Reference: https://onlinelibrary.wiley.com/iucr/itc/Cb/ch4o2v0001/sec4o2o2/
XRAY_WAVELENGTHS = {
    'Cu': {
        'K_alpha1': 1.540598,
        'K_alpha2': 1.544426,
        'K_beta': 1.392250,
        'K_alpha_avg': 1.541874,  # Weighted average
        'K_alpha2_K_alpha1_ratio': 0.5,
    },
    'Mo': {
        'K_alpha1': 0.709319,
        'K_alpha2': 0.713609,
        'K_beta': 0.632305,
        'K_alpha_avg': 0.711464,
        'K_alpha2_K_alpha1_ratio': 0.5,
    },
    'Co': {
        'K_alpha1': 1.788996,
        'K_alpha2': 1.792850,
        'K_beta': 1.620830,
        'K_alpha_avg': 1.790260,
        'K_alpha2_K_alpha1_ratio': 0.5,
    },
    'Fe': {
        'K_alpha1': 1.936046,
        'K_alpha2': 1.939980,
        'K_beta': 1.757462,
        'K_alpha_avg': 1.938013,
        'K_alpha2_K_alpha1_ratio': 0.5,
    },
    'Cr': {
        'K_alpha1': 2.289760,
        'K_alpha2': 2.293663,
        'K_beta': 2.084920,
        'K_alpha_avg': 2.291712,
        'K_alpha2_K_alpha1_ratio': 0.5,
    },
    'Ag': {
        'K_alpha1': 0.559420,
        'K_alpha2': 0.563813,
        'K_beta': 0.497082,
        'K_alpha_avg': 0.561617,
        'K_alpha2_K_alpha1_ratio': 0.5,
    },
}


class BrukerRAW4Parser:
    """
    Parser for Siemens/Bruker RAW v4 binary files (single-axis powder diffraction).

    IMPORTANT SCOPE LIMITATION:
    This parser is designed for SINGLE-AXIS theta-2theta powder diffraction scans.
    It extracts one scan axis (typically "Theta") and assumes a standard 1D scan geometry.

    NOT SUPPORTED (without validation):
    - Multi-axis scans (texture, pole figures, reciprocal space maps)
    - Non-standard scan geometries
    - Extraction of additional goniometer angles (omega, chi, phi) if present

    File Structure (RAW v4):
    - Header: 8 bytes ("RAW4.00\x00")
    - Metadata blocks with tags (USER, SITE, SAMPLEID, COMMENT, CREATOR)
    - Measurement parameters (start angle, step, count, scan axis name, etc.)
    - Intensity data (appears as interleaved float32 pairs in tested files)

    Note: Format is proprietary and not fully documented. Offsets and structure
    are reverse-engineered from example files and may not generalize to all RAW v4 variants.
    """

    HEADER_MAGIC = b'RAW4.00\x00'
    HEADER_SIZE = 8

    # Known metadata block tags
    METADATA_TAGS = {
        b'USER\x00\x00\x00\x00': 'user',
        b'SITE\x00\x00\x00\x00': 'site',
        b'SAMPLEID\x00\x00\x00': 'sample_id',
        b'COMMENT\x00': 'comment',
        b'CREATOR\x00': 'creator',
    }

    def __init__(self, filepath: str):
        """
        Initialize parser with RAW file path.

        Args:
            filepath: Path to .raw file
        """
        self.filepath = Path(filepath)
        self.metadata = {}
        self.scan_params = {}
        self.intensities = None
        self.angles = None
        self.logger = logging.getLogger(__name__)

    def parse(self) -> Dict[str, Any]:
        """
        Parse the RAW file and extract all data.

        Returns:
            Dictionary containing metadata, scan parameters, and intensity data

        Raises:
            ValueError: If file format is invalid

        Warnings:
            Logs warning about single-axis limitation - multi-axis scans may not be correctly parsed
        """
        # Log scope limitation warning
        self.logger.warning(
            'RAW v4 parser designed for single-axis theta-2theta scans. '
            'Multi-axis scans (texture, pole figures) may not be correctly parsed. '
            'Additional goniometer angles (omega, chi, phi) are not extracted.'
        )

        with open(self.filepath, 'rb') as f:
            # Read entire file
            full_data = f.read()

        # Verify header
        header = full_data[: self.HEADER_SIZE]
        if header != self.HEADER_MAGIC:
            raise ValueError(
                f'Invalid RAW file header. Expected {self.HEADER_MAGIC!r}, '
                f'got {header!r}'
            )

        # Read file info (next 20 bytes after header)
        # Contains: unknown flag (4), date (8), time (8)
        offset = self.HEADER_SIZE
        _ = struct.unpack('<I', full_data[offset : offset + 4])[
            0
        ]  # unknown_flag (unused)
        offset += 4

        # Date string (null-terminated)
        date_bytes = full_data[offset : offset + 8]
        self.metadata['date'] = date_bytes.rstrip(b'\x00').decode(
            'ascii', errors='ignore'
        )
        offset += 8

        # Time string (null-terminated)
        time_bytes = full_data[offset : offset + 8]
        self.metadata['time'] = time_bytes.rstrip(b'\x00').decode(
            'ascii', errors='ignore'
        )
        offset += 8

        # Parse metadata blocks (using data after header)
        remaining_data = full_data[offset:]
        self._parse_metadata_blocks(remaining_data)

        # Parse measurement parameters and intensity data (needs full file data for absolute offsets)
        self._parse_measurement_data(full_data)

        return {
            'metadata': self.metadata,
            'scan_params': self.scan_params,
            'intensities': self.intensities,
            'angles': self.angles,
        }

    def _parse_metadata_blocks(self, data: bytes):
        """
        Extract metadata blocks (USER, SITE, SAMPLEID, etc.) from binary data.

        Args:
            data: Raw binary data after header
        """
        # Search for known tags
        for tag_bytes, tag_name in self.METADATA_TAGS.items():
            pos = data.find(tag_bytes)
            if pos >= 0:
                # Read until next tag or null terminator
                # Metadata values follow the tag, typically null-terminated
                value_start = pos + len(tag_bytes)
                value_end = data.find(b'\x0a', value_start, value_start + 200)
                if value_end < 0:
                    value_end = data.find(b'\x00\x00', value_start, value_start + 200)

                if value_end > value_start:
                    value = data[value_start:value_end]
                    # Clean up value - strip leading and trailing null bytes, newlines, and spaces
                    value_str = value.strip(b'\x00\x0a\x20').decode(
                        'ascii', errors='ignore'
                    )
                    if value_str:
                        self.metadata[tag_name] = value_str

    def _parse_measurement_data(self, data: bytes):
        """
        Extract scan parameters and intensity array from RAW 4.00 format.

        RAW 4.00 format structure (reverse-engineered):
        - Number of points: uint16 at offset 0x01c3
        - Step size: float64 at offset 0x020f
        - Start angle: float64 at offset 0x04fb
        - Data section: Starts at offset 0x051f (1311 bytes from start)
        - Intensities: Interleaved float32 pairs [intensity, unknown, intensity, unknown, ...]
        - End angle: Calculated from start_angle + (num_points - 1) * step_size

        Args:
            data: Raw binary data (full file contents)
        """
        try:
            # Read number of points from fixed offset (uint16)
            num_points = struct.unpack('<H', data[0x01C3:0x01C5])[0]
            self.scan_params['num_points'] = num_points
            self.logger.info(f'Number of points: {num_points}')
        except Exception as e:
            self.logger.error(f'Could not read num_points at 0x01c3: {e}')
            num_points = None

        try:
            # Read step size from fixed offset (float64)
            step_size = struct.unpack('<d', data[0x020F:0x0217])[0]
            self.scan_params['step_size'] = step_size
            self.logger.info(f'Step size: {step_size:.10f}°')
        except Exception as e:
            self.logger.error(f'Could not read step_size at 0x020f: {e}')
            step_size = None

        try:
            # Read start angle from fixed offset (float64 / double precision)
            start_angle = struct.unpack('<d', data[0x04FB:0x0503])[0]
            self.scan_params['start_angle'] = start_angle
            self.logger.info(f'Start angle: {start_angle}°')
        except Exception as e:
            self.logger.error(f'Could not read start angle at 0x04fb: {e}')

        try:
            # Read scan axis name from fixed offset (null-terminated ASCII string)
            # At offset 0x04D0 we find the axis name (e.g., "Theta")
            axis_bytes = data[0x04D0:0x04E0]
            axis_name = (
                axis_bytes.split(b'\x00')[0].decode('ascii', errors='ignore').strip()
            )
            if axis_name:
                self.scan_params['scan_axis'] = axis_name
                self.logger.info(f'Scan axis: {axis_name}')
        except Exception as e:
            self.logger.error(f'Could not read scan axis at 0x04d0: {e}')
            start_angle = None

        try:
            # Read X-ray tube anode material from fixed offset (null-terminated ASCII string)
            # At offset 0x01A8 we find the anode material (e.g., "Cu", "Mo", "Co")
            # Note: There may be leading null bytes, so we filter them out
            anode_bytes = data[0x01A8:0x01B0]
            # Split on null and find first non-empty string
            anode_material = None
            for segment in anode_bytes.split(b'\x00'):
                if segment:
                    anode_material = segment.decode('ascii', errors='ignore').strip()
                    break

            if anode_material and anode_material in XRAY_WAVELENGTHS:
                self.metadata['anode_material'] = anode_material
                self.logger.info(f'X-ray tube anode: {anode_material}')

                # Add wavelength data from lookup table
                wavelengths = XRAY_WAVELENGTHS[anode_material]
                self.metadata['wavelength_kalpha1'] = wavelengths['K_alpha1']
                self.metadata['wavelength_kalpha2'] = wavelengths['K_alpha2']
                self.metadata['wavelength_kbeta'] = wavelengths['K_beta']
                self.metadata['wavelength_kalpha_avg'] = wavelengths['K_alpha_avg']
                self.metadata['kalpha2_kalpha1_ratio'] = wavelengths[
                    'K_alpha2_K_alpha1_ratio'
                ]
                self.logger.info(f'K-alpha1: {wavelengths["K_alpha1"]} Å')
                self.logger.info(f'K-alpha2: {wavelengths["K_alpha2"]} Å')
            elif anode_material:
                self.logger.warning(
                    f'Unknown anode material "{anode_material}" - wavelengths not available'
                )
                self.metadata['anode_material'] = anode_material
        except Exception as e:
            self.logger.warning(f'Could not read anode material at 0x01a8: {e}')

        # Calculate end angle if we have all required parameters
        if start_angle is not None and step_size is not None and num_points is not None:
            end_angle = start_angle + (num_points - 1) * step_size
            self.scan_params['end_angle'] = end_angle
            self.logger.info(f'Calculated end angle: {end_angle:.6f}°')

            # Generate angle array
            self.angles = [start_angle + i * step_size for i in range(num_points)]

        # Data section starts at fixed offset 0x051f (1311 bytes from start)
        data_start_offset = 0x051F

        # Read intensity values (every other float32, starting from data_offset)
        try:
            data_section = data[data_start_offset:]
            data_section_size = len(data) - data_start_offset
            num_floats = data_section_size // 4

            if num_floats > 0:
                # Unpack all float32 values
                all_floats = struct.unpack(
                    f'<{num_floats}f', data_section[: num_floats * 4]
                )

                # Extract every other float (the actual intensities)
                # Pattern discovered: [intensity1, unknown1, intensity2, unknown2, ...]
                # The unknown values appear to be related but we only need intensities
                self.intensities = [all_floats[i] for i in range(0, len(all_floats), 2)]

                self.logger.info(f'Extracted {len(self.intensities)} intensity values')
                if self.intensities:
                    self.logger.info(
                        f'Intensity range: {min(self.intensities):.2f} to {max(self.intensities):.2f}'
                    )

        except Exception as e:
            self.logger.error(f'Could not parse intensity data: {e}')
            self.intensities = []

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
        lines.append(f'# Start: {self.scan_params.get("start_angle", 0):.3f}°')
        lines.append(f'# Step: {self.scan_params.get("step_size", 0):.4f}°')
        lines.append(f'# Points: {len(self.intensities)}')
        lines.append('#')
        lines.append('# 2Theta(°)  Intensity(counts)')

        for angle, intensity in zip(self.angles, self.intensities):
            lines.append(f'{angle:.4f}  {intensity:.2f}')

        return '\n'.join(lines)


def read_bruker_raw4(filepath: str) -> Dict[str, Any]:
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


def convert_raw_to_xy(raw_filepath: str, xy_filepath: Optional[str] = None) -> str:
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
