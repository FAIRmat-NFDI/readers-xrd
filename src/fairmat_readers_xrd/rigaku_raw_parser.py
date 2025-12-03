"""
Pure Python parser for Rigaku RAW 4.00 X-ray diffraction files.

This module provides native Python parsing of Rigaku's proprietary binary .raw
format without requiring external tools, Wine, or .NET libraries.

Based on reverse engineering of RAW 4.00 file structure.

Author: Generated for NOMAD
License: MIT
"""

import struct
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
import logging


class RigakuRAW4Parser:
    """
    Parser for Rigaku RAW version 4.00 binary files.
    
    File Structure (RAW 4.00):
    - Header: 8 bytes ("RAW4.00\x00")
    - Metadata blocks with tags (USER, SITE, SAMPLEID, COMMENT, CREATOR)
    - Measurement parameters (start angle, step, count, etc.)
    - Intensity data (float array)
    """
    
    HEADER_MAGIC = b'RAW4.00\x00'
    HEADER_SIZE = 8
    
    # Known metadata block tags
    METADATA_TAGS = {
        b'USER\x00\x00\x00\x00': 'user',
        b'SITE\x00\x00\x00\x00': 'site',
        b'SAMPLEID': 'sample_id',
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
        """
        with open(self.filepath, 'rb') as f:
            # Read entire file
            full_data = f.read()
        
        # Verify header
        header = full_data[:self.HEADER_SIZE]
        if header != self.HEADER_MAGIC:
            raise ValueError(
                f'Invalid RAW file header. Expected {self.HEADER_MAGIC!r}, '
                f'got {header!r}'
            )
        
        # Read file info (next 20 bytes after header)
        # Contains: unknown flag (4), date (8), time (8)
        offset = self.HEADER_SIZE
        _ = struct.unpack('<I', full_data[offset:offset+4])[0]  # unknown_flag (unused)
        offset += 4
        
        # Date string (null-terminated)
        date_bytes = full_data[offset:offset+8]
        self.metadata['date'] = date_bytes.rstrip(b'\x00').decode('ascii', errors='ignore')
        offset += 8
        
        # Time string (null-terminated)
        time_bytes = full_data[offset:offset+8]
        self.metadata['time'] = time_bytes.rstrip(b'\x00').decode('ascii', errors='ignore')
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
                    # Clean up value
                    value_str = value.rstrip(b'\x00\x0a\x20').decode('ascii', errors='ignore')
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
            num_points = struct.unpack('<H', data[0x01c3:0x01c5])[0]
            self.scan_params['num_points'] = num_points
            self.logger.info(f"Number of points: {num_points}")
        except Exception as e:
            self.logger.error(f"Could not read num_points at 0x01c3: {e}")
            num_points = None
        
        try:
            # Read step size from fixed offset (float64)
            step_size = struct.unpack('<d', data[0x020f:0x0217])[0]
            self.scan_params['step_size'] = step_size
            self.logger.info(f"Step size: {step_size:.10f}°")
        except Exception as e:
            self.logger.error(f"Could not read step_size at 0x020f: {e}")
            step_size = None
        
        try:
            # Read start angle from fixed offset (float64 / double precision)
            start_angle = struct.unpack('<d', data[0x04fb:0x0503])[0]
            self.scan_params['start_angle'] = start_angle
            self.logger.info(f"Start angle: {start_angle}°")
        except Exception as e:
            self.logger.error(f"Could not read start angle at 0x04fb: {e}")
            start_angle = None
        
        # Calculate end angle if we have all required parameters
        if start_angle is not None and step_size is not None and num_points is not None:
            end_angle = start_angle + (num_points - 1) * step_size
            self.scan_params['end_angle'] = end_angle
            self.logger.info(f"Calculated end angle: {end_angle:.6f}°")
            
            # Generate angle array
            self.angles = [start_angle + i * step_size for i in range(num_points)]
        
        # Data section starts at fixed offset 0x051f (1311 bytes from start)
        data_start_offset = 0x051f
        
        # Read intensity values (every other float32, starting from data_offset)
        try:
            data_section = data[data_start_offset:]
            data_section_size = len(data) - data_start_offset
            num_floats = data_section_size // 4
            
            if num_floats > 0:
                # Unpack all float32 values
                all_floats = struct.unpack(f'<{num_floats}f', data_section[:num_floats*4])
                
                # Extract every other float (the actual intensities)
                # Pattern discovered: [intensity1, unknown1, intensity2, unknown2, ...]
                # The unknown values appear to be related but we only need intensities
                self.intensities = [all_floats[i] for i in range(0, len(all_floats), 2)]
                
                self.logger.info(f"Extracted {len(self.intensities)} intensity values")
                if self.intensities:
                    self.logger.info(
                        f"Intensity range: {min(self.intensities):.2f} to {max(self.intensities):.2f}"
                    )
                
        except Exception as e:
            self.logger.error(f"Could not parse intensity data: {e}")
            self.intensities = []

    def to_xy_format(self) -> str:
        """
        Convert parsed data to simple XY format (angle, intensity).
        
        Returns:
            String in XY format
        """
        if not self.intensities or not self.angles:
            raise ValueError('No data parsed. Call parse() first.')
        
        lines = ['# Rigaku RAW 4.00 converted to XY format']
        lines.append(f'# Sample: {self.metadata.get("sample_id", "Unknown")}')
        lines.append(f'# Date: {self.metadata.get("date", "")} {self.metadata.get("time", "")}')
        lines.append(f'# Start: {self.scan_params.get("start_angle", 0):.3f}°')
        lines.append(f'# Step: {self.scan_params.get("step_size", 0):.4f}°')
        lines.append(f'# Points: {len(self.intensities)}')
        lines.append('#')
        lines.append('# 2Theta(°)  Intensity(counts)')
        
        for angle, intensity in zip(self.angles, self.intensities):
            lines.append(f'{angle:.4f}  {intensity:.2f}')
        
        return '\n'.join(lines)


def read_rigaku_raw4(filepath: str) -> Dict[str, Any]:
    """
    Read a Rigaku RAW 4.00 file.
    
    Args:
        filepath: Path to .raw file
    
    Returns:
        Dictionary with 'metadata', 'scan_params', 'intensities', 'angles'
    
    Example:
        >>> data = read_rigaku_raw4('sample.raw')
        >>> print(data['metadata']['sample_id'])
        >>> print(len(data['intensities']))
    """
    parser = RigakuRAW4Parser(filepath)
    return parser.parse()


def convert_raw_to_xy(raw_filepath: str, xy_filepath: Optional[str] = None) -> str:
    """
    Convert Rigaku .raw file to simple XY format.
    
    Args:
        raw_filepath: Path to input .raw file
        xy_filepath: Optional output path. If None, replaces .raw with .xy
    
    Returns:
        Path to created XY file
    """
    parser = RigakuRAW4Parser(raw_filepath)
    parser.parse()
    
    xy_content = parser.to_xy_format()
    
    if xy_filepath is None:
        xy_filepath = str(Path(raw_filepath).with_suffix('.xy'))
    
    with open(xy_filepath, 'w') as f:
        f.write(xy_content)
    
    return xy_filepath
