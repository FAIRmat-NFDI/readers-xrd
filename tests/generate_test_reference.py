#!/usr/bin/env python3
"""
Generate reference JSON file for RAW parser testing.

Usage:
    python generate_test_reference.py input.raw [output.json]

This script reads a RAW file and generates a reference JSON that can be used
for regression testing. The JSON contains the expected output from read_rigaku_raw().
"""

import sys
import json
import argparse
from pathlib import Path

# Import parser (adjust path if needed)
try:
    from fairmat_readers_xrd import read_rigaku_raw
except ImportError:
    print('Error: Cannot import fairmat_readers_xrd')
    print("Make sure you're in the correct environment and the package is installed")
    sys.exit(1)

import pint


def convert_quantity_to_string(obj):
    """
    Convert pint Quantity objects to strings for JSON serialization.

    Args:
        obj: Object to convert (Quantity, dict, list, or other)

    Returns:
        JSON-serializable version of the object
    """
    if isinstance(obj, pint.Quantity):
        return str(obj)
    elif isinstance(obj, dict):
        return {k: convert_quantity_to_string(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_quantity_to_string(item) for item in obj]
    return obj


def generate_reference(raw_file: Path, output_file: Path = None):
    """
    Generate reference JSON from a RAW file.

    Args:
        raw_file: Path to input RAW file
        output_file: Path to output JSON file (default: raw_file.json)
    """
    if not raw_file.exists():
        print(f'Error: File not found: {raw_file}')
        sys.exit(1)

    if output_file is None:
        output_file = raw_file.with_suffix(raw_file.suffix + '.json')

    print(f'Reading RAW file: {raw_file}')
    try:
        data = read_rigaku_raw(str(raw_file))
    except Exception as e:
        print(f'Error reading RAW file: {e}')
        sys.exit(1)

    print(f'Successfully read RAW file')
    print(f"  - Data points: {len(data.get('2Theta', [[]])[0])}")
    print(f"  - Metadata keys: {list(data.get('metadata', {}).keys())}")

    # Convert to JSON-serializable format
    print('Converting to JSON format...')
    json_data = convert_quantity_to_string(data)

    # Save to file
    print(f'Writing reference to: {output_file}')
    with open(output_file, 'w') as f:
        json.dump(json_data, f, indent=2)

    print('Done!')
    print(f'\nTo use this reference in tests:')
    print(f'  1. Copy {output_file.name} to tests/data/')
    print(f'  2. Update test to compare against this reference')


def main():
    parser = argparse.ArgumentParser(
        description='Generate reference JSON for RAW parser testing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate reference from RAW file (output: input.raw.json)
  python generate_test_reference.py input.raw

  # Specify output file
  python generate_test_reference.py input.raw reference.json

  # Process file in tests directory
  python generate_test_reference.py tests/data/test_sample.raw
        """,
    )
    parser.add_argument('raw_file', type=Path, help='Input RAW file')
    parser.add_argument(
        'output_file',
        type=Path,
        nargs='?',
        help='Output JSON file (default: input.raw.json)',
    )

    args = parser.parse_args()

    generate_reference(args.raw_file, args.output_file)


if __name__ == '__main__':
    main()
