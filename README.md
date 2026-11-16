# fairmat-readers-xrd

Readers for X-ray diffraction files from PANalytical (`.xrdml`),
Rigaku (`.rasx`) and Bruker (`.brml`, RAW v4 `.raw`).

## Install and use

```bash
pip install fairmat-readers-xrd
```

```python
from fairmat_readers_xrd import read_file

data = read_file("measurement.xrdml")
intensity = data["intensity"]
```

Readers return a dictionary containing measurement arrays, physical quantities
with units, and available metadata. Original scan and point order is preserved.
Count-scale signals are dimensionless; explicitly declared rates use `count/s`.
Optional `countTime` values are nominal acquisition times in seconds, not total
measurement duration or detector live time. Missing metadata remains absent.

RAW support is limited to one RAW4 `Locked Coupled` powder range.
Other RAW versions, multiple ranges and unsupported modes are rejected.

## Development

Install an editable checkout with development dependencies:

```bash
pip install -e '.[dev]'
pytest
ruff check .
ruff format . --check
```

Include relevant tests when contributing. Open issues and pull requests on
[GitHub](https://github.com/FAIRmat-NFDI/readers-xrd).
