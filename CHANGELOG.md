# Changelog

## Unreleased

### Added

- Return nominal acquisition times, scan-boundary timestamps and supported
  user/software/optics metadata with source paths. Unreported values remain absent.

### Changed

- Preserve declared signal units: count-scale input stays dimensionless and
  explicit rates use `count/s`. Reject unknown units and mixed count/rate XRDML scans.
- **Breaking:** RAW4 reports `2Theta` coordinates, native wavelengths and nominal
  step times in seconds. Consumers using the former `Theta` label must adapt.
  Shared Kα references provide fallback values when needed.

### Fixed

- Respect RAW4 block lengths, record sizes and point counts; reject unsupported
  modes, multiple ranges and truncated input.
- Preserve repeated spectra and single-point shapes. Follow per-profile RASX axes
  and recorded BRML coordinates; omit incomplete optional arrays without inventing data.
