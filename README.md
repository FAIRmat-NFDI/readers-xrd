# fairmat-readers-xrd
A package developed by FAIRmat and collaborators that contains file readers for various
X-ray diffraction (XRD) file formats.

## Install and use
To use the readers, install the package using `pip` and import the `read_file`
function. It will check the file extension and, using an appropriate file
reader, return a Python dictionary containing the data.

```sh
pip install fairmat-readers-xrd
```
```py
from fairmat_readers_xrd import read_file

file_path = "<filepath>.xrdml"
data_dict =  read_file(file_path)
```

You can also import individual file readers.
```py
from fairmat_readers_xrd import read_panalytical_xrdml, read_bruker_raw

# Read PANalytical XRDML file
file_path = "<filepath>.xrdml"
data_dict = read_panalytical_xrdml(file_path)

# Read Bruker/Siemens RAW v4 file
raw_file_path = "<filepath>.raw"
raw_data_dict = read_bruker_raw(raw_file_path)
```

Currently, the following file extensions are supported:

| File Extension    | Corresponding Reader Function     | Vendor/Format                            | Notes                                    |
| ----------------- | --------------------------------- | ---------------------------------------- | ---------------------------------------- |
| `.xrdml`          | `read_panalytical_xrdml`          | PANalytical                              |                                          |
| `.rasx`           | `read_rigaku_rasx`                | Rigaku                                   |                                          |
| `.brml`           | `read_bruker_brml`                | Bruker                                   |                                          |
| `.raw`            | `read_bruker_raw`                 | Bruker/Siemens RAW v4                    | One Locked Coupled range; native wavelengths and times |

**RAW4 scope:** one `Locked Coupled` powder range. Other scan modes, multiple
ranges, and older RAW versions are rejected with an explanation. The reader
uses declared block lengths, point counts, and intensity-record strides rather
than absolute positions from one reference file. It reports the detector axis
as `2Theta`, reads native source wavelengths and tube settings, and converts
range step time from milliseconds to seconds. Zero wavelengths remain unreported;
a zero Kα2/Kα1 intensity ratio remains valid. No timezone is invented from the
RAW date/time strings. The block layout follows the primary
[GSAS-II RAW reader](https://gsas-ii.readthedocs.io/en/latest/_modules/GSASII/imports/G2pwd_BrukerRAW.html).

Kα1/Kα2 lookup values are shared with NOMAD Measurements through
`fairmat_readers_xrd.wavelengths.K_ALPHA_WAVELENGTHS`, in angstroms. They retain
the published precision of [International Tables for Crystallography (2006),
Vol. C, Table 4.2.2.1](https://onlinelibrary.wiley.com/iucr/itc/Cb/ch4o2v0001/sec4o2o2/).
These references supply fallback estimates; valid wavelengths and ratios read
from RAW4 or other formats take priority. A reported monochromated spectrum is
not replaced by a tabulated Kα doublet.

## Measurement arrays

Intensity and original counts retain every point and scan, including constant
spectra, identical rows, and repeated outer scan positions. Only shared axes and
optional metadata are compacted. Equivalent axis units are converted before
comparison and stacking; a single detector channel retains the RSM scan dimension.

XRDML may mix counts and already corrected intensities across scans. Corrected
intensity remains available in scan order; incomplete counts or attenuation-factor
arrays are omitted rather than filled with invented values. The original file
retains those partial metadata values. Partially missing geometric coordinates
are rejected before classification; they cannot become shared angles. Completely
absent optional axes remain optional.

RASX profile coordinates follow `ScanInformation.AxisName`. A pure Omega profile
produces `Omega`, while recorded fixed 2θ remains separate. `TwoTheta`,
`TwoThetaOmega`, and `TwoThetaTheta` use a 2θ abscissa. Coordinates and units are
read per profile before common values are compacted. No axis offsets are applied.
The Omega mapping is covered by a synthetic file; an independent real Omega
reference and confirmation of the Rigaku correction convention remain missing.

For BRML point scans, recorded `ScanAxes` coordinates take priority over a
reconstructed start/stop raster. Scanned Theta/Omega, Chi and Phi take priority
over fixed motor positions. Detector-channel coordinates without recorded point
positions still use the existing reconstruction; their channel-spacing convention
has not been independently established.

## Nominal counting times

`countTime` contains nominal integration/counting times in seconds per measurement
step or detector acquisition. Common values remain one-element vectors; varying
scan times retain scan order. Do not multiply them by the number of detector
channels or interpret them as total duration or measured detector live time.

| Format | Time source |
| --- | --- |
| XRDML | `commonCountingTime` |
| RASX STEP | `Speed` when `SpeedUnit` is a time unit |
| RASX CONTINUOUS | Absolute `Step / Speed`, with compatible units |
| BRML | `TimePerStepEffective`, interpreted in seconds |
| RAW4 | Range-header step time, converted from milliseconds |

Missing, nonpositive, nonfinite, or incompatible optional times are omitted.
If timing is incomplete across scans, omit the time array rather than assigning
another scan's time. This metadata does not change intensity normalization.
The Rigaku mode distinction follows the
[General Measurement manual, page 3](https://alexzevalkink.wordpress.com/wp-content/uploads/2017/04/meas_khp2_00068_eng.pdf#page=7).
The BRML choice follows the use of effective time in
[GSAS-II's BRML reader](https://gsas-ii.readthedocs.io/en/latest/_modules/GSASII/imports/G2pwd_BrukerBRML.html).

For XRDML, RASX, and BRML, metadata includes the first scan's start timestamp and
the last scan's end timestamp, when present. Their supplied timezone information
is retained. No duration or total counting time is calculated.

## Development
The package is still under development. To contribute, start with simply raising an
[issue](https://github.com/FAIRmat-NFDI/readers-xrd/issues) on the repo.

### Contributing code
In case you want to contribute a piece of code, or a whole new reader, clone the repo in
your local and start a virtual Python environment inside the directory.
```sh
git clone git@github.com:FAIRmat-NFDI/readers-xrd.git
cd readers-xrd
python3 -m venv .pyenv
source .pyenv/bin/activate
```

Then install the package in editable mode (-e flag), with `dev` dependencies.
You can also run the `pytest` to confirm all went well.
As a recommended practice, you can also install the pre-commit hook for linting (more on
this [here](#ruff-pre-commit-hook)).
```sh
python -m pip install --upgrade pip
pip install -e .[dev]
pytest
pre-commit install
```

Now you can start tinkering around on a new branch (or a [fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo)).
If you want to contribute the code back to the upstream, simply create a
[Pull Request (PR)](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request),
and we will get back to you. As a good practice, add tests for new readers and
make sure existing tests (`pytest`) are passing before creating PR.

### Ruff
We have placed GitHub actions for checking the linting. This can lead to your PR failing
the CI tests in case your code is not correctly linted or formatted.
Best way to avoid this is to use [Ruff](https://docs.astral.sh/ruff/tutorial/) for
automatic linting and formatting.
This is already part of the `dev` environment and you can run Ruff before committing:
```sh
ruff check --fix
ruff format
```
### Ruff pre-commit hook
Additionally, we also provide [pre-commit](https://pre-commit.com) hook for Ruff. This
will run Ruff every time you try to commit and raise errors (and provide fixes)
in case your code needs linting or formatting. You simply have to add these fixes and
commit again.

To use this pre-commit hook, make sure to run the following in the terminal after you
clone the repo:
```sh
pre-commit install
```
