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
from fairmat_readers_xrd import read_panalytical_xrdml

file_path = "<filepath>.xrdml"
data_dict =  read_panalytical_xrdml(file_path)
```

Currently, the following file extensions are supported:

| File Extension    | Corresponding Reader Function     |
| ----------------- | --------------------------------- |
| `.xrdml`          | `read_panalytical_xrdml`          |
| `.rasx`           | `read_rigaku_rasx`                |
| `.brml`           | `read_bruker_brml`                |

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
