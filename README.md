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
[issue](https://github.com/FAIRmat-NFDI/readers-xrd/issues) on the repo. In case you
want to contribute a piece of code, or a whole new reader, clone the repo in your local
and start a virtual Python environment inside the directory.
```sh
git clone git@github.com:FAIRmat-NFDI/readers-xrd.git
cd readers-xrd
python3.9 -m venv .pyenv
source .pyenv/bin/activate
```

Then install the package in editable mode (-e flag), with `dev` dependencies. 
You can also run the `pytest` to confirm all went well.
```sh
python -m pip install --upgrade pip
pip install -e .[dev]
pytest
```

Now you can start tinkering around. We recommend to create a new branch for this. 
If you want to contribute the code back to the upstream, simply create a 
[Pull Request (PR)](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request),
and we will get back to you. As a good practice, add tests for new readers and
make sure existing tests are passing before creating PR.

You can also use `Ruff` for automatic linting by running the following:
```sh
ruff format .
```
