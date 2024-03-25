# Developer guide

> Information: The build system changed in version 1.1.1 to poetry (https://python-poetry.org/).

## Installation of neofox via poetry

First setup poetry (see [Setup poetry](#setup-poetry)).

To install neofox into the virtual environment managed by poetry (if a `poetry.lock` file is given, the given versions are used. Otherwise the dependecies listed in pyproject.toml are used to generate the environment):

```
poetry install
```

If the **wheel** archive is needed, they can be generated using (generates them in the `dist` directory):

```
poetry build
```

To use the installed package use `poetry run` e.g.:

```
poetry run neofox --help
```

## Installation via pip

The wheel file has to be generated as described in [Installation of neofox via poetry](#installation-of-neofox-via-poetry) (`poetry build`).

```
pip install dist/neofox-x.y.z.whl
```

## Run integration tests

To run the integration tests make sure you have a file `.env` that contains the environment variables described in the configuration section.

Run the integration tests as follows:
```
python -m unittest discover neofox.tests.integration_tests
```

The integration tests run over some real datasets and they take some time to run.

The integration test that runs the whole program over a relevant dataset can be run as follows:
```
python -m unittest neofox.tests.integration_tests.test_neofox
```

### Regression tests

This last test (ie: `test_neofox`) writes its output to a file named `neofox/tests/resources/output_yyyymmddHHMMSS.txt`. If there is an existing file named `neofox/tests/resources/output_previous.txt` then it loads both files in memory and compares them. It outputs whether there are some lost or gained columns and for the common columns it evaluates if the values are the same. If they are the same the file `output_previous.txt` is overwritten by the new file, otherwise it outputs the details of the differing columns.

## Run unit tests

The unit tests do not have any dependency and they finish in seconds.

Run the unit tests as follows:
```
python -m unittest discover neofox.tests.unit_tests
```

## Logging

Logs are written to the standard error and to the output folder by default. Optionally they can be written to a file by setting the environment variable `NEOFOX_LOGFILE` pointing to the desired file.


## Build the protocol buffers models

The protocol buffers model rely on the betterproto library. Install it as follows:
```
pip install "betterproto[compiler]"
```

The models and the required scripts are in the folder `neofox/models`

Build the models into Python code with `make models`.

Build the HTML documentation with `make html` (this requires docker).

## Setup poetry

Use mamba to install poetry

```
mamba create -n poetry conda-forge::poetry
```

Required python version for the neofox package also installed via mamba

```
# installs python v3.7.12
mamba create -n python3.7 conda-forge::python=3.7
```

Set the python executable for poetry. This generates the python venv `~/.cache/pypoetry/virtualenvs/neofox...`.

```
poetry env use ~/.conda/envs/mamba/envs/python3.7/bin/python3.7
```

To get information on the currently used python environment variable run. This will also show the location of the virtual environment in which the package is or will be installed.

```
poetry env info
```
