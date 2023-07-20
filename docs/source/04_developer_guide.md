# Developer guide

## Build the package

To build the package just run:
```
python setup.py bdist_wheel
```

This will create an installable wheel file under `dist/neofox-x.y.z.whl`.

## Install the package

Install the wheel file as follows:
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

Build the HTML documentation with `make html` (this requires membership of the docker group).
