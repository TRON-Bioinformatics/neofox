import os

# TODO: change this import when we move to python3
# from unittest.mock import Mock
from mock import Mock


def _mock_file_existence(existing_files=[], unexisting_files=[]):
    original_os_path_exists = os.path.exists

    def side_effect(filename):
        if filename in existing_files:
            return True
        elif filename in unexisting_files:
            return False
        else:
            return original_os_path_exists(filename)

    os.path.exists = Mock(side_effect=side_effect)


def head(file_name, n=10):
    with open(file_name) as myfile:
        try:
            for x in range(n):
                print((next(myfile)))
        except StopIteration:
            pass


def print_and_delete(filename):
    head(filename)
    os.remove(filename)
