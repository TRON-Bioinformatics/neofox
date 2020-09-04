import os
import unittest
from unittest import TestCase

import neofox
from neofox.exceptions import NeofoxConfigurationException
from neofox.neofox import NeoFox


class TestPredictAllEpitopes(TestCase):

    def test_reference_environment_variable_is_required(self):
        # del os.environ[neofox.REFERENCE_FOLDER_ENV]
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(neoantigens=None, patient_id=None, patients=None, num_cpus=1)

    def test_empty_reference_folder_fails(self):
        os.environ[neofox.REFERENCE_FOLDER_ENV] = 'dummy'
        with self.assertRaises(NeofoxConfigurationException):
            NeoFox(neoantigens=None, patient_id=None, patients=None, num_cpus=1)


if __name__ == "__main__":
    unittest.main()
