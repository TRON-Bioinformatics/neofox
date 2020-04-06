import unittest
from unittest import TestCase
import os
import input
from input.predict_all_epitopes import Bunchepitopes
from input.exceptions import INPuTConfigurationException


class TestPredictAllEpitopes(TestCase):

    def test_reference_environment_variable_is_required(self):
        # del os.environ[input.REFERENCE_FOLDER_ENV]
        with self.assertRaises(INPuTConfigurationException):
            Bunchepitopes()

    def test_empty_reference_folder_fails(self):
        os.environ[input.REFERENCE_FOLDER_ENV] = 'dummy'
        with self.assertRaises(INPuTConfigurationException):
            Bunchepitopes()


if __name__ == "__main__":
    unittest.main()
