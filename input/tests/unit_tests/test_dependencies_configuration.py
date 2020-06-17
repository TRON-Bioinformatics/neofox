import os
import unittest
from unittest import TestCase

import input
import input.tests.unit_tests.tools as test_tools
from input.exceptions import INPuTConfigurationException
from input.references import DependenciesConfiguration


class TestDependenciesConfiguration(TestCase):

    def setUp(self):
        self.variables = {
            input.INPUT_BLASTP_ENV: '/path/to/blastp',
            input.INPUT_NETMHC2PAN_ENV: '/path/to/netmhc2pan',
            input.INPUT_NETMHCPAN_ENV: '/path/to/netmhcpan',
            input.INPUT_RSCRIPT_ENV: '/path/to/rscript',
            input.INPUT_MIXMHCPRED_ENV: '/path/to/mixmhcpred',
            input.INPUT_MIXMHC2PRED_ENV: '/path/to/mixmhc2pred'
        }
        self.non_existing = '/path/to/nothing'
        test_tools._mock_file_existence(
            existing_files=self.variables.values(),
            unexisting_files=[self.non_existing]
        )

    def _load_env_variables(self):
        for k, v in self.variables.items():
            os.environ[k] = v

    def test_not_provided_variable(self):
        self._load_env_variables()
        for v in self.variables.keys():
            del os.environ[v]
            with self.assertRaises(INPuTConfigurationException):
                DependenciesConfiguration()

    def test_empty_string_variable(self):
        self._load_env_variables()
        for v in self.variables.keys():
            os.environ[v] = ""
            with self.assertRaises(INPuTConfigurationException):
                DependenciesConfiguration()

    def test_non_existing_variable(self):
        self._load_env_variables()
        for v in self.variables.keys():
            os.environ[v] = self.non_existing
            with self.assertRaises(INPuTConfigurationException):
                DependenciesConfiguration()

    def test_all_resources_exist(self):
        self._load_env_variables()
        config = DependenciesConfiguration()
        self.assertTrue(config.blastp == self.variables[input.INPUT_BLASTP_ENV])
        self.assertTrue(config.mix_mhc2_pred == self.variables[input.INPUT_MIXMHC2PRED_ENV])
        self.assertTrue(config.mix_mhc_pred == self.variables[input.INPUT_MIXMHCPRED_ENV])
        self.assertTrue(config.rscript == self.variables[input.INPUT_RSCRIPT_ENV])
        self.assertTrue(config.net_mhc_pan == self.variables[input.INPUT_NETMHCPAN_ENV])
        self.assertTrue(config.net_mhc2_pan == self.variables[input.INPUT_NETMHC2PAN_ENV])


if __name__ == "__main__":
    unittest.main()
