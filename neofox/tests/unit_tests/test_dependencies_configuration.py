#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
import os
import unittest
from unittest import TestCase

import neofox
import neofox.tests.unit_tests.tools as test_tools
from neofox.exceptions import NeofoxConfigurationException
from neofox.references.references import DependenciesConfiguration


class TestDependenciesConfiguration(TestCase):

    def setUp(self):
        self.variables = {
            neofox.NEOFOX_BLASTP_ENV: '/path/to/blastp',
            neofox.NEOFOX_NETMHC2PAN_ENV: '/path/to/netmhc2pan',
            neofox.NEOFOX_NETMHCPAN_ENV: '/path/to/netmhcpan',
            neofox.NEOFOX_RSCRIPT_ENV: '/path/to/rscript',
            neofox.NEOFOX_MIXMHCPRED_ENV: '/path/to/mixmhcpred',
            neofox.NEOFOX_MIXMHC2PRED_ENV: '/path/to/mixmhc2pred'
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
            with self.assertRaises(NeofoxConfigurationException):
                DependenciesConfiguration()

    def test_empty_string_variable(self):
        self._load_env_variables()
        for v in self.variables.keys():
            os.environ[v] = ""
            with self.assertRaises(NeofoxConfigurationException):
                DependenciesConfiguration()

    def test_non_existing_variable(self):
        self._load_env_variables()
        for v in self.variables.keys():
            os.environ[v] = self.non_existing
            with self.assertRaises(NeofoxConfigurationException):
                DependenciesConfiguration()

    def test_all_resources_exist(self):
        self._load_env_variables()
        config = DependenciesConfiguration()
        self.assertTrue(config.blastp == self.variables[neofox.NEOFOX_BLASTP_ENV])
        self.assertTrue(config.mix_mhc2_pred == self.variables[neofox.NEOFOX_MIXMHC2PRED_ENV])
        self.assertTrue(config.mix_mhc_pred == self.variables[neofox.NEOFOX_MIXMHCPRED_ENV])
        self.assertTrue(config.rscript == self.variables[neofox.NEOFOX_RSCRIPT_ENV])
        self.assertTrue(config.net_mhc_pan == self.variables[neofox.NEOFOX_NETMHCPAN_ENV])
        self.assertTrue(config.net_mhc2_pan == self.variables[neofox.NEOFOX_NETMHC2PAN_ENV])


if __name__ == "__main__":
    unittest.main()
