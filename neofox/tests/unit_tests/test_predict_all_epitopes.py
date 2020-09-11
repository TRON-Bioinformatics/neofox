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
