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
import unittest
from unittest import TestCase

from neofox.helpers.runner import Runner


class TestRunner(TestCase):
    def setUp(self):
        self.runner = Runner()

    def test_runner(self):
        output, errors = self.runner.run_command(cmd=["python", "-V"])
        self.assertTrue("Python 3.7" in output or "Python 3.6" in output or "Python 3.8" in output)
        self.assertTrue(len(errors) == 0)

    def test_runner_failure(self):
        with self.assertRaises(Exception):
            self.runner.run_command(cmd=["nocommandwiththisname"])


if __name__ == "__main__":
    unittest.main()
