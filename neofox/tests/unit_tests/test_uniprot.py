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
from unittest import TestCase

from neofox.annotation_resources.uniprot.uniprot import Uniprot
import pkg_resources
import neofox.tests


class TestUniprot(TestCase):
    @classmethod
    def setUpClass(cls):
        cls.uniprot = Uniprot(
            pkg_resources.resource_filename(
                neofox.tests.__name__,
                "resources/uniprot_human_with_isoforms.first200linesfortesting.pickle",
            )
        )

    def test_sequence_not_in_uniprot(self):
        result = self.uniprot.is_sequence_not_in_uniprot("NOT_IN_UNIPROT")
        self.assertEqual(True, result)

    def test_sequence_in_uniprot(self):
        result = self.uniprot.is_sequence_not_in_uniprot("LLEKVKAHEIAWLHGTI")
        self.assertEqual(False, result)
