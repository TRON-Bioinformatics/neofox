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

from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import PredictedEpitope, Annotation


class EpitopeHelperTest(TestCase):
    def test_number_of_mismatches(self):
        # one mismatch
        result = EpitopeHelper.number_of_mismatches("FIAGLIAIV", "FIAGDIAIV")
        self.assertEqual(result, 1)
        # three mismatch
        result = EpitopeHelper.number_of_mismatches("FIAGLIAIV", "FLAGDIAIN")
        self.assertEqual(result, 3)
        # no mismatch
        result = EpitopeHelper.number_of_mismatches("FIAGLIAIV", "FIAGLIAIV")
        self.assertEqual(result, 0)
        # wt sequence shorter than mut sequence
        result = EpitopeHelper.number_of_mismatches("FIAGI", "FIAGLIAIV")
        self.assertEqual(result, 1)
        # mut sequence shorter than wt sequence
        result = EpitopeHelper.number_of_mismatches("FIAGLIAIV", "FIAGI")
        self.assertEqual(result, 1)

    def test_position_mutation(self):
        position = EpitopeHelper().position_of_mutation_epitope(
            PredictedEpitope(wild_type_peptide="AAAAAA", mutated_peptide="AAANAA"))
        self.assertEqual(position, 4)
        position = EpitopeHelper().position_of_mutation_epitope(
            PredictedEpitope(wild_type_peptide="AAAAAA", mutated_peptide="AAAAAA"))
        self.assertEqual(position, -1)
        position = EpitopeHelper().position_of_mutation_epitope(
            PredictedEpitope(wild_type_peptide="AAAAAA", mutated_peptide="AANNNN"))
        self.assertEqual(position, 6)

    def test_get_annotation_by_name(self):
        annotations = [Annotation(name='this', value='5'), Annotation(name='that', value='0')]
        self.assertEqual(EpitopeHelper.get_annotation_by_name(annotations=annotations, name='this'), '5')
        self.assertEqual(EpitopeHelper.get_annotation_by_name(annotations=annotations, name='that'), '0')
        try:
            EpitopeHelper.get_annotation_by_name(annotations=annotations, name='nothing')
            self.assertTrue(False)
        except ValueError:
            self.assertTrue(True)

    # TODO: test ther methods in the EpitopeHelper
