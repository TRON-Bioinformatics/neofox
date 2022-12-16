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

from neofox.model.factories import NOT_AVAILABLE_VALUE, NeoantigenFactory
from neofox.model.neoantigen import Neoantigen
from neofox.published_features.expression import Expression


class TestExpression(TestCase):

    def setUp(self) -> None:
        self.expression = Expression()

    def test_calculate_expression_mutation(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            rna_expression=12.0, dna_variant_allele_frequency=0.2, patient_identifier="patient1",
            mutated_xmer="DDDDD")
        result = self.expression.get_annotations(neoantigen=neoantigen)[0]
        self.assertGreater(float(result.value), 0.0)

        # no reads for mut
        neoantigen = NeoantigenFactory.build_neoantigen(
            rna_expression=12.0, dna_variant_allele_frequency=0.0, patient_identifier="patient1",
            mutated_xmer="DDDDD")
        result = self.expression.get_annotations(neoantigen=neoantigen)[0]
        self.assertEqual(result.value, "0")

        # no reads for mut/wt
        neoantigen = NeoantigenFactory.build_neoantigen(
            rna_expression=12.0, dna_variant_allele_frequency=-1, rna_variant_allele_frequency=-1, patient_identifier="patient1",
            mutated_xmer="DDDDD")
        result = self.expression.get_annotations(neoantigen=neoantigen)[0]
        self.assertEqual(result.value, NOT_AVAILABLE_VALUE)

        neoantigen = NeoantigenFactory.build_neoantigen(
            rna_expression=None, dna_variant_allele_frequency=-1, rna_variant_allele_frequency=-1,
            patient_identifier="patient1", mutated_xmer="DDDDD")
        result = self.expression.get_annotations(neoantigen=neoantigen)[0]
        self.assertEqual(result.value, NOT_AVAILABLE_VALUE)
