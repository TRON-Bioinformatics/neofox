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

from neofox.expression_imputation.expression_imputation import ExpressionAnnotator


class TestExpressionAnnotator(TestCase):

    def setUp(self) -> None:
        self.expression_annotator = ExpressionAnnotator()

    def test_expression_annotation(self):

        # test retrieval of cohort index
        self.assertEqual(22, self.expression_annotator.cohort_indices["SARC"])
        # test gene expression
        result = self.expression_annotator.get_gene_expression_annotation(
            gene_name="ATF2", tcga_cohort="SARC"
        )
        self.assertAlmostEqual(5.034754, result, places=5)
        result = self.expression_annotator.get_gene_expression_annotation(
            gene_name="NBPF24", tcga_cohort="KIPAN"
        )
        self.assertAlmostEqual(6.839773310, result, places=5)
        result = self.expression_annotator.get_gene_expression_annotation(
            gene_name="blabluplupp", tcga_cohort="KIPAN"
        )
        self.assertIsNone(result)
        result = self.expression_annotator.get_gene_expression_annotation(
            gene_name="NBPF24", tcga_cohort="blabluplupp"
        )
        self.assertIsNone(result)
        result = self.expression_annotator.get_gene_expression_annotation(
            gene_name="blabluplupp", tcga_cohort="blabluplupp"
        )
        self.assertIsNone(result)

    def test_empty_gene(self):
        result = self.expression_annotator.get_gene_expression_annotation(
            gene_name="", tcga_cohort="SARC"
        )
        self.assertIsNone(result)
        result = self.expression_annotator.get_gene_expression_annotation(
            gene_name=None, tcga_cohort="SARC"
        )
        self.assertIsNone(result)

    def test_empty_cohort(self):
        result = self.expression_annotator.get_gene_expression_annotation(
            gene_name="ATF2", tcga_cohort=None
        )
        self.assertIsNone(result)
        result = self.expression_annotator.get_gene_expression_annotation(
            gene_name="ATF2", tcga_cohort=""
        )
        self.assertIsNone(result)

