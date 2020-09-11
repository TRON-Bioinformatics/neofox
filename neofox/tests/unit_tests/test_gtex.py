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

from neofox.annotation_resources.gtex.gtex import GTEx


class TestGtex(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.gtex = GTEx()

    def test_get_metrics(self):
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("BRCA2", "skin")
        self.assertEqual(0.146363636363636, mean_expression)
        self.assertEqual(27.37, sum_expression)
        self.assertEqual(0.0950531054943768, sd_expression)
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("CFTR", "bladder")
        self.assertEqual(0.35454545454545505, mean_expression)
        self.assertEqual(3.9, sum_expression)
        self.assertEqual(1.11291835851839, sd_expression)

    def test_get_metrics_from_non_existing_tissue(self):
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("BRCA2", "chinichinchin")
        self.assertEqual(None, mean_expression)
        self.assertEqual(None, sum_expression)
        self.assertEqual(None, sd_expression)
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("BRCA2", None)
        self.assertEqual(None, mean_expression)
        self.assertEqual(None, sum_expression)
        self.assertEqual(None, sd_expression)

    def test_get_metrics_from_non_existing_gene(self):
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics("NOOOOOOOPE", "skin")
        self.assertEqual(None, mean_expression)
        self.assertEqual(None, sum_expression)
        self.assertEqual(None, sd_expression)
        mean_expression, sum_expression, sd_expression = self.gtex.get_metrics(None, "skin")
        self.assertEqual(None, mean_expression)
        self.assertEqual(None, sum_expression)
        self.assertEqual(None, sd_expression)


