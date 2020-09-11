#!/usr/bin/python
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

import math
from typing import List

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


class DifferentialExpression(object):

    @staticmethod
    def _fold_change(expression_tumor: float, expression_reference: float) -> float:
        """
        This function determines the classical fold change between tumour and reference transcript expression.
        Log2(expr in tumor / expr in reference)
        """
        fold_change = None
        try:
            fold_change = math.log(expression_tumor / expression_reference, 2)
        except (TypeError, ValueError, ZeroDivisionError) as e:
            pass
        return fold_change

    @staticmethod
    def _percentile_calc(expression_tumor: float, expression_reference_sum: float) -> float:
        """
        This function calculates the expression difference between tumour and reference data in form of a percentile
        value.
        expr in tumor * 100 / (sum of expr in ref tissue + 1)
        """
        percentile = None
        try:
            percentile = (expression_tumor * 100) / (expression_reference_sum + 1)
        except (TypeError, ValueError, ZeroDivisionError) as e:
            pass
        return percentile

    @staticmethod
    def _pepper_calc(expression_tumor: float, expression_reference: float, expression_reference_sd: float):
        """
        This function calculates the expression difference between tumour and reference data based on Pepper
        publication, in a z-score similar manner.
        expr in tumour - mean epxr in reference tissue / standard deviation of expression in refernce
        """
        pepper = None
        try:
            pepper = (expression_tumor - expression_reference) / expression_reference_sd
        except (TypeError, ValueError, ZeroDivisionError) as e:
            pass
        return pepper

    def get_annotations(self, expression_tumor, expression_reference, expression_reference_sum,
                        expression_reference_sd) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(name="Differential_expression_log2_foldchange", value=self._fold_change(
                expression_tumor=expression_tumor, expression_reference=expression_reference)),
            AnnotationFactory.build_annotation(name="Differential_expression_percentile",
                                               value=self._percentile_calc(
                                                   expression_tumor=expression_tumor,
                                                   expression_reference_sum=expression_reference_sum)),
            AnnotationFactory.build_annotation(name="Differential_expression_pepper", value=self._pepper_calc(
                expression_tumor=expression_tumor, expression_reference=expression_reference,
                expression_reference_sd=expression_reference_sd))
        ]
