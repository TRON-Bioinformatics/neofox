#!/usr/bin/env python
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
from typing import List

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory


class Expression:

    def __init__(self, transcript_expression, vaf_rna, tumor_content):
        self.expression = self._get_expression_annotation(transcript_expression, vaf_rna)
        self.expression_tumor_content_corrected = self._get_expression_tumor_content_corrected_annotation(
            self.expression, tumor_content)

    @staticmethod
    def _get_expression_annotation(transcript_expression: float, vaf_rna: float) -> float:
        """
        This function calculates the product of VAF in RNA and transcript expression
        to reflect the expression of the mutated transcript
        """
        expression_mut = None
        try:
            expression_mut = transcript_expression * vaf_rna if vaf_rna is not None and vaf_rna >= 0.0 else None
        except (TypeError, ValueError):
            pass
        return expression_mut

    @staticmethod
    def _get_expression_tumor_content_corrected_annotation(expression_mutation: float, tumor_content: float) -> float:
        """calculated expression of mutation corrected by tumour content"""
        expression_mut_tc = None
        try:
            expression_mut_tc = expression_mutation / tumor_content
        except (TypeError, ZeroDivisionError) as e:
            pass
        return expression_mut_tc

    def get_annotations(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(name="Expression_mutated_transcript", value=self.expression),
            AnnotationFactory.build_annotation(name="Expression_mutated_transcript_tumor_content",
                                               value=self.expression_tumor_content_corrected)
        ]

