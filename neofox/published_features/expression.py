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
from neofox.model.neoantigen import Annotation, Neoantigen
from neofox.model.factories import AnnotationFactory


class Expression:

    @staticmethod
    def _get_expression_annotation(
        transcript_gene_expression: float, vaf: float
    ) -> float:
        """
        This function calculates the product of VAF and transcript expression
        to reflect the expression of the mutated transcript
        """
        expression_mut = None
        try:
            expression_mut = (
                transcript_gene_expression * vaf
                if vaf is not None and vaf >= 0.0
                else None
            )
        except (TypeError, ValueError):
            pass
        return expression_mut

    def get_annotations(self, neoantigen: Neoantigen) -> List[Annotation]:

        return [
            AnnotationFactory.build_annotation(
                name="Mutated_rnaExpression_fromRNA", value=self._get_expression_annotation(
                    transcript_gene_expression=neoantigen.rna_expression, vaf=neoantigen.rna_variant_allele_frequency)),
            AnnotationFactory.build_annotation(
                name="Mutated_rnaExpression_fromDNA", value=self._get_expression_annotation(
                    transcript_gene_expression=neoantigen.rna_expression, vaf=neoantigen.dna_variant_allele_frequency)),
            AnnotationFactory.build_annotation(
                name="Mutated_imputedGeneExpression_fromRNA", value=self._get_expression_annotation(
                    transcript_gene_expression=neoantigen.imputed_gene_expression, vaf=neoantigen.rna_variant_allele_frequency)),
            AnnotationFactory.build_annotation(
                name="Mutated_imputedGeneExpression_fromDNA", value=self._get_expression_annotation(
                    transcript_gene_expression=neoantigen.imputed_gene_expression, vaf=neoantigen.dna_variant_allele_frequency))
        ]
