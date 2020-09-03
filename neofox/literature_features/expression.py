#!/usr/bin/env python
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
            AnnotationFactory.build_annotation(name="Expression_Mutated_Transcript", value=self.expression),
            AnnotationFactory.build_annotation(name="Expression_Mutated_Transcript_tumor_content",
                                               value=self.expression_tumor_content_corrected)
        ]

