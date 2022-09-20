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
import os
import pickle
import warnings
from typing import List
from neofox.model.validation import ModelValidator
from neofox.model.neoantigen import Annotation, Neoantigen, PredictedEpitope
from neofox.model.factories import AnnotationFactory
from neofox.published_features.Tcell_predictor.preprocess import Preprocessor
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)

CLASSIFIER_PICKLE = "Classifier.pickle"


class TcellPrediction:

    def __init__(self):
        # UserWarning: Trying to unpickle estimator DecisionTreeClassifier from version 0.19.0 when using version
        # 0.20.3. This might lead to breaking code or invalid results. Use at your own risk.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with open(
                os.path.join(
                    os.path.abspath(os.path.dirname(__file__)), CLASSIFIER_PICKLE
                ),
                'rb'
            ) as f:
                # this sets the n_jobs parameter otherwise inherited from the pickle file
                self.classifier = pickle.load(f).set_params(n_jobs=1)

        self.preprocessor = Preprocessor()

    def _wrapper_tcellpredictor(self, gene, epitope: PredictedEpitope):
        """
        wrapper function to determine
        """
        result = None
        has_gene = gene is not None and gene.strip() != ""
        if has_gene and len(epitope.mutated_peptide) == 9:
            mat = self.preprocessor.main(gene, epitope=epitope)
            scores = self.classifier.predict_proba(mat)
            result = "indefinable_by_TcellPredictor"
            if (
                    scores is not None
                    and len(scores) > 0
                    and scores[-1] is not None
                    and len(scores[-1]) > 0
            ):
                # it returns the last number from the latest entry in the list
                result = str(scores[-1][-1])
        return result

    def calculate_tcell_predictor_score(
        self, gene: str, epitope: PredictedEpitope
    ):
        """returns Tcell_predictor score given mps in dictionary format"""
        tcell_predictor_score = None
        if not ModelValidator.has_peptide_rare_amino_acids(epitope.mutated_peptide):
            tcell_predictor_score = self._wrapper_tcellpredictor(gene=gene, epitope=epitope)
        return tcell_predictor_score

    def get_annotations(
        self, neoantigen: Neoantigen, netmhcpan: BestAndMultipleBinder
    ) -> List[Annotation]:
        # TODO: this is difficult to extend to more complex mutations (eg: MNVs, indels) as only considers first mutated
        #  position
        tcell_predictor_score = None
        if neoantigen.wild_type_xmer and netmhcpan.best_ninemer_epitope_by_affinity.mutated_peptide:
            tcell_predictor_score = self.calculate_tcell_predictor_score(
                gene=neoantigen.gene,
                epitope=netmhcpan.best_ninemer_epitope_by_affinity)
        annotations = [
            AnnotationFactory.build_annotation(
                value=tcell_predictor_score,
                name="Tcell_predictor",
            )
        ]
        return annotations

    def get_annotations_epitope_mhci(self, epitope: PredictedEpitope, gene: str) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(
                value=self.calculate_tcell_predictor_score(epitope=epitope, gene=gene),
                name='Tcell_predictor')
            ]
