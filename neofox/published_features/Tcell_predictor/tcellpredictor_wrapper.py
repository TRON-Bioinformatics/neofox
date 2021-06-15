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
from neofox.helpers import intermediate_files
from neofox.model.conversion import ModelValidator
from neofox.model.neoantigen import Annotation, Neoantigen
from neofox.model.wrappers import AnnotationFactory
from neofox import AFFINITY_THRESHOLD_DEFAULT
from neofox.published_features.Tcell_predictor.preprocess import Preprocessor
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)

CLASSIFIER_PICKLE = "Classifier.pickle"


class TcellPrediction:

    def __init__(self, affinity_threshold=AFFINITY_THRESHOLD_DEFAULT):
        self.affinity_threshold = affinity_threshold
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
                self.classifier = pickle.load(f)

    def _triple_gen_seq_subst(self, gene, substitution, epitope, score):
        """
        extracts gene id, epitope sequence and substitution from epitope dictionary
        Tcell predictor works with 9mers only! --> extract for 9mers only
        """
        result = None
        has_gene = gene is not None and gene.strip() != ""
        if has_gene and len(epitope) == 9:
            if self.affinity_threshold is None or float(score) < self.affinity_threshold:
                result = [gene.replace(" ", ""), epitope, substitution]
        return result

    def _write_triple_to_file(self, triple, tmpfile_in):
        """
        writes triple (gene id, epitope sequence, substitution) to temporary file
        """
        with open(tmpfile_in, "w") as f:
            tripleString = " ".join(triple)
            f.write(tripleString + "\n")

    def _run_prediction(self, f_name):
        input_file = f_name
        mat = Preprocessor().main(input_file)
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

    def _wrapper_tcellpredictor(
        self, gene, substitution, epitope, score, tmpfile_in
    ):
        """
        wrapper function to determine
        """
        trp = self._triple_gen_seq_subst(
            gene=gene,
            substitution=substitution,
            epitope=epitope,
            score=score
        )
        pred_out = None
        if trp is not None:
            self._write_triple_to_file(trp, tmpfile_in)
            pred_out = self._run_prediction(tmpfile_in)
        return pred_out

    def _calculate_tcell_predictor_score(
        self, gene, substitution, epitope, score
    ):
        """returns Tcell_predictor score given mps in dictionary format"""
        tmp_tcellPredIN = intermediate_files.create_temp_file(
            prefix="tmp_TcellPredicIN_", suffix=".txt"
        )
        tcell_predictor_score = None
        if not ModelValidator.has_peptide_rare_amino_acids(epitope):
            tcell_predictor_score = self._wrapper_tcellpredictor(
                gene=gene, substitution=substitution, epitope=epitope, score=score, tmpfile_in=tmp_tcellPredIN, )
        return tcell_predictor_score

    def get_annotations(
        self, neoantigen: Neoantigen, netmhcpan: BestAndMultipleBinder
    ) -> List[Annotation]:
        # TODO: this is difficult to extend to more complex mutations (eg: MNVs, indels) as only considers first mutated
        #  position
        tcell_predictor_score = None
        if neoantigen.mutation.wild_type_xmer and netmhcpan.best_ninemer_epitope_by_affinity.peptide:
            mutation_position = neoantigen.mutation.position[0]
            wild_type_aminoacid = neoantigen.mutation.wild_type_xmer[
                mutation_position - 1
            ]  # it is 1-based
            mutated_aminoacid = neoantigen.mutation.mutated_xmer[mutation_position - 1]
            tcell_predictor_score = self._calculate_tcell_predictor_score(gene=neoantigen.gene,
                                                          substitution=wild_type_aminoacid + mutated_aminoacid,
                                                          epitope=netmhcpan.best_ninemer_epitope_by_affinity.peptide,
                                                          score=netmhcpan.best_ninemer_epitope_by_affinity.affinity_score)
        annotations = [
            AnnotationFactory.build_annotation(
                value=tcell_predictor_score,
                name="Tcell_predictor_score",
            )
        ]
        return annotations
