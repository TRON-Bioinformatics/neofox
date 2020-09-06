#!/usr/bin/env python
import os
import pickle
from typing import List

from neofox.helpers import intermediate_files
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.predictors.Tcell_predictor.preprocess import Preprocessor
from neofox.predictors.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder

CLASSIFIER_PICKLE = 'Classifier.pickle'


class TcellPrediction:

    def __init__(self):
        self.tcell_prediction_score = None
        self.tcell_prediction_score_9mer = None
        with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), CLASSIFIER_PICKLE), 'rb') as f:
            self.classifier = pickle.load(f)

    def _triple_gen_seq_subst(self, gene, substitution, epitope, score, threshold):
        """
        extracts gene id, epitope sequence and substitution from epitope dictionary
        Tcell predictor works with 9mers only! --> extract for 9mers only
        """
        result = None
        if str(len(epitope)) == str(9):
            if threshold is None or float(score) < threshold:
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
        if scores is not None and len(scores) > 0 and scores[-1] is not None and len(scores[-1]) > 0:
            # it returns the last number from the latest entry in the list
            result = str(scores[-1][-1])
        return result

    def _wrapper_tcellpredictor(self, gene, substitution, epitope, score, threshold, tmpfile_in):
        """
        wrapper function to determine
        """
        trp = self._triple_gen_seq_subst(
            gene=gene, substitution=substitution, epitope=epitope, score=score, threshold=threshold)
        pred_out = None
        if trp is not None:
            self._write_triple_to_file(trp, tmpfile_in)
            pred_out = self._run_prediction(tmpfile_in)
        return pred_out

    def _calculate_tcell_predictor_score(self, gene, substitution, epitope, score, threshold=None):
        ''' returns Tcell_predictor score given mps in dictionary format
                '''
        tmp_tcellPredIN = intermediate_files.create_temp_file(prefix="tmp_TcellPredicIN_", suffix=".txt")
        return self._wrapper_tcellpredictor(
            gene=gene, substitution=substitution, epitope=epitope, score=score, threshold=threshold,
            tmpfile_in=tmp_tcellPredIN)

    def get_annotations(self, gene, substitution, netmhcpan: BestAndMultipleBinder) -> List[Annotation]:

        return [
            AnnotationFactory.build_annotation(value=self._calculate_tcell_predictor_score(
                gene=gene, substitution=substitution, epitope=netmhcpan.mhcI_affinity_epitope_9mer,
                score=netmhcpan.mhcI_affinity_9mer),
                name="Tcell_predictor_score"),
            AnnotationFactory.build_annotation(value=self._calculate_tcell_predictor_score(
                gene=gene, substitution=substitution, epitope=netmhcpan.mhcI_affinity_epitope_9mer,
                score=netmhcpan.mhcI_affinity_9mer, threshold=500),
                name="Tcell_predictor_score_cutoff500nM")
        ]
