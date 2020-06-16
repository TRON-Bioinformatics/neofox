#!/usr/bin/env python

from logzero import logger

import input.Tcell_predictor.prediction as prediction
from input.helpers import intermediate_files


class TcellPrediction:

    def __init__(self, references):
        self.tcell_prediction_score = "NA"
        self.tcell_prediction_score_9mer = "NA"
        self.references = references

    def _triple_gen_seq_subst_for_prediction(self, gene, substitution, epitope, score, threshold):
        """
        extracts gene id, epitope sequence and substitution from epitope dictionary
        Tcell predictor works with 9mers only! --> extract for 9mers only
        """
        logger.debug("{} {} {} {} {}".format(gene, epitope, substitution, score, str(len(epitope))))
        result = (["NA", "NA", "NA"])
        if str(len(epitope)) == str(9):
            z = [gene.replace(" ", ""), epitope, substitution]
            if threshold is None:
                z = [gene.replace(" ", ""), epitope, substitution]
                result = (z)
            else:
                if float(score) < threshold:
                    z = [gene.replace(" ", ""), epitope, substitution]
                    result = (z)
        return result

    def _write_triple_to_file(self, triple, tmpfile_in):
        """
        writes triple (gene id, epitope sequence, substitution) to temporary file
        """
        with open(tmpfile_in, "w") as f:
            tripleString = " ".join(triple)
            f.write(tripleString + "\n")

    def _prediction_single_mps(self, tmpfile_in, tmpfile_out):
        """
        calls T cell predictor tool to perform predictions; returns
        """
        # TODO: the interface for prediction needs to be simplified, it does not make sense to write a file to
        # TODO: return a single score
        prediction.main(tmpfile_in, tmpfile_out, self.references)
        logger.debug(tmpfile_out)
        with open(tmpfile_out, "r") as f:
            l_prediction = f.readlines(0)
            logger.debug(l_prediction)
        try:
            score = l_prediction[-1].split(",")[-1].rstrip("\n")
        except IndexError:
            score = "indefinable_by_TcellPredictor"

        return score

    def _wrapper_tcellpredictor(self, gene, substitution, epitope, score, threshold, tmpfile_in, tmpfile_out):
        """
        wrapper function to determine
        """
        trp = self._triple_gen_seq_subst_for_prediction(
            gene=gene, substitution=substitution, epitope=epitope, score=score, threshold=threshold)
        logger.debug(trp)
        pred_out = "NA"
        if "NA" not in trp:
            self._write_triple_to_file(trp, tmpfile_in)
            pred_out = self._prediction_single_mps(tmpfile_in, tmpfile_out)
        return pred_out

    def calculate_tcell_predictor_score(self, gene, substitution, epitope, score, threshold=None):
        ''' returns Tcell_predictor score given mps in dictionary format
                '''
        tmp_tcellPredIN = intermediate_files.create_temp_file(prefix="tmp_TcellPredicIN_", suffix=".txt")
        tmp_tcellPredOUT = intermediate_files.create_temp_file(prefix="tmp_TcellPredicOUT_", suffix=".txt")
        return self._wrapper_tcellpredictor(
            gene=gene, substitution=substitution, epitope=epitope, score=score, threshold=threshold,
            tmpfile_in=tmp_tcellPredIN, tmpfile_out=tmp_tcellPredOUT)
