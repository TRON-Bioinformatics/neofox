#!/usr/bin/env python

import sys
import tempfile
import input.Tcell_predictor.prediction as prediction
from logzero import logger


class Tcellprediction:

    def __init__(self, references):
        self.TcellPrdictionScore = "NA"
        self.TcellPrdictionScore_9merPred = "NA"
        self.references = references

    def _triple_gen_seq_subst_for_prediction(self, props, all = True, affinity = False):
        """
        extracts gene id, epitope sequence and substitution from epitope dictionary
        Tcell predictor works with 9mers only! --> extract for 9mers only
        """
        if "gene.x" in props:
            gene = props["gene.x"]
        else:
            gene = props["gene"]
        subst = props["substitution"]
        if affinity:
            epi = props["best_affinity_epitope_netmhcpan4_9mer"]
            score = props["best_affinity_netmhcpan4_9mer"]
        else:
            epi = props["MHC_I_epitope_.best_prediction."]
            score = props["MHC_I_score_.best_prediction."]
        logger.debug("{} {} {} {} {}".format(gene, epi, subst, score, str(len(epi))))
        if str(len(epi)) == str(9):
            z = [gene.replace(" ", ""), epi, subst]
            if all:
                z = [gene.replace(" ", ""), epi, subst]
                return(z)
            else:
                if(affinity):
                    if float(score) < 500:
                        z = [gene.replace(" ", ""), epi, subst]
                        return(z)
                    else:
                        return(["NA", "NA", "NA"])
                else:
                    if float(score) < 2:
                        z = [gene.replace(" ", ""), epi, subst]
                        return(z)
                    else:
                        return(["NA", "NA", "NA"])
        else:
            return(["NA", "NA", "NA"])

    def _write_triple_to_file(self, triple, tmpfile_in):
        """
        writes triple (gene id, epitope sequence, substitution) to temporary file
        """
        with open(tmpfile_in,"w") as f:
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

    def _wrapper_tcellpredictor(self, props, tmpfile_in, tmpfile_out, all = True, affinity = False):
        """
        wrapper function to determine
        """
        trp = self._triple_gen_seq_subst_for_prediction(props, all, affinity)
        logger.debug(trp)
        pred_out = "NA"
        if "NA" not in trp:
            self._write_triple_to_file(trp, tmpfile_in)
            pred_out = self._prediction_single_mps(tmpfile_in, tmpfile_out)
        return pred_out

    def main(self, props):
        ''' returns Tcell_predictor score given mps in dictionary format
        '''
        tmp_tcellPredIN_file = tempfile.NamedTemporaryFile(prefix ="tmp_TcellPredicIN_", suffix = ".txt", delete = False)
        tmp_tcellPredIN = tmp_tcellPredIN_file.name
        tmp_tcellPredOUT_file = tempfile.NamedTemporaryFile(prefix ="tmp_TcellPredicOUT_", suffix = ".txt", delete = False)
        tmp_tcellPredOUT = tmp_tcellPredOUT_file.name
        # returns score for all epitopes --> no filtering based on mhc affinity here!
        self.TcellPrdictionScore = self._wrapper_tcellpredictor(props, tmp_tcellPredIN, tmp_tcellPredOUT)
        tmp_tcellPredIN_file = tempfile.NamedTemporaryFile(prefix ="tmp_TcellPredicIN_", suffix = ".txt", delete = False)
        tmp_tcellPredIN = tmp_tcellPredIN_file.name
        tmp_tcellPredOUT_file = tempfile.NamedTemporaryFile(prefix ="tmp_TcellPredicOUT_", suffix = ".txt", delete = False)
        tmp_tcellPredOUT = tmp_tcellPredOUT_file.name
        # returns score for all epitopes --> do filtering based on mhc affinity here (threshold 500 nM)!
        self.TcellPrdictionScore_9merPred = self._wrapper_tcellpredictor(
            props, tmp_tcellPredIN, tmp_tcellPredOUT, all=False, affinity=True)
