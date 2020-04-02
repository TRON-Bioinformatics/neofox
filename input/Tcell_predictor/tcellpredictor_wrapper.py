#!/usr/bin/env python

import os
import sys
import tempfile
import subprocess

# my_path = os.path.abspath(os.path.dirname(__file__))
# my_path2 = "/".join(my_path.split("/")[0:-1])
# sys.path.insert(0, my_path2)
# sys.path.insert(0, my_path)


class Tcellprediction:

    def __init__(self):
        self.TcellPrdictionScore = "NA"
        self.TcellPrdictionScore_9merPred = "NA"

    def _triple_gen_seq_subst_for_prediction(self, props, all = True, affinity = False):
        """ extracts gene id, epitope sequence and substitution from epitope dictionary
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
        print >> sys.stderr, gene, epi, subst, score, str(len(epi))
        if str(len(epi)) == str(9):
            z = [gene.replace(" ", ""), epi, subst]
            #return(z)
            #print("hal")
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
        '''writes triple (gene id, epitope sequence, substitution) to temporary file
        '''
        with open(tmpfile_in,"w") as f:
            tripleString = " ".join(triple)
            f.write(tripleString + "\n")

    def _prediction_single_mps(self, tmpfile_in, tmpfile_out):
        '''calls T cell predictor tool to perform predictions; returns
        '''
        # TODO: refactor this so we call directly once Python 3 migration is done
        pred_tool = "/".join([".", "prediction.py" ])
        cmd = " ".join(["/code/Anaconda/3/2018/bin/python", pred_tool, tmpfile_in, tmpfile_out])
        p = subprocess.Popen(cmd.split(" "),stderr=subprocess.PIPE,stdout=subprocess.PIPE)

        lines = p.stdout
        stdoutdata, stderrdata = p.communicate()
        print >> sys.stderr, stderrdata
        print >> sys.stderr, tmpfile_out
        with open(tmpfile_out, "r") as f:
            l_prediction =  f.readlines(0)
            print >> sys.stderr, l_prediction
        try:
            score = l_prediction[-1].split(",")[-1].rstrip("\n")
        except IndexError:
            score = "indefinable_by_TcellPredictor"

        return(score)

    def _wrapper_tcellpredictor(self, props, tmpfile_in, tmpfile_out, all = True, affinity = False):
        '''wrapper function to determine
        '''
        trp = self._triple_gen_seq_subst_for_prediction(props, all, affinity)
        print >> sys.stderr, trp
        #print trp
        if "NA" not in trp:
            self._write_triple_to_file(trp, tmpfile_in)
            pred_out = self._prediction_single_mps(tmpfile_in, tmpfile_out)
            return(pred_out)
        else:
            return("NA")

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
