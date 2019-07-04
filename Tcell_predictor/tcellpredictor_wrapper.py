#!/usr/bin/env python

#import Prediction2
import os
import sys
from Bio import SeqIO
from argparse import ArgumentParser
import tempfile
import subprocess

my_path = os.path.abspath(os.path.dirname(__file__))
my_path2 = "/".join(my_path.split("/")[0:-1])
sys.path.insert(0, my_path2)
sys.path.insert(0, my_path)


from helpers import data_import



class Tcellprediction:
    def __init__(self):
        self.TcellPrdictionScore = "NA"

    def triple_gen_seq_subst_for_prediction(self, props, all = False):
        """ extracts gene id, epitope sequence and substitution from epitope dictionary
        Tcell predictor works with 9mers only! --> extract for 9mers only
        """
        if "gene.x" in props:
            gene = props["gene.x"]
        else:
            gene = props["gene"]
        epi = props["MHC_I_epitope_.best_prediction."]
        subst = props["substitution"]
        length = props["MHC_I_peptide_length_.best_prediction."]
        score = props["MHC_I_score_.best_prediction."]
        #print gene, epi, subst, length, score
        if length == str(9):
            if all:
                z = [gene.replace(" ", ""), epi, subst]
                return(z)
            else:
                if float(score < 2):
                    z = [gene.replace(" ", ""), epi, subst]
                    return(z)
                else:
                    return(["NA", "NA", "NA"])
        else:
            return(["NA", "NA", "NA"])


    def write_triple_to_file(self, triple, tmpfile_in):
        '''writes triple (gene id, epitope sequence, substitution) to temporary file
        '''
        with open(tmpfile_in,"w") as f:
            tripleString = " ".join(triple)
            f.write(tripleString + "\n")


    def prediction_single_mps(self, tmpfile_in, tmpfile_out, path_to_Tcell_predictor):
        '''calls T cell predictor tool to perform predictions; returns
        '''
        pred_tool = "/".join([path_to_Tcell_predictor, "prediction.py" ])
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

    def wrapper_tcellpredictor(self, props, tmpfile_in, tmpfile_out, path_to_Tcell_predictor):
        '''wrapper function to determine
        '''
        trp = self.triple_gen_seq_subst_for_prediction(props, all = True)
        #print trp
        if "NA" not in trp:
            self.write_triple_to_file(trp, tmpfile_in)
            pred_out = self.prediction_single_mps(tmpfile_in, tmpfile_out, path_to_Tcell_predictor)
            return(pred_out)
        else:
            return("NA")


    def full_dataset(self, dat_tup, all = False):
        '''returns from icam output gene id, mhc I epitope and substitution for each mps
        pre-filtering for mps of length 9 and mhc I binding score < 2
        '''
        head = dat_tup[0]
        dat = dat_tup[1]
        dat_triple = []
        if "gene.x" in head:
            gene_col = head.index("gene.x")
        else:
            gene_col = head.index("gene")
        epi_col = head.index("MHC_I_epitope_.best_prediction.")
        subst_col = head.index("substitution")
        length_col = head.index("MHC_I_peptide_length_.best_prediction.")
        score_col = head.index("MHC_I_score_.best_prediction.")
        for ii,i in enumerate(dat):
            if i[length_col] == str(9):
                if all:
                    z = [i[gene_col].replace(" ", ""), i[epi_col], i[subst_col]]
                    dat_triple.append(z)
                else:
                    if float(i[score_col]) < 2:
                        z = [i[gene_col], i[epi_col], i[subst_col]]
                        dat_triple.append(z)
        return dat_triple

    def write_ouptut_to_file(self, epitope_data):
        '''
        This function prints output, semilicon separated --> txt file!!!!
        '''
        for ii,i in enumerate(epitope_data):
              print " ".join(i).lstrip(" ")

    def main(self, props):
        ''' returns Tcell_predictor score given mps in dictionary format
        '''
        path_to_Tcell_predictor = my_path
        tmp_tcellPredIN_file = tempfile.NamedTemporaryFile(prefix ="tmp_TcellPredicIN_", suffix = ".txt", delete = False)
        tmp_tcellPredIN = tmp_tcellPredIN_file.name
        tmp_tcellPredOUT_file = tempfile.NamedTemporaryFile(prefix ="tmp_TcellPredicOUT_", suffix = ".txt", delete = False)
        tmp_tcellPredOUT = tmp_tcellPredOUT_file.name
        self.TcellPrdictionScore = self.wrapper_tcellpredictor(props, tmp_tcellPredIN, tmp_tcellPredOUT, path_to_Tcell_predictor)





if __name__ == '__main__':
    # if full icam output table is passed to script
    '''
    f = sys.argv[1]
    dat = data_import.import_dat_icam(f, indel = False)
    #print dat
    #print full_dataset(dat)
    l = full_dataset(dat, all = True)
    write_ouptut_to_file(l)
    '''

    # test for input implementation
    import epitope
    import predict_all_epitopes

    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
    dat = data_import.import_dat_icam(file, False)
    if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
        dat = data_import.change_col_names(dat)

    path_to_Tcell_predictor = my_path

    for ii,i in enumerate(dat[1]):
        if ii < 10:
            print ii
            dict_epi = epitope.Epitope()
            dict_epi.init_properties(dat[0], dat[1][ii])
            #print dict_epi.properties
            tcellpred = Tcellprediction()

            tcellpred.main(dict_epi.properties)
            print tcellpred.TcellPrdictionScore
