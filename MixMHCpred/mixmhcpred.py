#!/usr/bin/env python

import subprocess
import os
import sys
import tempfile

my_path = os.path.abspath(os.path.dirname(__file__))
my_path2 = "/".join(my_path.split("/")[0:-1])
sys.path.insert(0, my_path2)
sys.path.insert(0, my_path)

from helpers import data_import
import timeit


class MixMHCpred:
    def __init__(self):
        self.all_peptides = "NA"
        self.all_scores = "NA"
        self.all_ranks = "NA"
        self.all_alleles = "NA"
        self.best_peptide = "NA"
        self.best_score = "NA"
        self.best_rank = "NA"
        self.best_allele = "NA"
        self.best_peptide_wt = "NA"
        self.best_score_wt = "NA"
        self.best_rank_wt = "NA"
        self.difference_score_mut_wt = "NA"

    def mut_position_xmer_seq(self, props):
        '''returns position of mutation in xmer sequence
        '''
        xmer_wt = props["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        xmer_mut = props["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        p1 = -1
        for i,aa in enumerate(xmer_mut):
            if aa != xmer_wt[i]:
                p1 = i + 1
        return str(p1)

    def generate_nmers(self, props, list_lengths, mut = True ):
        ''' generates peptides covering mutation of all lengts that are provided. Returns peptides as list
        '''
        xmer_wt = props["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        xmer_mut = props["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        list_peptides = []
        pos_mut = int(self.mut_position_xmer_seq(props))
        long_seq = xmer_mut if mut else xmer_wt
        for l in list_lengths:
            l = int(l)
            start_first = pos_mut - (l)
            starts = []
            for s in range(l):
                starts.append(int(start_first + s))
            ends = []
            [ends.append(int(s + (l))) for s in starts]
            for s,e in zip(starts, ends):
                list_peptides.append(long_seq[s:e])
        return list_peptides


    def generate_fasta(self, seqs, tmpfile):
        ''' Writes seqs given in seqs list into fasta file
        '''
        #fastafile = my_path + "/tmp.fasta"
        counter = 0
        with open(tmpfile,"w") as f:
            for seq in seqs:
                id = "".join([">seq", str(counter)])
                f.write(id + "\n")
                f.write(seq + "\n")
                counter += 1

    def get_hla_allels(self, props, hla_patient_dict):
        ''' returns hla allele of patients given in hla_file
        '''
        if "patient.id" in props:
            patientid = props["patient.id"]
        else:
            try:
                patientid = props["patient"]
            except KeyError:
                patientid = props["patient.x"]
        return hla_patient_dict[patientid]

    def mixmhcprediction(self, hla_alleles, tmpfasta, outtmp):
        ''' Performs MixMHCpred prediction for desired hla allele and writes result to temporary file.
        '''
        allels_for_prediction = []
        for allele in hla_alleles:
            allele = allele.replace("*", "")
            allele = allele.replace("HLA-", "")
            allels_for_prediction.append(allele)
        hla_allele = ",".join(allels_for_prediction)
        #
        cmd = "MixMHCpred -a " + hla_allele + " -i " + tmpfasta + " -o " + outtmp
        p = subprocess.Popen(cmd.split(" "),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        p_return = p.communicate()

    def read_mixmhcpred(self, outtmp):
        '''imports output of MixMHCpred prediction
        '''
        counter = 0
        header = []
        dat = []
        with open(outtmp) as f:
            for line in f:
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith("#"):
                        continue
                    if line.startswith("Peptide"):
                        counter += 1
                        header = line.split()
                        continue
                    line = line.split()
                    dat.append(line)
        return header, dat


    def extract_best_per_pep(self, pred_dat):
        '''extract info of best allele prediction for all potential ligands per muatation
        '''
        head = pred_dat[0]
        dat = pred_dat[1]
        peps = []
        scores = []
        alleles = []
        ranks = []
        pepcol = head.index("Peptide")
        scorecol = head.index("Score_bestAllele")
        allelecol = head.index("BestAllele")
        rankcol = head.index("%Rank_bestAllele")
        min_value = -1000000000000000000
        for ii,i in enumerate(dat):
            col_of_interest = [i[pepcol],i[scorecol], i[rankcol], i[allelecol]]
            # all potential peptides per mutation --> return ditionary
            peps.append(i[pepcol])
            scores.append(i[scorecol])
            ranks.append(i[rankcol])
            alleles.append(i[allelecol])
        return {"Peptide": peps, "Score_bestAllele": scores,"BestAllele": alleles, "%Rank_bestAllele": ranks }

    def extract_best_peptide_per_mutation(self, pred_dat):
        '''extract best predicted ligand per mutation
        '''
        head = pred_dat[0]
        dat = pred_dat[1]
        peps = []
        scores = []
        alleles = []
        ranks = []
        pepcol = head.index("Peptide")
        scorecol = head.index("Score_bestAllele")
        allelecol = head.index("BestAllele")
        rankcol = head.index("%Rank_bestAllele")
        min_value = -1000000000000000000
        for ii,i in enumerate(dat):
            col_of_interest = [str(i[pepcol]),str(i[scorecol]), str(i[rankcol]), str(i[allelecol])]
            # best ligand per mutation
            if i[scorecol] > min_value:
                min_value = i[scorecol]
                min_pep = col_of_interest
        head_new = ["Peptide", "Score_bestAllele", "%Rank_bestAllele", "BestAllele" ]
        return head_new, min_pep

    def add_best_epitope_info(self, epitope_tuple, column_name):
        '''returns desired information of prediction of best epitope from netmhcpan output;
        e.g. "%Rank": MHC I score, "HLA": HLA allele, "Icore": best epitope
        '''
        dat_head = epitope_tuple[0]
        dat = epitope_tuple[1]
        val = dat_head.index(column_name)
        try:
            return dat[val]
        except IndexError:
            return "NA"


    def extract_WT_for_best(self, props, best_mut_seq):
        '''extracts the corresponding WT epitope for best predicted mutated epitope
        '''
        xmer_wt = props["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        xmer_mut = props["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        start = xmer_mut.find(best_mut_seq)
        l = len(best_mut_seq)
        wt_epi = xmer_wt[start:(start+l)]
        return(wt_epi)

    def extract_WT_info(self, epitope_tuple, column_name):
        '''
        '''
        dat_head = epitope_tuple[0]
        dat = epitope_tuple[1][0]
        val = dat_head.index(column_name)
        try:
            return dat[val]
        except IndexError:
            return "NA"

    def difference_score(self, mut_score, wt_score):
        ''' calcualated difference in MixMHCpred scores between mutated and wt
        '''
        try:
            return str(float(mut_score) - float(wt_score))
        except ValueError:
            return "NA"



    def main(self, props_dict, dict_patient_hla):
        '''Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR
        '''
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix ="tmp_sequence_", suffix = ".fasta", delete = False)
        tmp_fasta = tmp_fasta_file.name
        tmp_prediction_file = tempfile.NamedTemporaryFile(prefix ="mixmhcpred", suffix = ".txt", delete = False)
        tmp_prediction = tmp_prediction_file.name
        seqs = self.generate_nmers(props_dict, [8,9,10,11,12])
        self.generate_fasta(seqs, tmp_fasta)
        alleles = self.get_hla_allels(props_dict, dict_patient_hla)
        self.mixmhcprediction(alleles, tmp_fasta, tmp_prediction)
        pred = self.read_mixmhcpred(tmp_prediction)
        pred_all = self.extract_best_per_pep(pred)
        pred_best = self.extract_best_peptide_per_mutation(pred)
        self.best_peptide = self.add_best_epitope_info(pred_best, "Peptide")
        self.best_score = self.add_best_epitope_info(pred_best, "Score_bestAllele")
        self.best_rank = self.add_best_epitope_info(pred_best, "%Rank_bestAllele")
        self.best_allele = self.add_best_epitope_info(pred_best, "BestAllele")
        self.all_peptides = "|".join(pred_all["Peptide"])
        self.all_scores = "|".join(pred_all["Score_bestAllele"])
        self.all_ranks = "|".join(pred_all["%Rank_bestAllele"])
        self.all_alleles = "|".join(pred_all["BestAllele"])
        # prediction of for wt epitope that correspond to best epitope
        wt = self.extract_WT_for_best(props_dict, self.best_peptide)
        wt_list = [wt]
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix ="tmp_sequence_wt_", suffix = ".fasta", delete = False)
        tmp_fasta = tmp_fasta_file.name
        tmp_prediction_file = tempfile.NamedTemporaryFile(prefix ="mixmhcpred_wt_", suffix = ".txt", delete = False)
        tmp_prediction = tmp_prediction_file.name
        self.generate_fasta(wt_list, tmp_fasta)
        self.mixmhcprediction(alleles, tmp_fasta, tmp_prediction)
        pred_wt = self.read_mixmhcpred(tmp_prediction)
        print >> sys.stderr, pred_wt
        self.best_peptide_wt = self.extract_WT_info(pred_wt, "Peptide")
        score_wt_of_interest = "_".join(["Score",self.best_allele])
        rank_wt_of_interest = "_".join(["%Rank",self.best_allele])
        self.best_score_wt = self.extract_WT_info(pred_wt, score_wt_of_interest)
        self.best_rank_wt = self.extract_WT_info(pred_wt, rank_wt_of_interest)
        # difference in scores between mut and wt
        self.difference_score_mut_wt = self.difference_score(self.best_score,self.best_score_wt)



if __name__ == '__main__':

    import epitope
    import predict_all_epitopes
    from datetime import datetime

    # test with ott data set
    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
    hla_file ="/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
    # test inest data set
    #file = "/flash/projects/WP3/AnFranziska/AnFranziska/head_seqs.txt"
    #hla_file = "/flash/projects/WP3/AnFranziska/AnFranziska/alleles.csv"
    dat = data_import.import_dat_icam(file, False)
    if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
        dat = data_import.change_col_names(dat)
    # available MHC alleles
    set_available_mhc = predict_all_epitopes.Bunchepitopes().add_available_hla_alleles()
    # hla allele of patients
    patient_hlaI = predict_all_epitopes.Bunchepitopes().add_patient_hla_I_allels(hla_file)
    patient_hlaII = predict_all_epitopes.Bunchepitopes().add_patient_hla_II_allels(hla_file)

    print patient_hlaI
    print patient_hlaII

    for ii,i in enumerate(dat[1]):
        if ii < 10:
            print ii
            dict_epi = epitope.Epitope()
            dict_epi.init_properties(dat[0], dat[1][ii])
            prediction = MixMHCpred()
            #print ii
            #print dict_epi.properties

            prediction.main(dict_epi.properties,  patient_hlaI)
            attrs = vars(prediction)
            print attrs


        #def wrapper(func, *args, **kwargs):
        #    def wrapped():
        #        return func(*args, **kwargs)
        #    return wrapped
        #wrapped = wrapper(prediction.main, dict_epi.properties, set_available_mhc, patient_hlaI)
        #prediction.generate_fasta(dict_epi.properties)
        #wrapped = wrapper(prediction.mhc_prediction, set_available_mhc, patient_hlaI)
        #print timeit.timeit(wrapped, number=1)
        #print timeit.timeit(wrapped, number=3)
        #print prediction.mhc_score
        #print prediction.epitope
        #print prediction.allele
        #print prediction.directed_to_TCR
    '''
    startTime1 = datetime.now()
    seqs = []
    col = dat[0].index("X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL.")
    for ii,i in enumerate(dat[1]):
        seqs.append(i[col])
    file_fasta = my_path + "/all_seqs.fasta"
    with open(file_fasta, "w") as f:
        counter = 0
        for seq in seqs:
            f.write(">"+ str(counter) + "\n")
            f.write(seq + "\n")
            counter += 1

    allels_for_prediction = []
    tmp_fasta = file_fasta
    hla_allele = 'HLA-A*23:01,HLA-A*68:02,HLA-B*14:02,HLA-B*49:01,HLA-C*07:01,HLA-C*08:02'
    cmd = "/code/netMHCpan-4.0/netMHCpan -a " + hla_allele + " -f " + tmp_fasta
    p = subprocess.Popen(cmd.split(" "),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    lines = p.stdout
    print cmd
    print lines
    fileout = my_path + "/netmhcpan_out.csv"
    counter = 0
    with open(fileout,"w") as f:
        for line in lines:
            line = line.rstrip().lstrip()
            if line:
                if line.startswith(("#", "-", "HLA", "Prot")):
                    continue
                if counter == 0 and line.startswith("Pos"):
                    counter += 1
                    line = line.split()
                    line = line[0:-1] if len(line) > 13 else line
                    f.write(";".join(line) + "\n")
                    continue
                elif counter >0 and line.startswith("Pos"):
                    continue
                line = line.split()
                line = line[0:-2] if len(line) > 13 else line
                line = ";".join(line)
                f.write(line + "\n")
    endTime1 = datetime.now()
    print >> sys.stderr, "MHC PREDICTION "+ str(startTime1) + "\nend: "+ str(endTime1) + "\nneeded: " + str(endTime1 - startTime1)
    '''










    '''
    # netmhcpan4 for example data set
    for ii,i in enumerate(dat[1]):
        #print dat[1][ii]
        dict_epi = epitope.Epitope()
        dict_epi.init_properties(dat[0], dat[1][ii])
        NetmhcpanBestPrediction.generate_fasta(dict_epi.properties)
        alleles = NetmhcpanBestPrediction.get_hla_allels(dict_epi.properties, hla_file)
        NetmhcpanBestPrediction.mhc_prediction(alleles, list_available_mhc)
        #print mut_position_xmer_seq(dict_epi.properties)
        dict_epi.properties["Position_Xmer_Seq"] = NetmhcpanBestPrediction.mut_position_xmer_seq(dict_epi.properties)
        preds = NetmhcpanBestPrediction.filter_binding_predictions(dict_epi.properties)
        #print preds
        best_epi = NetmhcpanBestPrediction.minimal_binding_score(preds)
        #print best_epi
        NetmhcpanBestPrediction.mutation_in_loop(dict_epi.properties, best_epi)
        dict_epi.properties["%Rank_netmhcpan4"] = NetmhcpanBestPrediction.add_best_epitope_info(best_epi, "%Rank")
        dict_epi.properties["HLA_allele_netmhcpan4"] = NetmhcpanBestPrediction.add_best_epitope_info(best_epi, "HLA")
        dict_epi.properties["Best_epitope_netmhcpan4"] = NetmhcpanBestPrediction.add_best_epitope_info(best_epi, "Icore")
        dict_epi.properties["directed_to_TCR"] = NetmhcpanBestPrediction.mutation_in_loop(dict_epi.properties, best_epi)
        #print dict_epi.properties
        print dict_epi.properties["Best_epitope_netmhcpan4"],dict_epi.properties["HLA_allele_netmhcpan4"], dict_epi.properties["%Rank_netmhcpan4"], dict_epi.properties["directed_to_TCR"]
    '''
