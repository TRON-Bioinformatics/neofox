#!/usr/bin/env python

import subprocess
import sys
import tempfile
from input.helpers import data_import
from input import MHC_I, MHC_II


class NetmhcIIpanBestPrediction:
    def __init__(self):
        self.mhcII_score = "NA"
        self.epitopeII = "NA"
        self.alleleII = "NA"
        self.affinityII = "NA"
        self.affinity_epitopeII = "NA"
        self.affinity_alleleII= "NA"

    def mhc_allele_in_netmhcpan_available(self, allele, set_available_mhc):
        '''checks if mhc prediction is possible for given hla allele
        '''
        return allele in set_available_mhc

    def generate_fasta(self, props, tmpfile, mut = True):
        ''' Writes 27mer to fasta file.
        '''
        if mut == True:
            seq = props["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        elif mut == False:
            seq = props["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        id = ">seq1"
        with open(tmpfile,"w") as f:
            f.write(id + "\n")
            f.write(seq + "\n")

    def get_hla_allels(self, props, hla_patient_dict):
        ''' returns hla allele of patients given in hla_file
        '''
        if "patient.id" in props:
            patientid = props["patient.id"]
            #print >> sys.stderr, patientid
        else:
            try:
                patientid = props["patient"]
            except KeyError:
                patientid = props["patient.x"]
        return hla_patient_dict[patientid]

    def generate_mhcII_alelles_combination_list(self, hla_alleles, set_available_mhc):
        ''' given list of HLA II alleles, returns list of HLA-DRB1 (2x), all possible HLA-DPA1/HLA-DPB1 (4x) and HLA-DQA1/HLA-DPQ1 (4x)
        '''
        allels_for_prediction = []
        dqa_alleles = []
        dpa_alleles = []
        dqb_alleles = []
        dpb_alleles = []
        for allele in hla_alleles:
            if allele.startswith("HLA-DRB1"):
                allele = allele.replace("HLA-","").replace("*", "_").replace(":", "")
                if self.mhc_allele_in_netmhcpan_available(allele, set_available_mhc):
                    allels_for_prediction.append(allele)
            else:
                allele = allele.replace("*", "").replace(":", "")
                if allele.startswith("HLA-DPA"):
                    dpa_alleles.append(allele)
                elif allele.startswith("HLA-DPB"):
                    dpb_alleles.append(allele)
                elif allele.startswith("HLA-DQA"):
                    dqa_alleles.append(allele)
                elif allele.startswith("HLA-DQB"):
                    dqb_alleles.append(allele)
        dp_alleles = [ "-".join([x,y.replace("HLA-","")]) for x in dpa_alleles for y in dpb_alleles]
        dq_alleles = [ "-".join([x,y.replace("HLA-","")]) for x in dqa_alleles for y in dqb_alleles]
        dp_dq_alleles = dp_alleles + dq_alleles
        #print >> sys.stderr, dp_dq_alleles
        for allele in dp_dq_alleles:
            if self.mhc_allele_in_netmhcpan_available(allele, set_available_mhc):
                allels_for_prediction.append(allele)
        return allels_for_prediction

    def mhcII_prediction(self, hla_alleles, set_available_mhc, tmpfasta, tmppred):
        ''' Performs netmhcIIpan prediction for desired hla alleles and writes result to temporary file.
        '''
        allels_for_prediction = self.generate_mhcII_alelles_combination_list(hla_alleles, set_available_mhc)
        hla_allele = ",".join(allels_for_prediction)
        tmp_folder = tempfile.mkdtemp(prefix ="tmp_netmhcIIpan_")
        print(tmp_folder, file=sys.stderr)
        cmd = "/code/net/MHCIIpan/3.2/netMHCIIpan -a " + hla_allele + " -f " + tmpfasta + " -tdir "+ tmp_folder + " -dirty"
        p = subprocess.Popen(cmd, stderr=subprocess.PIPE,stdout=subprocess.PIPE, shell=True)
        lines = p.stdout
        counter = 0
        with open(tmppred, "w") as f:
            for line in lines:
                line = line.decode('utf8').rstrip().lstrip()
                if line:
                    if line.startswith(("#", "-", "Number", "Temporary")):
                        continue
                    if counter == 0 and line.startswith("Seq"):
                        counter += 1
                        line = line.split()
                        line = line[0:-1] if len(line) > 12 else line
                        f.write(";".join(line) + "\n")
                        continue
                    elif counter >0 and line.startswith("Seq"):
                        continue
                    line = line.split()
                    line = line[0:-2] if len(line) > 11 else line
                    f.write(";".join(line) + "\n")

    def mut_position_xmer_seq(self, props):
        '''returns position of mutation in xmer sequence
        '''
        xmer_wt = props["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        xmer_mut = props["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        if len(xmer_wt) == len(xmer_mut):
            p1 = -1
            for i,aa in enumerate(xmer_mut):
                if aa != xmer_wt[i]:
                    p1 = i + 1
        else:
            p1 = 0
            # in case sequences do not have same length
            for a1, a2 in zip(xmer_wt, xmer_mut):
                if a1 == a2:
                    p1 += 1
        return str(p1)

    def epitope_covers_mutation(self, position_mutation, position_epitope, length_epitope):
        '''checks if predicted peptide (15mer) covers mutation
        '''
        cover = False
        if position_mutation != "-1":
            start = int(position_epitope) - 1
            end = start + int(length_epitope) - 1
            if int(position_mutation) >= start and int(position_mutation) <= end:
                cover = True
        return cover

    def filter_binding_predictions(self, props, tmppred):
        '''filters prediction file for predicted epitopes that cover mutations
        '''
        pos_xmer = props["Position_Xmer_Seq"]
        dat_prediction = data_import.import_dat_general(tmppred)
        dat = dat_prediction[1]
        dat_head = dat_prediction[0]
        dat_fil = []
        pos_epi = dat_head.index("Seq")
        epi = dat_head.index("Peptide")
        for ii,i in enumerate(dat):
            if self.epitope_covers_mutation(pos_xmer, i[pos_epi], len(i[epi])):
                dat_fil.append(dat[ii])
        return dat_head, dat_fil

    def minimal_binding_score(self, prediction_tuple, rank = True):
        '''reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used. if rank = False, Aff(nM) is used
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        if rank:
            mhc_sc = dat_head.index("%Rank")
        else:
            mhc_sc = dat_head.index("Affinity(nM)")
        epi = dat_head.index("Peptide")
        hla_allele = dat_head.index("Allele")
        max_score = float(1000000000)
        allele = "NA"
        epitope = "NA"
        row = []
        for ii,i in enumerate(dat):
            mhc_score = float(i[mhc_sc])
            if  mhc_score < max_score:
                max_score = mhc_score
                row = i
        return dat_head, row

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

    def Hamming_check_0_or_1(self, seq1, seq2):
        '''returns number of mismatches between 2 sequences
        '''
        errors = 0
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                errors += 1
                if errors >= 2:
                    return errors
        return errors

    def filter_for_WT_epitope_same_allele(self, prediction_tuple, mut_seq, mut_allele):
        '''returns wt epitope info for given mutated sequence. here best wt for same allele as mutated sequences
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        seq_col = dat_head.index("Peptide")
        allele_col = dat_head.index("Allele")
        wt_epi = "NA"
        for ii,i in enumerate(dat):
            wt_seq = i[seq_col]
            wt_allele = i[allele_col]
            if (len(wt_seq) == len(mut_seq)) and wt_allele == mut_allele:
                numb_mismatch = self.Hamming_check_0_or_1(mut_seq, wt_seq)
                if  numb_mismatch == 1:
                    wt_epi = i
        return(dat_head, wt_epi)

    def filter_for_WT_epitope(self, prediction_tuple, mut_seq):
        '''returns wt epitope info for given mutated sequence. best wt that is allowed to bind to any allele of patient
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        seq_col = dat_head.index("Peptide")
        allele_col = dat_head.index("Allele")
        wt_epi = []
        for ii,i in enumerate(dat):
            wt_seq = i[seq_col]
            wt_allele = i[allele_col]
            if (len(wt_seq) == len(mut_seq)):
                numb_mismatch = self.Hamming_check_0_or_1(mut_seq, wt_seq)
                if  numb_mismatch == 1:
                    wt_epi.append(i)
        dt = (dat_head, wt_epi)
        min = self.minimal_binding_score(dt)
        return(min)

    def main(self, props_dict, set_available_mhc, dict_patient_hla):
        '''Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR
        '''
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix ="tmp_singleseq_", suffix = ".fasta", delete = False)
        tmp_fasta = tmp_fasta_file.name
        tmp_prediction_file = tempfile.NamedTemporaryFile(prefix ="netmhcpanpred_", suffix = ".csv", delete = False)
        tmp_prediction = tmp_prediction_file.name
        #print tmp_prediction
        self.generate_fasta(props_dict, tmp_fasta)
        alleles = self.get_hla_allels(props_dict, dict_patient_hla)
        self.mhcII_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        props_dict["Position_Xmer_Seq"] = self.mut_position_xmer_seq(props_dict)
        preds = self.filter_binding_predictions(props_dict, tmp_prediction)
        best_epi =  self.minimal_binding_score(preds)
        best_epi_affinity =  self.minimal_binding_score(preds, rank = False)

        self.mhcII_score = self.add_best_epitope_info(best_epi, "%Rank")
        self.epitopeII = self.add_best_epitope_info(best_epi, "Peptide")
        self.alleleII = self.add_best_epitope_info(best_epi, "Allele")
        self.affinityII = self.add_best_epitope_info(best_epi_affinity, "Affinity(nM)")
        self.affinity_epitopeII = self.add_best_epitope_info(best_epi_affinity, "Peptide")
        self.affinity_alleleII = self.add_best_epitope_info(best_epi_affinity, "Allele")


if __name__ == '__main__':

    from input import predict_all_epitopes, epitope

    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
    #hla_file ="/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/20190730_alleles.csv"
    hla_file = "/projects/SUMMIT/WP1.2/dataset_annotation/Birmingham/20190821_alleles.csv"
    file = "/projects/SUMMIT/WP1.2/dataset_annotation/Birmingham/in_files/PtBI000048T_1PEB.transcript"
    #file = "/flash/projects/WP3/AnFranziska/AnFranziska/head_seqs.txt"
    #hla_file ="/flash/projects/WP3/AnFranziska/AnFranziska/alleles.csv"
    dat = data_import.import_dat_icam(file, False)
    if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
        dat = data_import.change_col_names(dat)
    if "patient.id" not in dat[0]:
        try:
            patient = file.split("/")[-3]
            if "Pt" not in patient:
                patient = file.split("/")[-1].split(".")[0]
        except IndexError:
            patient = file.split("/")[-1].split(".")[0]
        dat[0].append("patient.id")
        for ii,i in enumerate(dat[1]):
            dat[1][ii].append(str(patient))
    # available MHC alleles
    set_available_mhc = predict_all_epitopes.Bunchepitopes().add_available_hla_alleles()
    set_available_mhcII = predict_all_epitopes.Bunchepitopes().add_available_hla_alleles(mhc=MHC_II)
    #print >> sys.stderr, set_available_mhcII
    # hla allele of patients
    patient_hlaI = predict_all_epitopes.Bunchepitopes().add_patient_hla_I_allels(hla_file)
    patient_hlaII = predict_all_epitopes.Bunchepitopes().add_patient_hla_II_allels(hla_file)

    #print patient_hlaI
    #print patient_hlaII

    for ii,i in enumerate(dat[1]):
        if ii < 10:
            #print ii
            dict_epi = epitope.Epitope()
            dict_epi.init_properties(dat[0], dat[1][ii])
            prediction = NetmhcIIpanBestPrediction()
            #print ii
            #print dict_epi.properties

            prediction.main(dict_epi.properties, set_available_mhcII, patient_hlaII)
            attrs = vars(prediction)
            #print attrs