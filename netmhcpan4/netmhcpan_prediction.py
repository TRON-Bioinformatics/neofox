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


class NetmhcpanBestPrediction:
    def __init__(self):
        self.mhc_score = "NA"
        self.epitope = "NA"
        self.allele = "NA"
        self.directed_to_TCR = "NA"
        self.affinity = "NA"
        self.affinity_epitope = "NA"
        self.affinity_allele = "NA"
        self.affinity_directed_to_TCR = "NA"

    def mhc_allele_in_netmhcpan_available(self, allele, set_available_mhc):
        '''checks if mhc prediction is possible for given hla allele
        '''
        return allele in set_available_mhc

    def generate_fasta(self, props, tmpfile, mut = True):
        ''' Writes 27mer to fasta file.
        '''
        #fastafile = my_path + "/tmp.fasta"
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
        else:
            try:
                patientid = props["patient"]
            except KeyError:
                patientid = props["patient.x"]
        return hla_patient_dict[patientid]

    def mhc_prediction(self, hla_alleles, set_available_mhc, tmpfasta, tmppred):
        ''' Performs netmhcpan4 prediction for desired hla allele and writes result to temporary file.
        '''
        allels_for_prediction = []
        #tmp_fasta = "/".join([my_path, "tmp.fasta"])
        for allele in hla_alleles:
            allele = allele.replace("*", "")
            if self.mhc_allele_in_netmhcpan_available(allele, set_available_mhc):
                allels_for_prediction.append(allele)
        hla_allele = ",".join(allels_for_prediction)
        cmd = "/code/netMHCpan-4.0/netMHCpan -a " + hla_allele + " -f " + tmpfasta + " -BA"
        p = subprocess.Popen(cmd.split(" "),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        lines = lines = p.stdout
        #print >> sys.stderr, lines
        #stdoutdata, stderrdata = p.communicate()
        #lines = stdoutdata
        #print >> sys.stderr, "NETMHCPAN ERROR: ",stderrdata
        #fileout = my_path + "/netmhcpan_out.csv"
        counter = 0
        with open(tmppred,"w") as f:
            for line in lines:
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("#", "-", "HLA", "Prot")):
                        continue
                    if counter == 0 and line.startswith("Pos"):
                        counter += 1
                        line = line.split()
                        line = line[0:-1] if len(line) > 14 else line
                        f.write(";".join(line) + "\n")
                        continue
                    elif counter >0 and line.startswith("Pos"):
                        continue
                    line = line.split()
                    line = line[0:-2] if len(line) > 14 else line
                    line = ";".join(line)
                    f.write(line + "\n")
        #os.remove(tmp_fasta)

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

    def epitope_covers_mutation(self, position_mutation, position_epitope, length_epitope):
        '''checks if predicted epitope covers mutation
        '''
        cover = False
        if position_mutation != "-1":
            start = int(position_epitope)
            end = start + int(length_epitope) - 1
            if int(position_mutation) >= start and int(position_mutation) <= end:
                cover = True
        return cover

    def filter_binding_predictions(self, props, tmppred):
        '''filters prediction file for predicted epitopes that cover mutations
        '''
        pos_xmer = props["Position_Xmer_Seq"]
        #prediction_file = "/".join([my_path, "netmhcpan_out.csv"])
        dat_prediction = data_import.import_dat_general(tmppred)
        print >> sys.stderr, dat_prediction
        #print >> sys.stderr, dat_prediction
        #os.remove(prediction_file)
        dat = dat_prediction[1]
        dat_head = dat_prediction[0]
        dat_fil = []
        pos_epi = dat_head.index("Pos")
        epi = dat_head.index("Icore")
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
            mhc_sc = dat_head.index("Aff(nM)")
        print >> sys.stderr, mhc_sc
        epi = dat_head.index("Icore")
        hla_allele = dat_head.index("HLA")
        max_score = float(999)
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

    def mutation_in_loop(self, props, epitope_tuple):
        '''returns if mutation is directed to TCR (yes or no)
        '''
        pos_xmer = props["Position_Xmer_Seq"]
        dat_head = epitope_tuple[0]
        dat_epi = epitope_tuple[1]
        pos_epi = dat_head.index("Pos")
        del_pos = dat_head.index("Gp")
        del_len = dat_head.index("Gl")
        directed_to_TCR = "no"
        try:
            if del_pos > 0:
                pos = int(dat_epi[pos_epi])
                start = pos + int(dat_epi[del_pos]) - 1
                end = start + int(dat_epi[del_len])
                if int(pos_xmer) > start and int(pos_xmer) <= end:
                    directed_to_TCR = "yes"
            return directed_to_TCR
        except IndexError:
            return "NA"

    def main(self, props_dict, set_available_mhc, dict_patient_hla):
        '''Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR
        '''
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix ="tmp_singleseq_", suffix = ".fasta", delete = False)
        tmp_fasta = tmp_fasta_file.name
        tmp_prediction_file = tempfile.NamedTemporaryFile(prefix ="netmhcpanpred_", suffix = ".csv", delete = False)
        tmp_prediction = tmp_prediction_file.name
        self.generate_fasta(props_dict, tmp_fasta)
        alleles = self.get_hla_allels(props_dict, dict_patient_hla)
        self.mhc_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        props_dict["Position_Xmer_Seq"] = self.mut_position_xmer_seq(props_dict)
        preds = self.filter_binding_predictions(props_dict, tmp_prediction)
        best_epi =  self.minimal_binding_score(preds)
        best_epi_affinity =  self.minimal_binding_score(preds, rank = False)
        self.mhc_score = self.add_best_epitope_info(best_epi, "%Rank")
        self.epitope = self.add_best_epitope_info(best_epi, "Icore")
        self.allele = self.add_best_epitope_info(best_epi, "HLA")
        self.directed_to_TCR = self.mutation_in_loop(props_dict, best_epi)
        self.affinity = self.add_best_epitope_info(best_epi_affinity, "Aff(nM)")
        self.affinity_epitope = self.add_best_epitope_info(best_epi_affinity, "Icore")
        self.affinity_allele = self.add_best_epitope_info(best_epi_affinity, "HLA")
        self.affinity_directed_to_TCR =  self.mutation_in_loop(props_dict, best_epi_affinity)


if __name__ == '__main__':

    import epitope
    import predict_all_epitopes
    from datetime import datetime

    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_SD.csv"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/test_hugo.txt"
    #file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/riaz/output_tables_pre/test.txt"
    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_fulldat.txt"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_fulldat.txt"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/Datasets/201808_ivac_fullData_New_Analysis_correct_mergedepit_num_loop_features.csv"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
    #hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
    #hla_file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/indels/RB_0004_labHLA_V2.csv"
    #hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/hugo/output_tables/alleles.csv"
    #hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/riaz/output_tables_pre/alleles.csv"
    hla_file ="/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
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
            prediction = NetmhcpanBestPrediction()
            #print ii
            #print dict_epi.properties

            prediction.main(dict_epi.properties, set_available_mhc, patient_hlaI)
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
