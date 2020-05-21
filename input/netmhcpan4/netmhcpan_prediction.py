#!/usr/bin/env python

from input.helpers import data_import, properties_manager, intermediate_files


class NetMhcPanBestPrediction:

    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.mhc_score = "NA"
        self.epitope = "NA"
        self.allele = "NA"
        self.directed_to_TCR = "NA"
        self.affinity = "NA"
        self.affinity_epitope = "NA"
        self.affinity_allele = "NA"
        self.affinity_directed_to_TCR = "NA"
        self.mhcI_score_9mer = "NA"
        self.mhcI_score_allele_9mer = "NA"
        self.mhcI_score_epitope_9mer = "NA"
        self.mhcI_affinity_9mer = "NA"
        self.mhcI_affinity_allele_9mer = "NA"
        self.mhcI_affinity_epitope_9mer = "NA"

    def mhc_allele_in_netmhcpan_available(self, allele, set_available_mhc):
        '''checks if mhc prediction is possible for given hla allele
        '''
        return allele in set_available_mhc

    def mhc_prediction(self, hla_alleles, set_available_mhc, tmpfasta, tmppred):
        ''' Performs netmhcpan4 prediction for desired hla allele and writes result to temporary file.
        '''
        allels_for_prediction = []
        for allele in hla_alleles:
            allele = allele.replace("*", "")
            if self.mhc_allele_in_netmhcpan_available(allele, set_available_mhc):
                allels_for_prediction.append(allele)
        hla_allele = ",".join(allels_for_prediction)
        cmd = [
            self.configuration.net_mhc_pan,
            "-a", hla_allele,
            "-f", tmpfasta,
            "-BA"]
        lines, _ = self.runner.run_command(cmd)
        counter = 0
        with open(tmppred, "w") as f:
            for line in lines.splitlines():
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
                    elif counter > 0 and line.startswith("Pos"):
                        continue
                    line = line.split()
                    line = line[0:-2] if len(line) > 14 else line
                    line = ";".join(line)
                    f.write(line + "\n")

    def mut_position_xmer_seq(self, xmer_wt, xmer_mut):
        '''returns position of mutation in xmer sequence
        '''
        if len(xmer_wt) == len(xmer_mut):
            p1 = -1
            for i, aa in enumerate(xmer_mut):
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
        '''checks if predicted epitope covers mutation
        '''
        cover = False
        if position_mutation != "-1":
            start = int(position_epitope)
            end = start + int(length_epitope) - 1
            if int(position_mutation) >= start and int(position_mutation) <= end:
                cover = True
        return cover

    def filter_binding_predictions(self, position_xmer, tmppred):
        '''filters prediction file for predicted epitopes that cover mutations
        '''
        dat_prediction = data_import.import_dat_general(tmppred)
        dat = dat_prediction[1]
        dat_head = dat_prediction[0]
        dat_fil = []
        pos_epi = dat_head.index("Pos")
        epi = dat_head.index("Peptide")
        for ii, i in enumerate(dat):
            if self.epitope_covers_mutation(position_xmer, i[pos_epi], len(i[epi])):
                dat_fil.append(dat[ii])
        return dat_head, dat_fil

    def minimal_binding_score(self, prediction_tuple, rank=True):
        '''reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used. if rank = False, Aff(nM) is used
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        if rank:
            mhc_sc = dat_head.index("%Rank")
        else:
            mhc_sc = dat_head.index("Aff(nM)")
        max_score = float(1000000000000)
        row = []
        for ii, i in enumerate(dat):
            mhc_score = float(i[mhc_sc])
            if mhc_score < max_score:
                max_score = mhc_score
                row = i
        return dat_head, row

    def add_best_epitope_info(self, epitope_tuple, column_name):
        '''returns desired information of prediction of best epitope from netmhcpan output;
        e.g. "%Rank": MHC I score, "HLA": HLA allele, "Peptide": best epitope
        '''
        dat_head = epitope_tuple[0]
        dat = epitope_tuple[1]
        val = dat_head.index(column_name)
        try:
            return dat[val]
        except IndexError:
            return "NA"

    def mutation_in_loop(self, position_xmer, epitope_tuple):
        """
        returns if mutation is directed to TCR (yes or no)
        """
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
                if start < int(position_xmer) <= end:
                    directed_to_TCR = "yes"
            return directed_to_TCR
        except IndexError:
            return "NA"

    def filter_for_9mers(self, prediction_tuple):
        '''returns only predicted 9mers
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        seq_col = dat_head.index("Peptide")
        dat_9mers = []
        for ii, i in enumerate(dat):
            seq = i[seq_col]
            if len(seq) == 9:
                dat_9mers.append(i)
        return dat_head, dat_9mers

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

    def filter_for_WT_epitope(self, prediction_tuple, mut_seq, mut_allele):
        '''returns wt epitope info for given mutated sequence. best wt that is allowed to bind to any allele of patient
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        seq_col = dat_head.index("Peptide")
        allele_col = dat_head.index("HLA")
        wt_epi = []
        for ii, i in enumerate(dat):
            wt_seq = i[seq_col]
            wt_allele = i[allele_col]
            if (len(wt_seq) == len(mut_seq)):
                numb_mismatch = self.Hamming_check_0_or_1(mut_seq, wt_seq)
                if numb_mismatch == 1:
                    wt_epi.append(i)
        dt = (dat_head, wt_epi)
        min = self.minimal_binding_score(dt)
        return (min)


# if __name__ == '__main__':
#
#     from input import predict_all_epitopes, epitope
#
#     # test with ott data set
#     file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
#     hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
#     # test inest data set
#     # file = "/flash/projects/WP3/AnFranziska/AnFranziska/head_seqs.txt"
#     # hla_file = "/flash/projects/WP3/AnFranziska/AnFranziska/alleles.csv"
#     dat = data_import.import_dat_icam(file, False)
#     if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
#         dat = data_import.change_col_names(dat)
#     # available MHC alleles
#     set_available_mhc = predict_all_epitopes.Bunchepitopes().load_available_hla_alleles()
#     # hla allele of patients
#     patient_hlaI = predict_all_epitopes.Bunchepitopes().load_patient_hla_I_allels(hla_file)
#     patient_hlaII = predict_all_epitopes.Bunchepitopes().load_patient_hla_II_allels(hla_file)
#
#     print(patient_hlaI)
#     print(patient_hlaII)
#
#     for ii, i in enumerate(dat[1]):
#         if ii < 10:
#             print(ii)
#             dict_epi = epitope.Epitope()
#             dict_epi.init_properties(dat[0], dat[1][ii])
#             prediction = NetMhcPanBestPrediction()
#             prediction.main(dict_epi.properties, set_available_mhc, patient_hlaI)
#             attrs = vars(prediction)
#             print(attrs)
