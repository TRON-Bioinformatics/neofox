#!/usr/bin/env python

import os
import sys
import math

my_path = os.path.abspath(os.path.dirname(__file__))
my_path2 = "/".join(my_path.split("/")[0:-1])
sys.path.insert(0, my_path2)
sys.path.insert(0, my_path)


import netmhcpan_prediction

class MultipleBinding:
    def __init__(self):
        self.mean_type = ["arithmetic", "harmonic", "geometric"]
        self.score_all_epitopes = []
        self.score_top10 = []
        self.score_best_per_alelle = []
        self.number_strong_binders = ""
        self.number_weak_binders = ""
        self.epitope_seqs = ""
        self.epitope_scores = ""
        self.epitope_alleles = ""
        self.epitope_affinities = ""
        self.generator_rate = ""

    def calc_arimetric_mean(self, list_numbers):
        '''Calculates the arithmetic mean from a list of numbers
        '''
        sm = 0
        for num in list_numbers:
            sm = sm + float(num)
        try:
            return str(float(sm/len(list_numbers)))
        except ZeroDivisionError:
            return "NA"

    def calc_harmonic_mean(self, list_numbers):
        '''Calculates the harmonic mean from a list of numbers
        '''
        nums = [float(num)**-1 for num in list_numbers]
        sm = 0
        for num in nums:
            sm = sm + num
        try:
            return str(float(len(nums)/sm))
        except ZeroDivisionError:
            return "NA"

    def calc_geometric_mean_inefficient(self, list_numbers):
        '''Calculates the geometric mean from a list of numbers
        '''
        pr = 1
        for num in list_numbers:
            pr = float(num) * pr
        return str(pr**(len(list_numbers)**-1))

    def calc_geometric_mean(self, list_numbers):
        '''Calculates the geometric mean from a list of numbers; avoids product --> suitable for larger list of number
        '''
        sm = 0
        for num in list_numbers:
            num_log = math.log(float(num))
            sm = float(num_log) + sm
        try:
            num_log_mean = sm/len(list_numbers)
            num_log_mean_exp = math.exp(num_log_mean)
            return str(num_log_mean_exp)
        except ZeroDivisionError:
            return "NA"



    def wrapper_mean_calculation(self, list_numbers):
        '''returns list of arithmetic, harmonic and geometric mean from a list of numbers
        '''
        return [self.calc_arimetric_mean(list_numbers), self.calc_harmonic_mean(list_numbers), self.calc_geometric_mean(list_numbers)]

    def generate_epi_tuple(self, prediction_out, mhc = "mhcI"):
        '''Takes netmhcpan4 output or netmhcpanII output as input (parsed with Netmhc[II]panBestPrediction().filter_binding_predictions) and
        returns tuple of mhc binding rank scores, epitope and HLA allele for all predicted epitopes as list
        '''
        pred_data = prediction_out[1]
        list_of_tuples = []
        for ii,i in enumerate(pred_data):
            if mhc == "mhcII":
                # rank, affinity, epitope sequence, allele
                list_of_tuples.append((i[-2],i[-3], i[2], i[1]))
            else:
                # rank, affinity, epitope sequence, allele
                list_of_tuples.append((i[-1],i[-2], i[-5], i[1]))
        return list_of_tuples

    def extract_top10_epis(self, tuple_epis):
        '''this function sorts the predicted epitopes based on the mhc rank score and returns the top10 with lowest mhc binding score
        '''
        tuple_epis.sort(key=lambda x: float(x[0]))
        return tuple_epis[0:9]

    def check_for_homozygosity(self, patient_alleles):
        """ returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles. Otherwise retunrs empty list
        """
        homozygos_alleles = []
        [homozygos_alleles.append(allele) for allele in patient_alleles if patient_alleles.count(allele) > 1]
        return homozygos_alleles


    def extract_best_epi_per_alelle(self, tuple_epis, alleles):
        '''this function returns the predicted epitope with the lowest binding score for each patient allele
        '''
        homo_alleles = self.check_for_homozygosity(alleles)
        dict_allels = {}
        for allele in alleles:
            for epi in tuple_epis:
                if allele == epi[-1]:
                    if allele not in dict_allels:
                        dict_allels[allele] = [epi]
                    else:
                        dict_allels[allele].append(epi)
        best_epis_per_allele = []
        for allele in dict_allels:
            dict_allels[allele].sort(key=lambda x: float(x[0]))
            best_epis_per_allele.append(dict_allels[allele][0])
            if allele in homo_alleles:
                # append homozygous allleles two times
                best_epis_per_allele.append(dict_allels[allele][0])
        return best_epis_per_allele

    def scores_to_list(self, tuple_epis):
        '''Takes list of epitope tuple as input and returns a list of mhc rank scores of these tuples
        '''
        list_score = []
        for epi in tuple_epis:
            list_score.append(epi[0])
        return list_score

    def affinities_to_list(self, tuple_epis):
        '''Takes list of epitope tuple as input and returns a list of mhc rank scores of these tuples
        '''
        list_score = []
        for epi in tuple_epis:
            list_score.append(epi[1])
        return list_score

    def determine_number_of_binders(self, list_scores, threshold = 2):
        '''Determines the number of HLA binders per mutation based on a threshold. Default is set to 2, which is threshold for weak binding using netmhcpan4.
        '''
        number_binders = 0
        for score in list_scores:
            if float(score) < threshold:
                number_binders += 1
        return str(number_binders)

    def main(self, epi_dict, alleles, set_available_mhc):
        '''takes epitope dictionary as input and returns several scores that describe multiple binding.
        '''
        netmhcpan_prediction.NetmhcpanBestPrediction().generate_fasta(epi_dict)
        alleles = netmhcpan_prediction.NetmhcpanBestPrediction().get_hla_allels(epi_dict, patient_hlaI)
        netmhcpan_prediction.NetmhcpanBestPrediction().mhc_prediction(alleles, set_available_mhc)
        epi_dict["Position_Xmer_Seq"] = netmhcpan_prediction.NetmhcpanBestPrediction().mut_position_xmer_seq(epi_dict)
        preds = netmhcpan_prediction.NetmhcpanBestPrediction().filter_binding_predictions(epi_dict)
        list_tups = MultipleBinding().generate_epi_tuple(preds)
        self.epitope_scores = "/".join([tup[0] for tup in list_tups])
        self.epitope_affinities = "/".join([tup[1] for tup in list_tups])
        self.epitope_seqs = "/".join([tup[2] for tup in list_tups])
        self.epitope_alleles = "/".join([tup[3] for tup in list_tups])
        top10 = self.extract_top10_epis(list_tups)
        #print top10
        best_per_alelle = self.extract_best_epi_per_alelle(list_tups, alleles)
        all = self.scores_to_list(list_tups)
        all_affinities = self.affinities_to_list(list_tups)
        top10 = self.scores_to_list(top10)
        #print top10
        self.score_top10 = self.wrapper_mean_calculation(top10)
        best_per_alelle = self.scores_to_list(best_per_alelle)
        self.score_all_epitopes = self.wrapper_mean_calculation(all)
        self.score_best_per_alelle = self.wrapper_mean_calculation(best_per_alelle)
        self.number_strong_binders = self.determine_number_of_binders(all, 0.5)
        self.number_weak_binders = self.determine_number_of_binders(all, 2)
        self.generator_rate = self.determine_number_of_binders(list_scores = all_affinities, threshold = 50)


if __name__ == '__main__':

    import epitope
    import predict_all_epitopes
    from helpers import data_import
    from netmhcpan_prediction import NetmhcpanBestPrediction


    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_SD.csv"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/20170713_IS_IM_data.complete.update_Dv10.csv.annotation.csv_v2.csv"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_fulldat.txt"
    hla_file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/indels/RB_0004_labHLA_V2.csv"
    dat = data_import.import_dat_icam(file, False)
    # available MHC alleles
    set_available_mhc = predict_all_epitopes.Bunchepitopes().add_available_hla_alleles()

    # hla allele of patients
    patient_hlaI = predict_all_epitopes.Bunchepitopes().add_patient_hla_I_allels(hla_file)
    patient_hlaII = predict_all_epitopes.Bunchepitopes().add_patient_hla_II_allels(hla_file)

    Allepit = {}
    for ii,i in enumerate(dat[1]):
        #print ii
        dict_epi = epitope.Epitope()
        dict_epi.init_properties(dat[0], dat[1][ii])
        x =  MultipleBinding()
        x.main(dict_epi.properties, patient_hlaI, set_available_mhc)
        for sc, mn in zip(x.score_all_epitopes, x.mean_type):
            dict_epi.add_features(sc, "MB_score_all_epitopes_" + mn)
        for sc, mn in zip(x.score_top10, x.mean_type):
            dict_epi.add_features(sc, "MB_score_top10_" + mn)
        for sc, mn in zip(x.score_best_per_alelle, x.mean_type):
            dict_epi.add_features(sc, "MB_score_best_per_alelle_" + mn)
        dict_epi.add_features(x.epitope_scores, "MB_epitope_scores")
        dict_epi.add_features(x.epitope_seqs, "MB_epitope_sequences")
        dict_epi.add_features(x.epitope_alleles, "MB_alleles")
        dict_epi.add_features(x.number_strong_binders, "MB_number_of_strong_binders")
        dict_epi.add_features(x.number_weak_binders, "MB_number_of_weak_binders")
        z = dict_epi.properties
        for key in z:
            if key not in Allepit:
                Allepit[key] = [z[key]]
            else:
                Allepit[key].append(z[key])
    predict_all_epitopes.Bunchepitopes().write_to_file(Allepit)
