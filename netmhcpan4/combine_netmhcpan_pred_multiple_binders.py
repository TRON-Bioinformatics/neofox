#!/usr/bin/env python

import os
import sys
import tempfile

my_path = os.path.abspath(os.path.dirname(__file__))
my_path2 = "/".join(my_path.split("/")[0:-1])
sys.path.insert(0, my_path2)
sys.path.insert(0, my_path)


import netmhcpan_prediction
import multiple_binders


class Bestandmultiplebinder:
    def __init__(self):
        self.mean_type = ["arithmetic", "harmonic", "geometric"]
        self.MHC_score_all_epitopes = []
        self.MHC_score_top10 = []
        self.MHC_score_best_per_alelle = []
        self.MHC_number_strong_binders = ""
        self.MHC_number_weak_binders = ""
        self.MHC_epitope_seqs = ""
        self.MHC_epitope_scores = ""
        self.MHC_epitope_alleles = ""
        self.best4_mhc_score = "NA"
        self.best4_mhc_epitope = "NA"
        self.best4_mhc_allele = "NA"
        self.directed_to_TCR = "NA"
        self.best4_affinity = "NA"
        self.best4_affinity_epitope = "NA"
        self.best4_affinity_allele = "NA"
        self.best4_affinity_directed_to_TCR = "NA"
        self.epitope_affinities = ""
        self.generator_rate = ""
        # WT features
        self.MHC_epitope_scores_WT = ""
        self.MHC_epitope_seqs_WT = ""
        self.MHC_epitope_alleles_WT = []
        self.MHC_score_top10_WT = []
        self.MHC_score_all_epitopes_WT = []
        self.MHC_score_best_per_alelle_WT = []
        self.MHC_number_strong_binders_WT = ""
        self.MHC_number_weak_binders_WT = ""
        self.best4_mhc_score_WT = "NA"
        self.best4_mhc_epitope_WT = "NA"
        self.best4_mhc_allele_WT = "NA"
        self.best4_affinity_WT = "NA"
        self.best4_affinity_epitope_WT = "NA"
        self.best4_affinity_allele_WT = "NA"
        self.epitope_affinities_WT = ""
        self.generator_rate_WT = ""




    def main(self, epi_dict, patient_hlaI, set_available_mhc):
        '''predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        '''
        ### PREDICTION FOR MUTATED SEQUENCE
        xmer_mut = epi_dict["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        print >> sys.stderr, "MUT seq: " + xmer_mut
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix ="tmp_singleseq_", suffix = ".fasta", delete = False)
        tmp_fasta = tmp_fasta_file.name
        tmp_prediction_file = tempfile.NamedTemporaryFile(prefix ="netmhcpanpred_", suffix = ".csv", delete = False)
        tmp_prediction = tmp_prediction_file.name
        np = netmhcpan_prediction.NetmhcpanBestPrediction()
        mb = multiple_binders.MultipleBinding()
        np.generate_fasta(epi_dict, tmp_fasta, mut = True)
        alleles = np.get_hla_allels(epi_dict, patient_hlaI)
        np.mhc_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        epi_dict["Position_Xmer_Seq"] = np.mut_position_xmer_seq(epi_dict)
        preds = np.filter_binding_predictions(epi_dict, tmp_prediction)
        #print >> sys.stderr, preds
        # multiple binding
        list_tups = mb.generate_epi_tuple(preds)
        self.MHC_epitope_scores = "/".join([tup[0] for tup in list_tups])
        self.epitope_affinities = "/".join([tup[1] for tup in list_tups])
        self.MHC_epitope_seqs = "/".join([tup[2] for tup in list_tups])
        self.MHC_epitope_alleles = "/".join([tup[3] for tup in list_tups])
        top10 = mb.extract_top10_epis(list_tups)
        best_per_alelle = mb.extract_best_epi_per_alelle(list_tups, alleles)
        all = mb.scores_to_list(list_tups)
        all_affinities = mb.affinities_to_list(list_tups)
        top10 = mb.scores_to_list(top10)
        self.MHC_score_top10 = mb.wrapper_mean_calculation(top10)
        best_per_alelle = mb.scores_to_list(best_per_alelle)
        self.MHC_score_all_epitopes = mb.wrapper_mean_calculation(all)
        self.MHC_score_best_per_alelle = mb.wrapper_mean_calculation(best_per_alelle)
        self.MHC_number_strong_binders = mb.determine_number_of_binders(all, 1)
        self.MHC_number_weak_binders = mb.determine_number_of_binders(all, 2)
        # best prediction
        best_epi =  np.minimal_binding_score(preds)
        self.best4_mhc_score =np.add_best_epitope_info(best_epi, "%Rank")
        self.best4_mhc_epitope = np.add_best_epitope_info(best_epi, "Icore")
        self.best4_mhc_allele = np.add_best_epitope_info(best_epi, "HLA")
        self.directed_to_TCR = np.mutation_in_loop(epi_dict, best_epi)
        best_epi_affinity =  np.minimal_binding_score(preds, rank = False)
        self.best4_affinity = np.add_best_epitope_info(best_epi_affinity, "Aff(nM)")
        self.best4_affinity_epitope = np.add_best_epitope_info(best_epi_affinity, "Icore")
        self.best4_affinity_allele = np.add_best_epitope_info(best_epi_affinity, "HLA")
        self.best4_affinity_directed_to_TCR =  np.mutation_in_loop(epi_dict, best_epi_affinity)
        # multiple binding based on affinity
        self.generator_rate = mb.determine_number_of_binders(list_scores = all_affinities, threshold = 50)

        ### PREDICTION FOR WT SEQUENCE
        xmer_wt = epi_dict["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        #print >> sys.stderr, "WT seq: " + xmer_wt
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix ="tmp_singleseq_", suffix = ".fasta", delete = False)
        tmp_fasta = tmp_fasta_file.name
        tmp_prediction_file = tempfile.NamedTemporaryFile(prefix ="netmhcpanpred_", suffix = ".csv", delete = False)
        tmp_prediction = tmp_prediction_file.name
        np = netmhcpan_prediction.NetmhcpanBestPrediction()
        mb = multiple_binders.MultipleBinding()
        np.generate_fasta(epi_dict, tmp_fasta, mut = False)
        np.mhc_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        preds = np.filter_binding_predictions(epi_dict, tmp_prediction)
        # multiple binding
        list_tups = mb.generate_epi_tuple(preds)
        self.MHC_epitope_scores_WT = "/".join([tup[0] for tup in list_tups])
        self.epitope_affinities_WT = "/".join([tup[1] for tup in list_tups])
        self.MHC_epitope_seqs_WT = "/".join([tup[2] for tup in list_tups])
        self.MHC_epitope_alleles_WT = "/".join([tup[3] for tup in list_tups])
        top10 = mb.extract_top10_epis(list_tups)
        best_per_alelle = mb.extract_best_epi_per_alelle(list_tups, alleles)
        all = mb.scores_to_list(list_tups)
        all_affinities = mb.affinities_to_list(list_tups)
        top10 = mb.scores_to_list(top10)
        self.MHC_score_top10_WT = mb.wrapper_mean_calculation(top10)
        best_per_alelle = mb.scores_to_list(best_per_alelle)
        self.MHC_score_all_epitopes_WT = mb.wrapper_mean_calculation(all)
        self.MHC_score_best_per_alelle_WT = mb.wrapper_mean_calculation(best_per_alelle)
        self.MHC_number_strong_binders_WT = mb.determine_number_of_binders(all, 1)
        self.MHC_number_weak_binders_WT = mb.determine_number_of_binders(all, 2)
        # best prediction
        best_epi =  np.minimal_binding_score(preds)
        self.best4_mhc_score_WT =np.add_best_epitope_info(best_epi, "%Rank")
        self.best4_mhc_epitope_WT = np.add_best_epitope_info(best_epi, "Icore")
        self.best4_mhc_allele_WT = np.add_best_epitope_info(best_epi, "HLA")
        best_epi_affinity =  np.minimal_binding_score(preds, rank = False)
        self.best4_affinity_WT =np.add_best_epitope_info(best_epi_affinity, "Aff(nM)")
        self.best4_affinity_epitope_WT = np.add_best_epitope_info(best_epi_affinity, "Icore")
        self.best4_affinity_allele_WT = np.add_best_epitope_info(best_epi_affinity, "HLA")
        self.generator_rate_WT = mb.determine_number_of_binders(list_scores = all_affinities, threshold = 50)
        print >> sys.stderr, "WT: " + self.generator_rate_WT +"; MUT: "+ self.generator_rate




if __name__ == '__main__':

    import epitope
    import predict_all_epitopes
    from helpers import data_import

    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_SD.csv"
    #file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/riaz/output_tables_pre/test.txt"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/20170713_IS_IM_data.complete.update_Dv10.csv.annotation.csv_v2.csv"
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_fulldat.txt"
    #hla_file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/indels/RB_0004_labHLA_V2.csv"
    #hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/riaz/output_tables_pre/alleles.csv"
    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
    hla_file ="/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
    dat = data_import.import_dat_icam(file, False)
    if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
        dat = data_import.change_col_names(dat)
    # available MHC alleles
    set_available_mhc = predict_all_epitopes.Bunchepitopes().add_available_hla_alleles()
    # hla allele of patients
    patient_hlaI = predict_all_epitopes.Bunchepitopes().add_patient_hla_I_allels(hla_file)
    patient_hlaII = predict_all_epitopes.Bunchepitopes().add_patient_hla_II_allels(hla_file)


    Allepit = {}
    for ii,i in enumerate(dat[1]):
        if ii < 10:
            dict_epi = epitope.Epitope()
            dict_epi.init_properties(dat[0], dat[1][ii])
            x =  Bestandmultiplebinder()
            x.main(dict_epi.properties, patient_hlaI, set_available_mhc)
            #print x.MHC_epitope_scores_WT
            #print x.MHC_epitope_seqs_WT
            #print x.MHC_epitope_seqs
            attrs = vars(x)
            print attrs
            '''
            for sc, mn in zip(x.MHC_score_all_epitopes, x.mean_type):
                dict_epi.add_features(sc, "MB_score_all_epitopes_" + mn)
            for sc, mn in zip(x.MHC_score_top10, x.mean_type):
                dict_epi.add_features(sc, "MB_score_top10_" + mn)
            for sc, mn in zip(x.MHC_score_best_per_alelle, x.mean_type):
                dict_epi.add_features(sc, "MB_score_best_per_alelle_" + mn)
            dict_epi.add_features(x.MHC_epitope_scores, "MB_epitope_scores")
            dict_epi.add_features(x.MHC_epitope_seqs, "MB_epitope_sequences")
            dict_epi.add_features(x.MHC_epitope_alleles, "MB_alleles")
            dict_epi.add_features(x.MHC_number_strong_binders, "MB_number_of_strong_binders")
            dict_epi.add_features(x.MHC_number_weak_binders, "MB_number_of_weak_binders")
            dict_epi.add_features(x.best4_mhc_score, "best4_mhc_score")
            dict_epi.add_features(x.best4_mhc_epitope, "best4_mhc_epitope")
            dict_epi.add_features(x.best4_mhc_allele, "best4_mhc_allele")
            dict_epi.add_features(x.directed_to_TCR, "directed_to_TCR")
            z = dict_epi.properties
            for key in z:
                if key not in Allepit:
                    Allepit[key] = [z[key]]
                else:
                    Allepit[key].append(z[key])
            '''
    #predict_all_epitopes.Bunchepitopes().write_to_file(Allepit)
