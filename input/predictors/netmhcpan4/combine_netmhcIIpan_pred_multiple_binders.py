#!/usr/bin/env python

from logzero import logger

import input.predictors.netmhcpan4.netmhcIIpan_prediction as netmhcIIpan_prediction
from input import MHC_II
from input.helpers import intermediate_files
from input.predictors.netmhcpan4 import multiple_binders


class BestAndMultipleBinderMhcII:

    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.mean_type = ["arithmetic", "harmonic", "geometric"]
        self.MHCII_score_all_epitopes = ["NA", "NA", "NA"]
        self.MHCII_score_top10 = ["NA", "NA", "NA"]
        self.MHCII_score_best_per_alelle = ["NA", "NA", "NA"]
        self.MHCII_number_strong_binders = "NA"
        self.MHCII_number_weak_binders = "NA"
        self.MHCII_epitope_seqs = "NA"
        self.MHCII_epitope_scores = "NA"
        self.MHCII_epitope_alleles = "NA"
        self.best_mhcII_pan_score = "NA"
        self.best_mhcII_pan_epitope = "-"
        self.best_mhcII_pan_allele = "NA"
        self.best_mhcII_pan_position = "NA"
        self.best_mhcII_pan_affinity = "NA"
        self.best_mhcII_pan_affinity_epitope = "-"
        self.best_mhcII_pan_affinity_allele = "NA"
        self.best_mhcII_pan_affinity_position = "NA"
        # WT features
        self.MHCII_epitope_scores_WT = "Na"
        self.MHCII_epitope_seqs_WT = "NA"
        self.MHCII_epitope_alleles_WT = "NA"
        self.MHCII_score_top10_WT = ["NA", "NA", "NA"]
        self.MHCII_score_all_epitopes_WT = ["NA", "NA", "NA"]
        self.MHCII_score_best_per_alelle_WT = ["NA", "NA", "NA"]
        self.MHCII_number_strong_binders_WT = "NA"
        self.MHCII_number_weak_binders_WT = "NA"
        self.best_mhcII_pan_score_WT = "NA"
        self.best_mhcII_pan_epitope_WT = "-"
        self.best_mhcII_pan_allele_WT = "NA"
        self.best_mhcII_affinity_WT = "NA"
        self.best_mhcII_affinity_epitope_WT = "-"
        self.best_mhcII_affinity_allele_WT = "NA"

    def MHCII_MB_score_best_per_allele(self, tuple_best_per_allele):
        '''returns list of multiple binding scores for mhcII considering best epitope per allele, applying different types of means (harmonic ==> PHRB-II, Marty et al).
        2 copies of DRA - DRB1 --> consider this gene 2x when averaging mhcii binding scores
        '''
        number_alleles = len(tuple_best_per_allele)
        multbind = multiple_binders.MultipleBinding()
        tuple_best_per_allele_new = list(tuple_best_per_allele)
        logger.debug(tuple_best_per_allele)
        logger.debug(len(tuple_best_per_allele))
        for best_epi in tuple_best_per_allele:
            if best_epi[-1].startswith("DRB1"):
                tuple_best_per_allele_new.append(best_epi)
        logger.debug(tuple_best_per_allele_new)
        logger.debug(len(tuple_best_per_allele_new))
        if len(tuple_best_per_allele_new) == 12:
            # 12 genes gene copies should be included into PHBR_II
            best_scores_allele = multbind.scores_to_list(tuple_best_per_allele_new)
            return multbind.get_means(best_scores_allele)
        else:
            return ["NA", "NA", "NA"]

    def main(self, sequence, sequence_reference, alleles, set_available_mhc):
        '''predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        '''
        ### PREDICTION FOR MUTATED SEQUENCE
        logger.info("MUT seq MHC II: {}".format(sequence))
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        logger.debug(tmp_prediction)
        np = netmhcIIpan_prediction.NetMhcIIPanPredictor(runner=self.runner, configuration=self.configuration)
        mb = multiple_binders.MultipleBinding()
        tmp_fasta = intermediate_files.create_temp_fasta([sequence], prefix="tmp_singleseq_")
        alleles_formated = np.generate_mhcII_alelles_combination_list(alleles, set_available_mhc)
        logger.debug(alleles_formated)
        np.mhcII_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        position_xmer_sequence = np.mut_position_xmer_seq(xmer_wt=sequence_reference, xmer_mut=sequence)
        try:
            preds = np.filter_binding_predictions(position_xmer_sequence, tmp_prediction)
            # multiple binding
            list_tups = mb.generate_epi_tuple(preds, mhc=MHC_II)
            self.MHCII_epitope_scores = "/".join([tup[0] for tup in list_tups])
            self.MHCII_epitope_seqs = "/".join([tup[2] for tup in list_tups])
            self.MHCII_epitope_alleles = "/".join([tup[3] for tup in list_tups])
            top10 = mb.extract_top10_epis(list_tups)
            best_per_alelle = mb.extract_best_epi_per_alelle(list_tups, alleles_formated)
            all = mb.scores_to_list(list_tups)
            top10 = mb.scores_to_list(top10)
            self.MHCII_score_top10 = mb.get_means(top10)
            self.MHCII_score_all_epitopes = mb.get_means(all)
            self.MHCII_score_best_per_alelle = self.MHCII_MB_score_best_per_allele(best_per_alelle)
            self.MHCII_number_strong_binders = mb.determine_number_of_binders(all, 2)
            self.MHCII_number_weak_binders = mb.determine_number_of_binders(all, 10)
            logger.debug(self.MHCII_number_weak_binders)
            # best prediction
            best_epi = np.minimal_binding_score(preds)
            self.best_mhcII_pan_score = np.add_best_epitope_info(best_epi, "%Rank")
            self.best_mhcII_pan_epitope = np.add_best_epitope_info(best_epi, "Peptide")
            self.best_mhcII_pan_allele = np.add_best_epitope_info(best_epi, "Allele")
            self.best_mhcII_pan_position = np.add_best_epitope_info(best_epi, "Seq")
            best_epi_affinity = np.minimal_binding_score(preds, rank=False)
            self.best_mhcII_pan_affinity = np.add_best_epitope_info(best_epi_affinity, "Affinity(nM)")
            self.best_mhcII_pan_affinity_epitope = np.add_best_epitope_info(best_epi_affinity, "Peptide")
            self.best_mhcII_pan_affinity_allele = np.add_best_epitope_info(best_epi_affinity, "Allele")
            self.best_mhcII_pan_affinity_position = np.add_best_epitope_info(best_epi_affinity, "Seq")
        except IndexError:
            # if input sequence shorter than 15 aa
            pass

        ### PREDICTION FOR WT SEQUENCE
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        logger.debug(tmp_prediction)
        np = netmhcIIpan_prediction.NetMhcIIPanPredictor(runner=self.runner, configuration=self.configuration)
        mb = multiple_binders.MultipleBinding()
        tmp_fasta = intermediate_files.create_temp_fasta([sequence_reference], prefix="tmp_singleseq_")
        np.mhcII_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        try:
            preds = np.filter_binding_predictions(position_xmer_sequence, tmp_prediction)
            # multiple binding
            list_tups = mb.generate_epi_tuple(preds, mhc=MHC_II)
            self.MHCII_epitope_scores_WT = "/".join([tup[0] for tup in list_tups])
            self.epitope_affinities__mhcII_pan_WT = "/".join([tup[1] for tup in list_tups])
            self.MHCII_epitope_seqs_WT = "/".join([tup[2] for tup in list_tups])
            self.MHCII_epitope_alleles_WT = "/".join([tup[3] for tup in list_tups])
            top10 = mb.extract_top10_epis(list_tups)
            best_per_alelle = mb.extract_best_epi_per_alelle(list_tups, alleles_formated)
            all = mb.scores_to_list(list_tups)
            all_affinities = mb.affinities_to_list(list_tups)
            top10 = mb.scores_to_list(top10)
            self.MHCII_score_top10_WT = mb.get_means(top10)
            self.MHCII_score_all_epitopes_WT = mb.get_means(all)
            self.MHCII_score_best_per_alelle_WT = self.MHCII_MB_score_best_per_allele(best_per_alelle)
            self.MHCII_number_strong_binders_WT = mb.determine_number_of_binders(all, 1)
            self.MHCII_number_weak_binders_WT = mb.determine_number_of_binders(all, 2)
            # best prediction
            best_epi = np.filter_for_wt_epitope_position(preds, self.best_mhcII_pan_epitope,
                                                         position_epi_xmer=self.best_mhcII_pan_position)
            self.best_mhcII_pan_score_WT = np.add_best_epitope_info(best_epi, "%Rank")
            self.best_mhcII_pan_epitope_WT = np.add_best_epitope_info(best_epi, "Peptide")
            self.best_mhcII_pan_allele_WT = np.add_best_epitope_info(best_epi, "Allele")
            best_epi_affinity = np.filter_for_wt_epitope_position(preds, self.best_mhcII_pan_affinity_epitope,
                                                                  position_epi_xmer=self.best_mhcII_pan_affinity_position)
            self.best_mhcII_affinity_WT = np.add_best_epitope_info(best_epi_affinity, "Affinity(nM)")
            self.best_mhcII_affinity_epitope_WT = np.add_best_epitope_info(best_epi_affinity, "Peptide")
            self.best_mhcII_affinity_allele_WT = np.add_best_epitope_info(best_epi_affinity, "Allele")
        except IndexError:
            # if input sequence shorter than 15 aa
            pass


# if __name__ == '__main__':
#
#     from input import predict_all_epitopes, epitope
#     from input.helpers import data_import
#
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
#     # hla_file ="/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/20190730_alleles.csv"
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/test_ott_head_pt10.txt"
#     # hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/20190819_alleles_extended.csv"
#     # hla_file = "/projects/SUMMIT/WP1.2//Literature_Cohorts/data_analysis/cohorts/hugo/icam_hugo/20190819_alleles_extended.csv"
#     # file = "/projects/SUMMIT/WP1.2//Literature_Cohorts/data_analysis/cohorts/hugo/icam_hugo/Pt10/scratch/Pt10_mut_set.txt.transcript.squish.somatic.freq"
#     # file = "/projects/SUMMIT/WP1.2/dataset_annotation/Birmingham/in_files/PtBI000048T_1PEB.transcript"
#     # hla_file = "/projects/SUMMIT/WP1.2/dataset_annotation/Birmingham/20190822_alleles.csv"
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/test_ott_head_pt10.txt"
#     # hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/20190819_alleles_extended.csv"
#     # file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/rizvi/output_tables/freq_expression_input1.1/PtTU0428.rizvi.freq.expr.imput.txt"
#     file = "/projects/SUMMIT/WP1.2/input/development/netmhcIIpan/PtCU9061.test.txt"
#     hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/rizvi/icam_rizvi/20190819_alleles_extended.csv"
#     dat = data_import.import_dat_icam(file, False)
#     if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
#         dat = data_import.change_col_names(dat)
#     if "patient.id" not in dat[0]:
#         try:
#             patient = file.split("/")[-3]
#             if "Pt" not in patient:
#                 patient = file.split("/")[-1].split(".")[0]
#         except IndexError:
#             patient = file.split("/")[-1].split(".")[0]
#         dat[0].append("patient.id")
#         for ii, i in enumerate(dat[1]):
#             dat[1][ii].append(str(patient))
#     # available MHC alleles
#     set_available_mhc = predict_all_epitopes.Bunchepitopes().load_available_hla_alleles(mhc=MHC_II)
#     # hla allele of patients
#     patient_hlaI = predict_all_epitopes.Bunchepitopes().load_patient_hla_I_allels(hla_file)
#     patient_hlaII = predict_all_epitopes.Bunchepitopes().load_patient_hla_II_allels(hla_file)
#
#     Allepit = {}
#     for ii, i in enumerate(dat[1]):
#         if ii < 10:
#             dict_epi = epitope.Epitope()
#             dict_epi.init_properties(dat[0], dat[1][ii])
#             x = BestAndMultipleBinderMhcII()
#             x.main(dict_epi.properties, patient_hlaII, set_available_mhc)
#             # print x.MHC_epitope_scores_WT
#             # print x.MHC_epitope_seqs_WT
#             # print x.MHC_epitope_seqs
#             attrs = vars(x)
#             # print "score"
#             # print x.best_mhcII_pan_epitope
#             # print x.best_mhcII_pan_epitope_WT
#             # print "affinity"
#             # print x.best_mhcII_pan_affinity_epitope
#             # print x.best_mhcII_affinity_epitope_WT
#
#             print("PHBR")
#             print(x.MHCII_score_best_per_alelle)
#             # print attrs
