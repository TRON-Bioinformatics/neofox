#!/usr/bin/env python

from logzero import logger

import input.predictors.netmhcpan4.multiple_binders as multiple_binders
import input.predictors.netmhcpan4.netmhcpan_prediction as netmhcpan_prediction
from input.helpers import intermediate_files


class BestAndMultipleBinder:
    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
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
        self.best4_mhc_position = "NA"
        self.directed_to_TCR = "NA"
        self.best4_affinity = "NA"
        self.best4_affinity_epitope = "-"
        self.best4_affinity_allele = "NA"
        self.best4_affinity_position = "NA"
        self.best4_affinity_directed_to_TCR = "NA"
        self.epitope_affinities = ""
        self.generator_rate = ""
        self.mhcI_score_9mer = "NA"
        self.mhcI_score_allele_9mer = "NA"
        self.mhcI_score_position_9mer = "NA"
        self.mhcI_score_epitope_9mer = "-"
        self.mhcI_affinity_9mer = "NA"
        self.mhcI_affinity_allele_9mer = "NA"
        self.mhcI_affinity_position_9mer = "NA"
        self.mhcI_affinity_epitope_9mer = "-"
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
        self.mhcI_score_9mer_WT = "NA"
        self.mhcI_score_allele_9mer_WT = "NA"
        self.mhcI_score_epitope_9mer_WT = "NA"
        self.mhcI_affinity_9mer_WT = "NA"
        self.mhcI_affinity_allele_9mer_WT = "NA"
        self.mhcI_affinity_epitope_9mer_WT = "NA"

    def mhc_mb_score_best_per_allele(self, tuple_best_per_allele):
        '''returns list of multiple binding scores for mhcII considering best epitope per allele, applying different types of means (harmonic ==> PHRB-II, Marty et al).
        2 copies of DRA - DRB1 --> consider this gene 2x when averaging mhcii binding scores
        '''
        multbind = multiple_binders.MultipleBinding()
        tuple_best_per_allele_new = list(tuple_best_per_allele)
        if len(tuple_best_per_allele_new) == 6:
            mean_best_per_per_allele = multbind.get_means(tuple_best_per_allele_new)
        else:
            mean_best_per_per_allele = ["NA", "NA", "NA"]
        return mean_best_per_per_allele

    def main(self, xmer_wt, xmer_mut, alleles, set_available_mhc):
        """
        predicts MHC epitopes; returns on one hand best binder and on the other hand multiple binder analysis is performed
        """
        ### PREDICTION FOR MUTATED SEQUENCE
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        np = netmhcpan_prediction.NetMhcPanPredictor(runner=self.runner, configuration=self.configuration)
        mb = multiple_binders.MultipleBinding()
        tmp_fasta = intermediate_files.create_temp_fasta(sequences=[xmer_mut], prefix="tmp_singleseq_")
        # print alleles
        np.mhc_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)

        position_xmer = np.mut_position_xmer_seq(xmer_mut=xmer_mut, xmer_wt=xmer_wt)
        preds = np.filter_binding_predictions(position_xmer=position_xmer, tmppred=tmp_prediction)

        # multiple binding
        list_tups = mb.generate_epi_tuple(preds)
        self.MHC_epitope_scores = "/".join([tup[0] for tup in list_tups])
        self.epitope_affinities = "/".join([tup[1] for tup in list_tups])
        self.MHC_epitope_seqs = "/".join([tup[2] for tup in list_tups])
        self.MHC_epitope_alleles = "/".join([tup[3] for tup in list_tups])
        top10 = mb.extract_top10_epis(list_tups)
        best_per_alelle = mb.extract_best_epi_per_alelle(list_tups, alleles)
        logger.debug("sdfsd")
        logger.debug(alleles)
        all = mb.scores_to_list(list_tups)
        all_affinities = mb.affinities_to_list(list_tups)
        top10 = mb.scores_to_list(top10)
        self.MHC_score_top10 = mb.get_means(top10)
        best_per_alelle = mb.scores_to_list(best_per_alelle)
        self.MHC_score_all_epitopes = mb.get_means(all)
        self.MHC_score_best_per_alelle = self.mhc_mb_score_best_per_allele(best_per_alelle)
        self.MHC_number_strong_binders = mb.determine_number_of_binders(all, 1)
        self.MHC_number_weak_binders = mb.determine_number_of_binders(all, 2)
        # best prediction
        best_epi = np.minimal_binding_score(preds)
        self.best4_mhc_score = np.add_best_epitope_info(best_epi, "%Rank")
        self.best4_mhc_epitope = np.add_best_epitope_info(best_epi, "Peptide")
        self.best4_mhc_allele = np.add_best_epitope_info(best_epi, "HLA")
        self.best4_mhc_position = np.add_best_epitope_info(best_epi, "Pos")
        self.directed_to_TCR = np.mutation_in_loop(position_xmer_list=position_xmer, epitope_tuple=best_epi)
        best_epi_affinity = np.minimal_binding_score(preds, rank=False)
        self.best4_affinity = np.add_best_epitope_info(best_epi_affinity, "Aff(nM)")
        self.best4_affinity_epitope = np.add_best_epitope_info(best_epi_affinity, "Peptide")
        self.best4_affinity_allele = np.add_best_epitope_info(best_epi_affinity, "HLA")
        self.best4_affinity_position = np.add_best_epitope_info(best_epi_affinity, "Pos")
        self.best4_affinity_directed_to_TCR = np.mutation_in_loop(
            position_xmer_list=position_xmer, epitope_tuple=best_epi_affinity)
        # multiple binding based on affinity
        self.generator_rate = mb.determine_number_of_binders(list_scores=all_affinities, threshold=50)
        # best predicted epitope of length 9
        preds_9mer = np.filter_for_9mers(preds)
        best_9mer = np.minimal_binding_score(preds_9mer)
        best_9mer_affinity = np.minimal_binding_score(preds_9mer, rank=False)
        self.mhcI_score_9mer = np.add_best_epitope_info(best_9mer, "%Rank")
        self.mhcI_score_allele_9mer = np.add_best_epitope_info(best_9mer, "HLA")
        self.mhcI_score_position_9mer = np.add_best_epitope_info(best_9mer, "Pos")
        self.mhcI_score_epitope_9mer = np.add_best_epitope_info(best_9mer, "Peptide")
        self.mhcI_affinity_9mer = np.add_best_epitope_info(best_9mer_affinity, "Aff(nM)")
        self.mhcI_affinity_allele_9mer = np.add_best_epitope_info(best_9mer_affinity, "HLA")
        self.mhcI_affinity_position_9mer = np.add_best_epitope_info(best_9mer_affinity, "Pos")
        self.mhcI_affinity_epitope_9mer = np.add_best_epitope_info(best_9mer_affinity, "Peptide")

        ### PREDICTION FOR WT SEQUENCE
        tmp_prediction = intermediate_files.create_temp_file(prefix="netmhcpanpred_", suffix=".csv")
        logger.debug(tmp_prediction)
        np = netmhcpan_prediction.NetMhcPanPredictor(runner=self.runner, configuration=self.configuration)
        mb = multiple_binders.MultipleBinding()
        tmp_fasta = intermediate_files.create_temp_fasta(sequences=[xmer_wt],
                                                         prefix="tmp_singleseq_")
        np.mhc_prediction(alleles, set_available_mhc, tmp_fasta, tmp_prediction)
        preds = np.filter_binding_predictions(position_xmer=position_xmer, tmppred=tmp_prediction)
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
        self.MHC_score_top10_WT = mb.get_means(top10)
        best_per_alelle = mb.scores_to_list(best_per_alelle)
        self.MHC_score_all_epitopes_WT = mb.get_means(all)
        self.MHC_score_best_per_alelle_WT = self.mhc_mb_score_best_per_allele(best_per_alelle)
        self.MHC_number_strong_binders_WT = mb.determine_number_of_binders(all, 1)
        self.MHC_number_weak_binders_WT = mb.determine_number_of_binders(all, 2)
        # best prediction
        best_epi = np.filter_for_WT_epitope_position(preds, self.best4_mhc_epitope,
                                            position_epi_xmer=self.best4_mhc_position)
        self.best4_mhc_score_WT = np.add_best_epitope_info(best_epi, "%Rank")
        self.best4_mhc_epitope_WT = np.add_best_epitope_info(best_epi, "Peptide")
        self.best4_mhc_allele_WT = np.add_best_epitope_info(best_epi, "HLA")

        best_epi_affinity = np.filter_for_WT_epitope_position(preds, self.best4_affinity_epitope,
                                                     position_epi_xmer=self.best4_affinity_position)
        self.best4_affinity_WT = np.add_best_epitope_info(best_epi_affinity, "Aff(nM)")
        self.best4_affinity_epitope_WT = np.add_best_epitope_info(best_epi_affinity, "Peptide")
        self.best4_affinity_allele_WT = np.add_best_epitope_info(best_epi_affinity, "HLA")
        self.generator_rate_WT = mb.determine_number_of_binders(list_scores=all_affinities, threshold=50)
        logger.info("WT: {}; MUT: {}".format(self.generator_rate_WT, self.generator_rate))
        # best predicted epitope of length 9
        preds_9mer = np.filter_for_9mers(preds)
        best_9mer = np.filter_for_WT_epitope_position(preds_9mer, self.mhcI_score_epitope_9mer,
                                             position_epi_xmer=self.mhcI_score_position_9mer)
        best_9mer_affinity = np.filter_for_WT_epitope_position(preds_9mer, mut_seq=self.mhcI_affinity_epitope_9mer,
                                                      position_epi_xmer=self.mhcI_affinity_position_9mer)
        self.mhcI_score_9mer_WT = np.add_best_epitope_info(best_9mer, "%Rank")
        self.mhcI_score_allele_9mer_WT = np.add_best_epitope_info(best_9mer, "HLA")
        self.mhcI_score_epitope_9mer_WT = np.add_best_epitope_info(best_9mer, "Peptide")
        self.mhcI_affinity_9mer_WT = np.add_best_epitope_info(best_9mer_affinity, "Aff(nM)")
        self.mhcI_affinity_allele_9mer_WT = np.add_best_epitope_info(best_9mer_affinity, "HLA")
        self.mhcI_affinity_epitope_9mer_WT = np.add_best_epitope_info(best_9mer_affinity, "Peptide")


# if __name__ == '__main__':
#
#     from input import predict_all_epitopes, epitope
#     from input.helpers import data_import
#
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_SD.csv"
#     # file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/riaz/output_tables_pre/test.txt"
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/20170713_IS_IM_data.complete.update_Dv10.csv.annotation.csv_v2.csv"
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT/nonprogramm_files/test_fulldat.txt"
#     # hla_file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/indels/RB_0004_labHLA_V2.csv"
#     # hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/riaz/output_tables_pre/alleles.csv"
#     # file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/test_ott_head_pt10.txt"
#     # hla_file ="/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
#     # hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/20190819_alleles_extended.csv"
#     hla_file = "/projects/SUMMIT/WP1.2/dataset_annotation/Birmingham/20190822_alleles.csv"
#     file = "/projects/SUMMIT/WP1.2/dataset_annotation/Birmingham/in_files/PtBI000048T_1PEB.transcript"
#     # hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/hugo/icam_hugo/20190816_alleles_extended.csv"
#     # file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/hugo/icam_hugo/Pt12/scratch/Pt12_mut_set.txt.transcript.squish.somatic.freq"
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
#     set_available_mhc = predict_all_epitopes.Bunchepitopes().load_available_hla_alleles()
#     # hla allele of patients
#     patient_hlaI = predict_all_epitopes.Bunchepitopes().load_patient_hla_I_allels(hla_file)
#     patient_hlaII = predict_all_epitopes.Bunchepitopes().load_patient_hla_II_allels(hla_file)
#
#     Allepit = {}
#     for ii, i in enumerate(dat[1]):
#         if ii < 10:
#             dict_epi = epitope.Epitope()
#             dict_epi.init_properties(dat[0], dat[1][ii])
#             x = Bestandmultiplebinder()
#             x.main(dict_epi.properties, patient_hlaI, set_available_mhc)
#             print(dict_epi.properties["patient.id"])
#             attrs = vars(x)
#             print(x.mhcI_affinity_epitope_9mer)
#             print(x.MHC_score_best_per_alelle)
#
#             '''
#             for sc, mn in zip(x.MHC_score_all_epitopes, x.mean_type):
#                 dict_epi.add_features(sc, "MB_score_all_epitopes_" + mn)
#             for sc, mn in zip(x.MHC_score_top10, x.mean_type):
#                 dict_epi.add_features(sc, "MB_score_top10_" + mn)
#             for sc, mn in zip(x.MHC_score_best_per_alelle, x.mean_type):
#                 dict_epi.add_features(sc, "MB_score_best_per_alelle_" + mn)
#             dict_epi.add_features(x.MHC_epitope_scores, "MB_epitope_scores")
#             dict_epi.add_features(x.MHC_epitope_seqs, "MB_epitope_sequences")
#             dict_epi.add_features(x.MHC_epitope_alleles, "MB_alleles")
#             dict_epi.add_features(x.MHC_number_strong_binders, "MB_number_of_strong_binders")
#             dict_epi.add_features(x.MHC_number_weak_binders, "MB_number_of_weak_binders")
#             dict_epi.add_features(x.best4_mhc_score, "best4_mhc_score")
#             dict_epi.add_features(x.best4_mhc_epitope, "best4_mhc_epitope")
#             dict_epi.add_features(x.best4_mhc_allele, "best4_mhc_allele")
#             dict_epi.add_features(x.directed_to_TCR, "directed_to_TCR")
#             z = dict_epi.properties
#             for key in z:
#                 if key not in Allepit:
#                     Allepit[key] = [z[key]]
#                 else:
#                     Allepit[key].append(z[key])
#             '''
#     # predict_all_epitopes.Bunchepitopes().write_to_file(Allepit)
