#!/usr/bin/env python

import os
import sys
import tempfile

from helpers import data_import
import FeatureLiterature
from self_similarity import self_similarity
from neoantigen_fitness import neoantigen_fitness
from new_features import differential_expression
from new_features import amino_acid_frequency_scores as freq_score
from new_features import conservation_scores
from aa_index import aa_index
from netmhcpan4 import combine_netmhcpan_pred_multiple_binders as mhcprediction



class Epitope:
    def __init__(self):
        self.properties = {}

    def init_properties(self, col_nam, prop_list):
        """Initiates epitope property storage in a dictionary
        """
        for nam,char in zip(col_nam, prop_list):
            self.properties[nam] = char

    def add_features(self, new_feature, new_feature_nam):
        """Adds new features to already present epitope properties, stored in form of a dictioninary
        """
        self.properties[new_feature_nam] = new_feature

    def write_to_file(self):
        print ";".join([self.properties[key] for key in self.properties])

    def main(self, col_nam, prop_list, db, ref_dat, aa_freq_dict, nmer_freq_dict, aaindex1_dict, aaindex2_dict, prov_matrix, set_available_mhc, patient_hlaI):
        """ Calculate new epitope features and add to dictonary that stores all properties
        """
        self.init_properties(col_nam, prop_list)
        # selfsimilarity
        self.add_features(self_similarity.selfsimilarity(self.properties,"mhcI"), "Selfsimilarity_mhcI")
        self.add_features(self_similarity.selfsimilarity(self.properties,"mhcII"), "Selfsimilarity_mhcII")
        self.add_features(self_similarity.improved_binder(self.properties,"mhcI"), "ImprovedBinding_mhcI")
        self.add_features(self_similarity.improved_binder(self.properties,"mhcII"), "ImprovedBinding_mhcII")
        self.add_features(self_similarity.position_of_mutation_epitope(self.properties,"mhcI"), "pos_MUT_MHCI")
        self.add_features(self_similarity.position_of_mutation_epitope(self.properties,"mhcII"), "pos_MUT_MHCII")
        self.add_features(self_similarity.position_in_anchor_position(self.properties), "Mutation_in_anchor")
        # neoantigen fitness
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix ="tmpseq", suffix = ".fasta", delete = False)
        tmp_fasta = tmp_fasta_file.name
        self.add_features(neoantigen_fitness.wrap_pathogensimilarity(self.properties, "mhcI", tmp_fasta), "Pathogensimiliarity_mhcI")
        self.add_features(neoantigen_fitness.wrap_pathogensimilarity(self.properties, "mhcII", tmp_fasta), "Pathogensimiliarity_mhcII")
        self.add_features(neoantigen_fitness.amplitude_mhc(self.properties, "mhcI"), "Amplitude_mhcI")
        self.add_features(neoantigen_fitness.amplitude_mhc(self.properties, "mhcII"), "Amplitude_mhcII")
        self.add_features(neoantigen_fitness.recognition_potential(self.properties, "mhcI"), "Recognition_Potential_mhcI")
        self.add_features(neoantigen_fitness.recognition_potential(self.properties, "mhcII"), "Recognition_Potential_mhcII")
        # IEDB immunogenicity
        self.add_features(FeatureLiterature.calc_IEDB_immunogenicity(self.properties, "mhcI"), "IEDB_Immunogenicity_mhcI")
        self.add_features(FeatureLiterature.calc_IEDB_immunogenicity(self.properties, "mhcII"), "IEDB_Immunogenicity_mhcII")
        # differential agretopicity index
        self.add_features(FeatureLiterature.dai(self.properties, "mhcI"), "DAI_mhcI")
        self.add_features(FeatureLiterature.dai(self.properties, "mhcII"), "DAI_mhcII")
        # priority score
        self.add_features(FeatureLiterature.rna_expression_mutation(self.properties), "Expression_Mutated_Transcript")
        self.add_features(FeatureLiterature.number_of_mismatches(self.properties, "mhcI"), "Number_of_mismatches_mhcI")
        self.add_features(FeatureLiterature.number_of_mismatches(self.properties, "mhcII"), "Number_of_mismatches_mhcII")
        if "mutation_found_in_proteome" not in self.properties:
            self.add_features(FeatureLiterature.match_in_proteome(self.properties, db), "mutation_found_in_proteome")
        self.add_features(FeatureLiterature.calc_priority_score(self.properties), "Priority_score")
        self.add_features(FeatureLiterature.wt_mut_aa(self.properties, "mut"), "MUT_AA")
        self.add_features(FeatureLiterature.wt_mut_aa(self.properties, "wt"), "WT_AA")
        # differential expression
        self.add_features(differential_expression.add_rna_reference(self.properties, ref_dat, 0), "mean_ref_expression")
        self.add_features(differential_expression.add_rna_reference(self.properties, ref_dat, 1), "sd_ref_expression")
        self.add_features(differential_expression.add_rna_reference(self.properties, ref_dat, 2), "sum_ref_expression")
        self.add_features(differential_expression.fold_change(self.properties), "log2_fc_tumour_ref")
        self.add_features(differential_expression.percentile_calc(self.properties), "percentile_tumour_ref")
        self.add_features(differential_expression.pepper_calc(self.properties), "DE_pepper")
        # amino acid frequency
        self.add_features(freq_score.freq_aa(self.properties, aa_freq_dict), "Frequency_mutated_AA")
        self.add_features(freq_score.freq_prod_4mer(self.properties, aa_freq_dict), "Product_Frequency_4mer")
        self.add_features(freq_score.freq_4mer(self.properties, nmer_freq_dict), "Frequency_of_4mer")
        for k in aaindex1_dict:
            z = FeatureLiterature.add_aa_index1(self.properties, "wt", k, aaindex1_dict[k])
            self.add_features(z[1],z[0])
            z = FeatureLiterature.add_aa_index1(self.properties, "mut", k, aaindex1_dict[k])
            self.add_features(z[1],z[0])
        for k in aaindex2_dict:
            try:
                z = FeatureLiterature.add_aa_index2(self.properties, k, aaindex2_dict[k])
                self.add_features(z[1],z[0])
            except:
                print aaindex2_dict[k], wt, mut
        # PROVEAN score
        self.add_features(conservation_scores.add_ucsc_id_to_dict(self.properties), "UCSC_ID_position")
        self.add_features(conservation_scores.add_provean_score_from_matrix(self.properties, prov_matrix), "PROVEAN_score")
        pred = mhcprediction.Bestandmultiplebinder()
        pred.main(self.properties, patient_hlaI, set_available_mhc)
        # netmhcpan4 MUT
        self.add_features(pred.best4_mhc_score, "best%Rank_netmhcpan4")
        self.add_features(pred.best4_mhc_epitope, "best_epitope_netmhcpan4")
        self.add_features(pred.best4_mhc_allele, "bestHLA_allele_netmhcpan4")
        self.add_features(pred.directed_to_TCR, "directed_to_TCR")
        # multiplexed representation MUT
        for sc, mn in zip(pred.MHC_score_all_epitopes, pred.mean_type):
            self.add_features(sc, "MB_score_all_epitopes_" + mn)
        for sc, mn in zip(pred.MHC_score_top10, pred.mean_type):
            self.add_features(sc, "MB_score_top10_" + mn)
        for sc, mn in zip(pred.MHC_score_best_per_alelle, pred.mean_type):
            self.add_features(sc, "MB_score_best_per_alelle_" + mn)
        self.add_features(pred.MHC_epitope_scores, "MB_epitope_scores")
        self.add_features(pred.MHC_epitope_seqs, "MB_epitope_sequences")
        self.add_features(pred.MHC_epitope_alleles, "MB_alleles")
        self.add_features(pred.MHC_number_strong_binders, "MB_number_pep_MHCscore<1")
        self.add_features(pred.MHC_number_weak_binders, "MB_number_pep_MHCscore<2")
        # multiplexed representation WT
        self.add_features(pred.MHC_epitope_scores_WT, "MB_epitope_WT_scores")
        self.add_features(pred.MHC_epitope_seqs_WT, "MB_epitope_WT_sequences")
        self.add_features(pred.MHC_epitope_alleles_WT, "MB_alleles_WT")
        for sc, mn in zip(pred.MHC_score_top10_WT, pred.mean_type):
            self.add_features(sc, "MB_score_WT_top10_" + mn)
        for sc, mn in zip(pred.MHC_score_all_epitopes_WT, pred.mean_type):
            self.add_features(sc, "MB_score_WT_all_epitopes_" + mn)
        for sc, mn in zip(pred.MHC_score_best_per_alelle_WT, pred.mean_type):
            self.add_features(sc, "MB_score_WT_best_per_alelle_" + mn)
        self.add_features(pred.MHC_number_strong_binders_WT, "MB_number_pep_WT_MHCscore<1")
        self.add_features(pred.MHC_number_weak_binders_WT, "MB_number_pep_WT_MHCscore<2")
        self.add_features(FeatureLiterature.dai(self.properties, "mhcI", True), "DAI_mhcI_MB")
        self.add_features(FeatureLiterature.calc_priority_score(self.properties, True), "Priority_score_MB")
        self.add_features(FeatureLiterature.diff_number_binders(self.properties, "1"), "Diff_numb_epis_<1")
        self.add_features(FeatureLiterature.diff_number_binders(self.properties, "2"), "Diff_numb_epis_<2")
        self.add_features(FeatureLiterature.ratio_number_binders(self.properties, "1"), "Ratio_numb_epis_<1")
        self.add_features(FeatureLiterature.ratio_number_binders(self.properties, "2"), "Ratio_numb_epis_<2")
        self.add_features(neoantigen_fitness.amplitude_mhc(self.properties, "mhcI", multiple_binding=True), "Amplitude_mhcI_MB")
        # priority score using multiplexed representation score


        return self.properties


if __name__ == '__main__':
    import sys
    from datetime import datetime
    import predict_all_epitopes as predfeatall

    startTime = datetime.now()
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/nonprogramm_files/test_SD.csv"
    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/nonprogramm_files/test_fulldat.txt"
    indel = False
    fasta_proteome = "/projects/data/human/2018_uniprot_with_isoforms/uniprot_human_with_isoforms.fasta"
    ref_file = "/projects/CM27_IND_patients/GTEX_normal_tissue_data/Skin .csv"

    #predfeatallBunchepitopes

    z = Epitope().main(dat[0], dat[1][ii], self.proteome_dictionary, self.rna_reference, self.aa_frequency, self.fourmer_frequency, self.aa_index1_dict, self.aa_index2_dict, self.provean_matrix, self.hla_available_alleles, self.patient_hla_I_allels)

    Bunchepitopes().main(file, indel, fasta_proteome, ref_file)
    endTime = datetime.now()
    print >> sys.stderr, "start: "+ str(startTime) + "\nend: "+ str(endTime) + "\nneeded: " + str(endTime - startTime)
    #print dat
    #x = Epitope()
    #x = Epitope(dat[1][1], dat[0])
    #print vars(x)
    #print dat[1][1][1]
    #print dat[0][1]


    #for ii,i in enumerate(dat[1]):
    #    Epitope().main(dat[0],dat[1][ii])
        #print x.tricks


    #x.main(dat[0], dat[1][1])
    #print x.tricks
    #print x.tricks["transcript_position"]
    #print dir(x)
