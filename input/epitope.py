#!/usr/bin/env python

import tempfile

from logzero import logger

from input import FeatureLiterature
from input import MHC_I, MHC_II
from input.MixMHCpred.mixmhcpred import MixMHCpred
from input.MixMHCpred.mixmhc2pred import MixMhc2Pred
from input.Tcell_predictor.tcellpredictor_wrapper import Tcellprediction
from input.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator
from input.neoag.neoag_gbm_model import NeoagCalculator
from input.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from input.netmhcIIpan.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from input.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from input.new_features import amino_acid_frequency_scores as freq_score, differential_expression, conservation_scores
from input.self_similarity import self_similarity
from input.vaxrank import vaxrank
from input.helpers.runner import Runner


class Epitope:

    def __init__(self, runner, references, configuration, provean_annotator):
        """
        :type runner: input.helpers.runner.Runner
        :type references: input.references.ReferenceFolder
        :type configuration: input.references.DependenciesConfiguration
        :type provean_annotator: input.new_features.conservation_scores.ProveanAnnotator
        """
        self.references = references
        self.provean_annotator = provean_annotator
        self.properties = {}
        self.dissimilarity_calculator = DissimilarityCalculator(runner=runner, configuration=configuration)
        self.neoantigen_fitness_calculator = NeoantigenFitnessCalculator(runner=runner, configuration=configuration)
        self.neoag_calculator = NeoagCalculator(runner=runner, configuration=configuration)
        self.predII = BestAndMultipleBinderMhcII(runner=runner, configuration=configuration)
        self.predpresentation2 = MixMhc2Pred(runner=runner, configuration=configuration)
        self.pred = BestAndMultipleBinder(runner=runner, configuration=configuration)
        self.predpresentation = MixMHCpred(runner=runner, configuration=configuration)
        self.tcellpredict = Tcellprediction(references=self.references)

    def init_properties(self, col_nam, prop_list):
        """Initiates epitope property storage in a dictionary
        """
        for nam, char in zip(col_nam, prop_list):
            self.properties[nam] = char

    def add_features(self, new_feature, new_feature_nam):
        """Adds new features to already present epitope properties, stored in form of a dictioninary
        """
        self.properties[new_feature_nam] = new_feature

    def write_to_file(self):
        print(";".join([self.properties[key] for key in self.properties]))

    def main(self, col_nam, prop_list, db, ref_dat, aa_freq_dict, nmer_freq_dict, aaindex1_dict, aaindex2_dict,
             set_available_mhc, set_available_mhcII, patient_hlaI, patient_hlaII, tumour_content,
             list_HLAII_MixMHC2pred, rna_avail):
        """ Calculate new epitope features and add to dictonary that stores all properties
        """
        self.init_properties(col_nam, prop_list)
        logger.info(self.properties["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."])

        self.add_features(self_similarity.position_of_mutation_epitope(self.properties, MHC_I), "pos_MUT_MHCI")
        self.add_features(self_similarity.position_of_mutation_epitope(self.properties, MHC_II), "pos_MUT_MHCII")
        self.add_features(self_similarity.position_in_anchor_position(self.properties), "Mutation_in_anchor")

        # differential agretopicity index
        self.add_features(FeatureLiterature.dai(self.properties, MHC_I), "DAI_mhcI")
        self.add_features(FeatureLiterature.dai(self.properties, MHC_II), "DAI_mhcII")
        # expression
        self.add_features(FeatureLiterature.rna_expression_mutation(self.properties, rna_avail=rna_avail),
                          "Expression_Mutated_Transcript")
        self.add_features(FeatureLiterature.expression_mutation_tc(self.properties, tumour_content=tumour_content),
                          "Expression_Mutated_Transcript_tumor_content")

        # differential expression
        self.add_features(differential_expression.add_rna_reference(self.properties, ref_dat, 0), "mean_ref_expression")
        self.add_features(differential_expression.add_rna_reference(self.properties, ref_dat, 1), "sd_ref_expression")
        self.add_features(differential_expression.add_rna_reference(self.properties, ref_dat, 2), "sum_ref_expression")
        self.add_features(differential_expression.fold_change(self.properties), "log2_fc_tumour_ref")
        self.add_features(differential_expression.percentile_calc(self.properties), "percentile_tumour_ref")
        self.add_features(differential_expression.pepper_calc(self.properties), "DE_pepper")
        # amino acid frequency
        self.add_features(FeatureLiterature.wt_mut_aa(self.properties, "mut"), "MUT_AA")
        self.add_features(FeatureLiterature.wt_mut_aa(self.properties, "wt"), "WT_AA")
        self.add_features(freq_score.freq_aa(self.properties, aa_freq_dict), "Frequency_mutated_AA")
        self.add_features(freq_score.freq_prod_4mer(self.properties, aa_freq_dict), "Product_Frequency_4mer")
        self.add_features(freq_score.freq_4mer(self.properties, nmer_freq_dict), "Frequency_of_4mer")
        # amino acid index
        for k in aaindex1_dict:
            z = FeatureLiterature.add_aa_index1(self.properties, "wt", k, aaindex1_dict[k])
            self.add_features(z[1], z[0])
            z = FeatureLiterature.add_aa_index1(self.properties, "mut", k, aaindex1_dict[k])
            self.add_features(z[1], z[0])
        for k in aaindex2_dict:
            try:
                z = FeatureLiterature.add_aa_index2(self.properties, k, aaindex2_dict[k])
                self.add_features(z[1], z[0])
            except:
                print(aaindex2_dict[k], wt, mut)

        # PROVEAN score
        ucsc_id = self.provean_annotator.build_ucsc_id_plus_position(
            substitution=self.properties["substitution"], ucsc_id=self.properties["UCSC_transcript"])
        self.add_features(ucsc_id, "UCSC_ID_position")
        self.add_features(self.provean_annotator.get_provean_annotation(
            mutated_aminoacid=self.properties['MUT_AA'], ucsc_id_position=ucsc_id),
            "PROVEAN_score")
        self.pred.main(self.properties, patient_hlaI, set_available_mhc)

        # netmhcpan4 MUT rank score
        self.add_features(self.pred.best4_mhc_score, "best%Rank_netmhcpan4")
        self.add_features(self.pred.best4_mhc_epitope, "best_epitope_netmhcpan4")
        self.add_features(self.pred.best4_mhc_allele, "bestHLA_allele_netmhcpan4")
        self.add_features(self.pred.directed_to_TCR, "directed_to_TCR")

        # netmhcpan4 mut affinity
        self.add_features(self.pred.best4_affinity, "best_affinity_netmhcpan4")
        self.add_features(self.pred.best4_affinity_epitope, "best_affinity_epitope_netmhcpan4")
        self.add_features(self.pred.best4_affinity_allele, "bestHLA_allele_affinity_netmhcpan4")
        self.add_features(self.pred.best4_affinity_directed_to_TCR, "affinity_directed_to_TCR")

        # netMHCpan MUT best 9mer score
        self.add_features(self.pred.mhcI_score_9mer, "best%Rank_netmhcpan4_9mer")
        self.add_features(self.pred.mhcI_score_epitope_9mer, "best_epitope_netmhcpan4_9mer")
        self.add_features(self.pred.mhcI_score_allele_9mer, "bestHLA_allele_netmhcpan4_9mer")

        # netmhcpan4 mut best 9mer affinity
        self.add_features(self.pred.mhcI_affinity_9mer, "best_affinity_netmhcpan4_9mer")
        self.add_features(self.pred.mhcI_affinity_allele_9mer, "bestHLA_allele_affinity_netmhcpan4_9mer")
        self.add_features(self.pred.mhcI_affinity_epitope_9mer, "best_affinity_epitope_netmhcpan4_9mer")

        # multiplexed representation MUT
        for sc, mn in zip(self.pred.MHC_score_all_epitopes, self.pred.mean_type):
            self.add_features(sc, "MB_score_all_epitopes_" + mn)
        for sc, mn in zip(self.pred.MHC_score_top10, self.pred.mean_type):
            self.add_features(sc, "MB_score_top10_" + mn)
        for sc, mn in zip(self.pred.MHC_score_best_per_alelle, self.pred.mean_type):
            self.add_features(sc, "MB_score_best_per_alelle_" + mn)

        # rename MB_score_best_per_alelle_harmonic to PHBR (described in Marty et al)
        self.properties["PHBR-I"] = self.properties.pop("MB_score_best_per_alelle_harmonic")
        self.add_features(self.pred.MHC_epitope_scores, "MB_epitope_scores")
        self.add_features(self.pred.MHC_epitope_seqs, "MB_epitope_sequences")
        self.add_features(self.pred.MHC_epitope_alleles, "MB_alleles")
        self.add_features(self.pred.MHC_number_strong_binders, "MB_number_pep_MHCscore<1")
        self.add_features(self.pred.MHC_number_weak_binders, "MB_number_pep_MHCscore<2")

        # generator rate
        self.add_features(self.pred.epitope_affinities, "MB_affinities")
        self.add_features(self.pred.generator_rate, "Generator_rate")

        # multiplexed representation WT
        self.add_features(self.pred.MHC_epitope_scores_WT, "MB_epitope_WT_scores")
        self.add_features(self.pred.MHC_epitope_seqs_WT, "MB_epitope_WT_sequences")
        self.add_features(self.pred.MHC_epitope_alleles_WT, "MB_alleles_WT")
        for sc, mn in zip(self.pred.MHC_score_top10_WT, self.pred.mean_type):
            self.add_features(sc, "MB_score_WT_top10_" + mn)
        for sc, mn in zip(self.pred.MHC_score_all_epitopes_WT, self.pred.mean_type):
            self.add_features(sc, "MB_score_WT_all_epitopes_" + mn)
        for sc, mn in zip(self.pred.MHC_score_best_per_alelle_WT, self.pred.mean_type):
            self.add_features(sc, "MB_score_WT_best_per_alelle_" + mn)
        self.properties["PHBR-I_WT"] = self.properties.pop("MB_score_WT_best_per_alelle_harmonic")
        self.add_features(self.pred.MHC_number_strong_binders_WT, "MB_number_pep_WT_MHCscore<1")
        self.add_features(self.pred.MHC_number_weak_binders_WT, "MB_number_pep_WT_MHCscore<2")

        # generator rate
        self.add_features(self.pred.epitope_affinities_WT, "MB_affinities_WT")
        self.add_features(self.pred.generator_rate_WT, "Generator_rate_WT")
        self.add_features(FeatureLiterature.dai(self.properties, MHC_I, multiple_binding=True), "DAI_mhcI_MB")

        # netmhcpan4 wt affinity
        self.add_features(self.pred.best4_affinity_WT, "best_affinity_netmhcpan4_WT")
        self.add_features(self.pred.best4_affinity_epitope_WT, "best_affinity_epitope_netmhcpan4_WT")
        self.add_features(self.pred.best4_affinity_allele_WT, "bestHLA_allele_affinity_netmhcpan4_WT")

        # netmhcpan4 mut rank score
        self.add_features(self.pred.best4_mhc_score_WT, "best%Rank_netmhcpan4_WT")
        self.add_features(self.pred.best4_mhc_epitope_WT, "best_epitope_netmhcpan4_WT")
        self.add_features(self.pred.best4_mhc_allele_WT, "bestHLA_allele_netmhcpan4_WT")

        # netMHCpan MUT best 9mer score
        self.add_features(self.pred.mhcI_score_9mer_WT, "best%Rank_netmhcpan4_9mer_WT")
        self.add_features(self.pred.mhcI_score_epitope_9mer_WT, "best_epitope_netmhcpan4_9mer_WT")
        self.add_features(self.pred.mhcI_score_allele_9mer_WT, "bestHLA_allele_netmhcpan4_9mer_Wt")

        # netmhcpan4 mut best 9mer affinity
        self.add_features(self.pred.mhcI_affinity_9mer_WT, "best_affinity_netmhcpan4_9mer_WT")
        self.add_features(self.pred.mhcI_affinity_allele_9mer_WT, "bestHLA_allele_affinity_netmhcpan4_9mer_WT")
        self.add_features(self.pred.mhcI_affinity_epitope_9mer_WT, "best_affinity_epitope_netmhcpan4_9mer_WT")

        # multiplex representation
        self.add_features(FeatureLiterature.diff_number_binders(self.properties, mhc=MHC_I, threshold="1"),
                          "Diff_numb_epis_<1")
        self.add_features(FeatureLiterature.diff_number_binders(self.properties, mhc=MHC_I, threshold="2"),
                          "Diff_numb_epis_<2")
        self.add_features(FeatureLiterature.ratio_number_binders(self.properties, mhc=MHC_I, threshold="1"),
                          "Ratio_numb_epis_<1")
        self.add_features(FeatureLiterature.ratio_number_binders(self.properties, mhc=MHC_I, threshold="2"),
                          "Ratio_numb_epis_<2")
        self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            self.properties, MHC_I, multiple_binding=True), "Amplitude_mhcI_MB")

        # position of mutation
        self.add_features(self_similarity.position_of_mutation_epitope_affinity(self.properties),
                          "pos_MUT_MHCI_affinity_epi")

        # position of mutation
        self.add_features(self_similarity.position_of_mutation_epitope_affinity(self.properties, nine_mer=True),
                          "pos_MUT_MHCI_affinity_epi_9mer")
        self.add_features(self_similarity.position_in_anchor_position(self.properties, netMHCpan=True),
                          "Mutation_in_anchor_netmhcpan")
        self.add_features(self_similarity.position_in_anchor_position(self.properties, nine_mer=True),
                          "Mutation_in_anchor_netmhcpan_9mer")

        # selfsimilarity
        self.add_features(self_similarity.get_self_similarity(
            props=self.properties, mhc=MHC_I, references=self.references), "Selfsimilarity_mhcI")
        self.add_features(self_similarity.get_self_similarity(
            props=self.properties, mhc=MHC_II, references=self.references), "Selfsimilarity_mhcII")
        self.add_features(self_similarity.is_improved_binder(self.properties, MHC_I), "ImprovedBinding_mhcI")
        self.add_features(self_similarity.is_improved_binder(self.properties, MHC_II), "ImprovedBinding_mhcII")
        self.add_features(self_similarity.selfsimilarity_of_conserved_binder_only(self.properties),
                          "Selfsimilarity_mhcI_conserved_binder")

        # neoantigen fitness
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix="tmpseq", suffix=".fasta", delete=False)
        tmp_fasta = tmp_fasta_file.name
        self.add_features(
            self.neoantigen_fitness_calculator.wrap_pathogensimilarity(
                props=self.properties, mhc=MHC_I, fastafile=tmp_fasta, iedb=self.references.iedb),
            "Pathogensimiliarity_mhcI")
        self.add_features(
            self.neoantigen_fitness_calculator.wrap_pathogensimilarity(
                props=self.properties, mhc=MHC_II, fastafile=tmp_fasta, iedb=self.references.iedb),
            "Pathogensimiliarity_mhcII")
        self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            self.properties, MHC_I), "Amplitude_mhcI")
        self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            self.properties, MHC_II), "Amplitude_mhcII")
        self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
            self.properties, MHC_I), "Recognition_Potential_mhcI")
        self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
            self.properties, MHC_II), "Recognition_Potential_mhcII")

        # T cell predictor
        self.tcellpredict.main(self.properties)
        self.add_features(self.tcellpredict.TcellPrdictionScore, "Tcell_predictor_score")
        self.add_features(self.tcellpredict.TcellPrdictionScore_9merPred, "Tcell_predictor_score_9mersPredict")

        # DAI with affinity values
        self.add_features(FeatureLiterature.dai(self.properties, MHC_I, multiple_binding=False, affinity=True),
                          "DAI_affinity")

        # DAI wiht rank scores by netmhcpan4
        self.add_features(
            FeatureLiterature.dai(self.properties, MHC_I, multiple_binding=False, affinity=False, netmhcscore=True),
            "DAI_rank_netmhcpan4")

        # Amplitude with affinity values
        self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            self.properties, MHC_I, False, True), "Amplitude_mhcI_affinity")

        # Amplitude with rank by netmhcpan4
        self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                self.properties, mhc=MHC_I, multiple_binding=False, affinity=False, netmhcscore=True),
            "Amplitude_mhcI_rank_netmhcpan4")

        # Amplitude based on best affinity prediction restricted to 9mers
        self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            self.properties, mhc=MHC_I, multiple_binding=False, nine_mer=True),
            "Amplitude_mhcI_affinity_9mer_netmhcpan4")
        self.add_features(
            self.neoantigen_fitness_calculator.wrap_pathogensimilarity(
                props=self.properties, mhc=MHC_I, fastafile=tmp_fasta, iedb=self.references.iedb, nine_mer=True),
            "Pathogensimiliarity_mhcI_9mer")
        self.add_features(
            self.neoantigen_fitness_calculator.wrap_pathogensimilarity(
                props=self.properties, mhc=MHC_I, fastafile=tmp_fasta, iedb=self.references.iedb, affinity=True),
            "Pathogensimiliarity_mhcI_affinity_nmers")

        # recogntion potential with amplitude by affinity and netmhcpan4 score
        self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
            self.properties, MHC_I, affinity=True), "Recognition_Potential_mhcI_affinity")
        self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
            self.properties, MHC_I, affinity=False, netmhcscore=True), "Recognition_Potential_mhcI_rank_netmhcpan4")

        # recogntion potential with amplitude by affinity and only 9mers considered --> value as published!!
        self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
            self.properties, MHC_I, nine_mer=True), "Recognition_Potential_mhcI_9mer_affinity")
        self.add_features(FeatureLiterature.classify_adn_cdn(self.properties, mhc=MHC_I, category="CDN"), "CDN_mhcI")
        self.add_features(FeatureLiterature.classify_adn_cdn(self.properties, mhc=MHC_II, category="CDN"), "CDN_mhcII")
        self.add_features(FeatureLiterature.classify_adn_cdn(self.properties, mhc=MHC_I, category="ADN"), "ADN_mhcI")
        self.add_features(FeatureLiterature.classify_adn_cdn(self.properties, mhc=MHC_II, category="ADN"), "ADN_mhcII")

        # netMHCIIpan predictions
        self.predII.main(self.properties, patient_hlaII, set_available_mhcII)

        # netmhcpan4 MUT scores
        self.add_features(self.predII.best_mhcII_pan_score, "best%Rank_netmhcIIpan")
        self.add_features(self.predII.best_mhcII_pan_epitope, "best_epitope_netmhcIIpan")
        self.add_features(self.predII.best_mhcII_pan_allele, "bestHLA_allele_netmhcIIpan")

        # netmhcpan4 mut affinity
        self.add_features(self.predII.best_mhcII_pan_affinity, "best_affinity_netmhcIIpan")
        self.add_features(self.predII.best_mhcII_pan_affinity_epitope, "best_affinity_epitope_netmhcIIpan")
        self.add_features(self.predII.best_mhcII_pan_affinity_allele, "bestHLA_allele_affinity_netmhcIIpan")

        # multiplexed representation MUT MHC II
        for sc, mn in zip(self.predII.MHCII_score_all_epitopes, self.predII.mean_type):
            self.add_features(sc, "MB_score_MHCII_all_epitopes_" + mn)
        for sc, mn in zip(self.predII.MHCII_score_top10, self.predII.mean_type):
            self.add_features(sc, "MB_score_MHCII_top10_" + mn)
        for sc, mn in zip(self.predII.MHCII_score_best_per_alelle, self.predII.mean_type):
            self.add_features(sc, "MB_score_MHCII_best_per_alelle_" + mn)

        # rename MB_score_best_per_alelle_harmonic to PHBR (described in Marty et al)
        self.properties["PHBR-II"] = self.properties.pop("MB_score_MHCII_best_per_alelle_harmonic")
        self.add_features(self.predII.MHCII_epitope_scores, "MB_mhcII_epitope_scores")
        self.add_features(self.predII.MHCII_epitope_seqs, "MB_mhcII_epitope_sequences")
        self.add_features(self.predII.MHCII_epitope_alleles, "MB_mhcII_alleles")
        self.add_features(self.predII.MHCII_number_strong_binders, "MB_number_pep_MHCIIscore<2")
        self.add_features(self.predII.MHCII_number_weak_binders, "MB_number_pep_MHCIIscore<10")

        # netmhcIIpan WT scores
        self.add_features(self.predII.best_mhcII_pan_score_WT, "best%Rank_netmhcIIpan_WT")
        self.add_features(self.predII.best_mhcII_pan_epitope_WT, "best_epitope_netmhcIIpan_WT")
        self.add_features(self.predII.best_mhcII_pan_allele_WT, "bestHLA_allele_netmhcIIpan_Wt")

        # netmhcIIpan wt affinity
        self.add_features(self.predII.best_mhcII_affinity_WT, "best_affinity_netmhcIIpan_WT")
        self.add_features(self.predII.best_mhcII_affinity_epitope_WT, "best_affinity_epitope_netmhcIIpan_WT")
        self.add_features(self.predII.best_mhcII_affinity_allele_WT, "bestHLA_allele_affinity_netmhcIIpan_WT")

        # multiplexed representation WT MHC II
        for sc, mn in zip(self.predII.MHCII_score_all_epitopes_WT, self.predII.mean_type):
            self.add_features(sc, "MB_score_MHCII_all_epitopes_WT_" + mn)
        for sc, mn in zip(self.predII.MHCII_score_top10_WT, self.predII.mean_type):
            self.add_features(sc, "MB_score_MHCII_top10_WT_" + mn)
        for sc, mn in zip(self.predII.MHCII_score_best_per_alelle_WT, self.predII.mean_type):
            self.add_features(sc, "MB_score_MHCII_best_per_alelle_WT_" + mn)

        # rename MB_score_best_per_alelle_harmonic to PHBR (described in Marty et al)
        if "MB_score_MHCII_best_per_alelle_WT_harmonic" in self.properties:
            self.properties["PHBR-II_WT"] = self.properties.pop("MB_score_MHCII_best_per_alelle_WT_harmonic")
        self.add_features(self.predII.MHCII_epitope_scores_WT, "MB_mhcII_epitope_scores_WT")
        self.add_features(self.predII.MHCII_epitope_seqs_WT, "MB_mhcII_epitope_sequences_WT")
        self.add_features(self.predII.MHCII_epitope_alleles_WT, "MB_mhcII_alleles_WT")
        self.add_features(self.predII.MHCII_number_strong_binders_WT, "MB_number_pep_MHCIIscore<2_WT")
        self.add_features(self.predII.MHCII_number_weak_binders_WT, "MB_number_pep_MHCIIscore<10_WT")

        # dai mhc II affinity
        self.add_features(FeatureLiterature.dai(self.properties, MHC_II, multiple_binding=False, affinity=True),
                          "DAI_mhcII_affinity")
        self.add_features(
            FeatureLiterature.dai(self.properties, MHC_II, multiple_binding=False, affinity=True, affin_filtering=True),
            "DAI_mhcII_affinity_aff_filtered")

        # dai mhc II netMHCIIpan score
        self.add_features(FeatureLiterature.dai(self.properties, MHC_II, multiple_binding=False, affinity=False),
                          "DAI_mhcII_netmhcIIpan")

        # dai multiple binding mhc II
        self.add_features(FeatureLiterature.dai(self.properties, MHC_II, multiple_binding=True), "DAI_mhcII_MB")

        # difference number of binders
        self.add_features(FeatureLiterature.diff_number_binders(self.properties, mhc=MHC_II, threshold="2"),
                          "Diff_numb_epis_mhcII<2")
        self.add_features(FeatureLiterature.diff_number_binders(self.properties, mhc=MHC_II, threshold="10"),
                          "Diff_numb_epis_mhcII<10")
        self.add_features(FeatureLiterature.ratio_number_binders(self.properties, mhc=MHC_II, threshold="2"),
                          "Ratio_numb_epis_mhcII<2")
        self.add_features(FeatureLiterature.ratio_number_binders(self.properties, mhc=MHC_II, threshold="10"),
                          "Ratio_numb_epis_mhcII<10")

        # amplitude affinity mhc II
        self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            self.properties, MHC_II, False, True), "Amplitude_mhcII_affinity")

        # amplitude multiple binding mhc II
        self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            self.properties, MHC_II, True, False), "Amplitude_mhcII_mb")

        # amplitude rank score mhc II
        self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
            self.properties, mhc=MHC_II, multiple_binding=False, affinity=False, netmhcscore=True),
            "Amplitude_mhcII_rank_netmhcpan4")
        logger.info("Amplitude mhc II: {}".format(self.properties["Amplitude_mhcII_rank_netmhcpan4"]))

        # priority score
        self.add_features(FeatureLiterature.number_of_mismatches(self.properties, MHC_I), "Number_of_mismatches_mhcI")
        self.add_features(FeatureLiterature.number_of_mismatches(self.properties, MHC_II), "Number_of_mismatches_mhcII")
        if "mutation_found_in_proteome" not in self.properties:
            self.add_features(FeatureLiterature.match_in_proteome(self.properties, db), "mutation_found_in_proteome")
        self.add_features(FeatureLiterature.calc_priority_score(self.properties), "Priority_score")

        # priority score using multiplexed representation score
        self.add_features(FeatureLiterature.calc_priority_score(self.properties, True), "Priority_score_MB")

        # neoag immunogenicity model
        self.add_features(self.neoag_calculator.wrapper_neoag(self.properties), "neoag_immunogencity")

        # IEDB immunogenicity only for epitopes with affinity < 500 nM (predicted with netMHCpan) --> in publications
        self.add_features(FeatureLiterature.calc_IEDB_immunogenicity(self.properties, MHC_I),
                          "IEDB_Immunogenicity_mhcI")
        self.add_features(FeatureLiterature.calc_IEDB_immunogenicity(self.properties, MHC_II),
                          "IEDB_Immunogenicity_mhcII")
        self.add_features(FeatureLiterature.calc_IEDB_immunogenicity(self.properties, MHC_I, affin_filtering=True),
                          "IEDB_Immunogenicity_mhcI_affinity_filtered")

        # MixMHCpred
        self.predpresentation.main(self.properties, patient_hlaI)
        self.add_features(self.predpresentation.all_peptides, "MixMHCpred_all_peptides")
        self.add_features(self.predpresentation.all_scores, "MixMHCpred_all_scores")
        self.add_features(self.predpresentation.all_ranks, "MixMHCpred_all_ranks")
        self.add_features(self.predpresentation.all_alleles, "MixMHCpred_all_alleles")
        self.add_features(self.predpresentation.best_peptide, "MixMHCpred_best_peptide")
        self.add_features(self.predpresentation.best_score, "MixMHCpred_best_score")
        self.add_features(self.predpresentation.best_rank, "MixMHCpred_best_rank")
        self.add_features(self.predpresentation.best_allele, "MixMHCpred_best_allele")
        self.add_features(self.predpresentation.best_peptide_wt, "MixMHCpred_best_peptide_wt")
        self.add_features(self.predpresentation.best_score_wt, "MixMHCpred_best_score_wt")
        self.add_features(self.predpresentation.best_rank_wt, "MixMHCpred_best_rank_wt")
        self.add_features(self.predpresentation.difference_score_mut_wt, "MixMHCpred_difference_score_mut_wt")

        # MixMHC2pred
        self.predpresentation2.main(self.properties, patient_hlaII, list_HLAII_MixMHC2pred)
        self.add_features(self.predpresentation2.all_peptides, "MixMHC2pred_all_peptides")
        self.add_features(self.predpresentation2.all_ranks, "MixMHC2pred_all_ranks")
        self.add_features(self.predpresentation2.all_alleles, "MixMHC2pred_all_alleles")
        self.add_features(self.predpresentation2.best_peptide, "MixMHC2pred_best_peptide")
        self.add_features(self.predpresentation2.best_rank, "MixMHC2pred_best_rank")
        self.add_features(self.predpresentation2.best_allele, "MixMHC2pred_best_allele")
        self.add_features(self.predpresentation2.best_peptide_wt, "MixMHC2pred_best_peptide_wt")
        self.add_features(self.predpresentation2.best_rank_wt, "MixMHC2pred_best_rank_wt")
        self.add_features(self.predpresentation2.difference_score_mut_wt, "MixMHC2pred_difference_rank_mut_wt")

        # dissimilarity to self-proteome

        # neoantigen fitness
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix="tmpseq", suffix=".fasta", delete=False)
        tmp_fasta = tmp_fasta_file.name
        self.add_features(self.dissimilarity_calculator.calculate_dissimilarity(
            self.properties, tmp_fasta, self.references), "dissimilarity")
        self.add_features(self.dissimilarity_calculator.calculate_dissimilarity(
            self.properties, tmp_fasta, self.references, filter_binder=True), "dissimilarity_filter500")

        # vaxrank
        vaxrankscore = vaxrank.VaxRank()
        vaxrankscore.main(self.properties)
        self.add_features(vaxrankscore.total_binding_score, "vaxrank_binding_score")
        self.add_features(vaxrankscore.ranking_score, "vaxrank_total_score")

        return self.properties
