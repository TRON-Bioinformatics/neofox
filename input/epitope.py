#!/usr/bin/env python

from logzero import logger

from input import FeatureLiterature
from input import MHC_I, MHC_II
from input.MixMHCpred.mixmhcpred import MixMHCpred
from input.MixMHCpred.mixmhc2pred import MixMhc2Pred
from input.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from input.dissimilarity_garnish.dissimilaritycalculator import DissimilarityCalculator
from input.helpers import properties_manager
from input.neoag.neoag_gbm_model import NeoagCalculator
from input.neoantigen_fitness.neoantigen_fitness import NeoantigenFitnessCalculator
from input.netmhcpan4.combine_netmhcIIpan_pred_multiple_binders import BestAndMultipleBinderMhcII
from input.netmhcpan4.combine_netmhcpan_pred_multiple_binders import BestAndMultipleBinder
from input.new_features import amino_acid_frequency_scores as freq_score, differential_expression
from input.self_similarity import self_similarity
from input.vaxrank import vaxrank


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
            self.tcell_predictor = TcellPrediction(references=self.references)

        def init_properties(self, col_nam, prop_list):
            """Initiates epitope property storage in a dictionary
            """
            properties = {}
            for nam, char in zip(col_nam, prop_list):
                properties[nam] = char
            return properties

        def add_features(self, new_feature, new_feature_nam):
            """Adds new features to already present epitope properties, stored in form of a dictioninary
            """
            self.properties[new_feature_nam] = new_feature

        def write_to_file(self):
            print(";".join([self.properties[key] for key in self.properties]))

        def main(self, col_nam, prop_list, db, ref_dat, aa_freq_dict, nmer_freq_dict, aaindex1_dict, aaindex2_dict,
                 set_available_mhc, set_available_mhcII, patient_hlaI, patient_hlaII, tumour_content, rna_avail):
            """ Calculate new epitope features and add to dictonary that stores all properties
            """
            self.properties = self.init_properties(col_nam, prop_list)
            xmer_wt = self.properties["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
            xmer_mut = self.properties["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
            logger.info(xmer_mut)

            gene = properties_manager.get_gene(properties=self.properties)
            vaf_tumor = self.properties.get("VAF_in_tumor", "NA")
            vaf_rna = self.properties.get("VAF_in_RNA", "NA")
            transcript_expr = self.properties["transcript_expression"]
            alleles = properties_manager.get_hla_allele(self.properties, patient_hlaI)
            alleles_hlaII = properties_manager.get_hla_allele(self.properties, patient_hlaII)

            #TODO: remove --> IEDB
            wild_type_mhci, mutation_mhci = properties_manager.get_epitopes(
                properties=self.properties, mhc=MHC_I)
            wild_type_mhcii, mutation_mhcii = properties_manager.get_epitopes(
                properties=self.properties, mhc=MHC_II)
            '''
            self.add_features(self_similarity.position_of_mutation_epitope(
                wild_type=wild_type_mhci, mutation=mutation_mhci), "pos_MUT_MHCI")
            self.add_features(self_similarity.position_of_mutation_epitope(
                wild_type=wild_type_mhcii, mutation=mutation_mhcii), "pos_MUT_MHCII")
            self.add_features(self_similarity.position_in_anchor_position(
                position_mhci=self.properties["pos_MUT_MHCI"],
                peptide_length=self.properties["MHC_I_peptide_length_.best_prediction."]
            ), "Mutation_in_anchor")
            self.add_differential_agretopicity_index_features(mutation_mhci, mutation_mhcii, wild_type_mhci,
                                                              wild_type_mhcii)
            '''

            # MHC binding independent features
            self.add_expression_features(rna_avail, tumour_content)
            self.add_differential_expression_features(gene, ref_dat)
            self.add_aminoacid_frequency_features(aa_freq_dict, mutation_mhci, nmer_freq_dict)
            self.add_aminoacid_index_features(aaindex1_dict, aaindex2_dict)
            self.add_provean_score_features()

            # HLA I predictions: NetMHCpan
            self.pred.main(xmer_mut=xmer_mut, xmer_wt=xmer_wt, alleles=alleles, set_available_mhc=set_available_mhc)
            self.add_netmhcpan4_features()
            self.add_netmhcpan4_WT_features()
            self.add_multiple_binding_features()
            self.add_multiple_binding_numdiff()
            wild_type_netmhcpan4, mutation_netmhcpan4 = properties_manager.get_netmhcpan4_epitopes(
                properties=self.properties)
            wild_type_netmhcpan4_9mer, mutation_netmhcpan4_9mer = properties_manager.get_netmhcpan4_epitopes(
                properties=self.properties, nine_mer=True)
            wild_type_netmhcpan4_affinity, mutation_netmhcpan4_affinity = properties_manager.get_scores_netmhcpan4_affinity(
                properties=self.properties, mhc=MHC_I)
            wild_type_netmhcpan4_rank, mutation_netmhcpan4_rank = properties_manager \
                .get_scores_netmhcpan4_ranks(properties=self.properties, mhc=MHC_I)
            wild_type_multiple_binding_score, mutation_multiple_binding_score = properties_manager. \
                get_scores_multiple_binding(self.properties, mhc=MHC_I)

            self.add_multiple_binding_scorediff(mut_score=mutation_multiple_binding_score, wt_score=wild_type_multiple_binding_score)
            # position of mutation
            self.add_position_mutation(epi_wt=wild_type_netmhcpan4, epi_mut=mutation_netmhcpan4,
                                       epi_wt_9mer=wild_type_netmhcpan4_9mer, epi_mut_9mer=mutation_netmhcpan4_9mer)
            # mutation in anchor
            self.add_mutation_in_anchor()

            #TODO: change to netmhcpan
            self.add_self_similarity_features(mutation_mhci, mutation_mhcii, wild_type_mhci, wild_type_mhcii)
            # TODO: remove --> IEDB
            self.add_neoantigen_fitness_features(mutation_mhci, mutation_mhcii)

            self.add_tcell_predictor_features(gene)

            #TODO:summarise into function
            # DAI with affinity values
            self.add_features(
                FeatureLiterature.dai(score_mutation=mutation_netmhcpan4_affinity, score_wild_type=wild_type_netmhcpan4_affinity),
                "DAI_affinity")
            # DAI wiht rank scores by netmhcpan4
            self.add_features(
                FeatureLiterature.dai(score_mutation=mutation_netmhcpan4_rank, score_wild_type=wild_type_netmhcpan4_rank),
                "DAI_rank_netmhcpan4")

            #TODO:summarise into function
            # Amplitude with affinity values
            self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=mutation_netmhcpan4_affinity, score_wild_type=wild_type_netmhcpan4_affinity,
                apply_correction=True), "Amplitude_mhcI_affinity")

            # Amplitude with rank by netmhcpan4
            self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=mutation_netmhcpan4_rank, score_wild_type=wild_type_netmhcpan4_rank),
                "Amplitude_mhcI_rank_netmhcpan4")

            # Amplitude based on best affinity prediction restricted to 9mers
            self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=self.properties["best_affinity_netmhcpan4_9mer"],
                score_wild_type = self.properties["best_affinity_netmhcpan4_9mer_WT"],
                apply_correction=True), "Amplitude_mhcI_affinity_9mer_netmhcpan4")

            #TODO:summarise into function
            self.add_features(
                self.neoantigen_fitness_calculator.wrap_pathogen_similarity(
                    mutation=mutation_netmhcpan4_9mer, iedb=self.references.iedb),
                "Pathogensimiliarity_mhcI_9mer")
            self.add_features(
                self.neoantigen_fitness_calculator.wrap_pathogen_similarity(
                    mutation=mutation_netmhcpan4, iedb=self.references.iedb),
                "Pathogensimiliarity_mhcI_affinity_nmers")

            # recogntion potential with amplitude by affinity and netmhcpan4 score
            self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
                amplitude=self.properties["Amplitude_mhcI_affinity"],
                pathogen_similarity=self.properties["Pathogensimiliarity_mhcI_affinity_nmers"],
                mutation_in_anchor=self.properties["Mutation_in_anchor_netmhcpan"]),
                "Recognition_Potential_mhcI_affinity")
            self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
                amplitude=self.properties["Amplitude_mhcI_rank_netmhcpan4"],
                pathogen_similarity=self.properties["Pathogensimiliarity_mhcI"],
                mutation_in_anchor=self.properties["Mutation_in_anchor_netmhcpan"]),
                "Recognition_Potential_mhcI_rank_netmhcpan4")

            # recogntion potential with amplitude by affinity and only 9mers considered --> value as published!!
            self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
                amplitude=self.properties["Amplitude_mhcI_affinity_9mer_netmhcpan4"],
                pathogen_similarity=self.properties["Pathogensimiliarity_mhcI_9mer"],
                mutation_in_anchor=self.properties["Mutation_in_anchor_netmhcpan_9mer"],
                mhc_affinity_mut=float(self.properties["best_affinity_netmhcpan4_9mer"])),
                "Recognition_Potential_mhcI_9mer_affinity")

            score_mutation_mhci = self.properties["best_affinity_netmhcpan4"]
            amplitude_mhci = self.properties["Amplitude_mhcI_affinity"]
            bdg_cutoff_classical_mhci = 50
            bdg_cutoff_alternative_mhci = 5000
            amplitude_cutoff_mhci = 10
            self.add_features(FeatureLiterature.classify_adn_cdn(
                score_mutation=score_mutation_mhci, amplitude=amplitude_mhci,
                bdg_cutoff_classical=bdg_cutoff_classical_mhci, bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                amplitude_cutoff=amplitude_cutoff_mhci, category="CDN"), "CDN_mhcI")

            #TODO:change to netmhcIIpan --> prediction later
            score_mutation_mhcii = self.properties["MHC_II_score_.best_prediction."]
            amplitude_mhcii = self.properties["Amplitude_mhcII"]
            bdg_cutoff_classical_mhcii = 1
            bdg_cutoff_alternative_mhcii = 4
            amplitude_cutoff_mhcii = 4
            self.add_features(FeatureLiterature.classify_adn_cdn(
                score_mutation=score_mutation_mhcii, amplitude=amplitude_mhcii,
                bdg_cutoff_classical=bdg_cutoff_classical_mhcii, bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii,
                amplitude_cutoff=amplitude_cutoff_mhcii, category="CDN"), "CDN_mhcII")

            self.add_features(FeatureLiterature.classify_adn_cdn(
                score_mutation=score_mutation_mhci, amplitude=amplitude_mhci,
                bdg_cutoff_classical=bdg_cutoff_classical_mhci, bdg_cutoff_alternative=bdg_cutoff_alternative_mhci,
                amplitude_cutoff=amplitude_cutoff_mhci,
                category="ADN"), "ADN_mhcI")
            self.add_features(FeatureLiterature.classify_adn_cdn(
                score_mutation=score_mutation_mhcii, amplitude=amplitude_mhcii,
                bdg_cutoff_classical=bdg_cutoff_classical_mhcii, bdg_cutoff_alternative=bdg_cutoff_alternative_mhcii,
                amplitude_cutoff=amplitude_cutoff_mhcii, category="ADN"), "ADN_mhcII")



            # netMHCIIpan predictions
            self.predII.main(sequence=xmer_mut, sequence_reference=xmer_wt, alleles=alleles_hlaII,
                             set_available_mhc=set_available_mhcII)
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
            wild_type_netmhciipan4, mutation_netmhciipan4 = properties_manager.get_scores_netmhcpan4_affinity(
                properties=self.properties, mhc=MHC_II)
            self.add_features(
                FeatureLiterature.dai(score_mutation=mutation_netmhciipan4, score_wild_type=wild_type_netmhciipan4),
                "DAI_mhcII_affinity")
            self.add_features(
                FeatureLiterature.dai(score_mutation=mutation_netmhciipan4, score_wild_type=wild_type_netmhciipan4,
                                      affin_filtering=True),
                "DAI_mhcII_affinity_aff_filtered")

            # dai mhc II netMHCIIpan score
            self.add_features(FeatureLiterature.dai(score_mutation=mutation_mhcii, score_wild_type=wild_type_mhcii),
                              "DAI_mhcII_netmhcIIpan")

            # dai multiple binding mhc II
            wild_type_multiple_binding_ii, mutation_multiple_binding_ii = properties_manager. \
                get_scores_multiple_binding(self.properties, mhc=MHC_II)
            self.add_features(
                FeatureLiterature.dai(score_mutation=mutation_multiple_binding_ii,
                                      score_wild_type=wild_type_multiple_binding_ii),
                "DAI_mhcII_MB")

            # difference number of binders
            for threshold in [2, 10]:
                num_mutation = self.properties["MB_number_pep_MHCIIscore<{}".format(threshold)]
                num_wild_type = self.properties["MB_number_pep_MHCIIscore<{}_WT".format(threshold)]
                self.add_features(FeatureLiterature.diff_number_binders(
                    num_mutation=num_mutation,num_wild_type=num_wild_type), "Diff_numb_epis_mhcII<{}".format(threshold))
                self.add_features(FeatureLiterature.ratio_number_binders(
                    num_mutation=num_mutation,num_wild_type=num_wild_type), "Ratio_numb_epis_mhcII<{}".format(threshold))

            # amplitude affinity mhc II
            self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=self.properties["best_affinity_netmhcIIpan"],
                score_wild_type = self.properties["best_affinity_netmhcIIpan_WT"],
                apply_correction=True), "Amplitude_mhcII_affinity")

            # amplitude multiple binding mhc II
            self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=self.properties["MB_score_MHCII_top10_harmonic"],
                score_wild_type = self.properties["MB_score_MHCII_top10_WT_harmonic"]),
                "Amplitude_mhcII_mb")

            # amplitude rank score mhc II
            self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=self.properties["best%Rank_netmhcIIpan"],
                score_wild_type = self.properties["best%Rank_netmhcIIpan_WT"]),
                "Amplitude_mhcII_rank_netmhcpan4")
            logger.info("Amplitude mhc II: {}".format(self.properties["Amplitude_mhcII_rank_netmhcpan4"]))

            # priority score
            self.add_features(FeatureLiterature.number_of_mismatches(
                epitope_wild_type=self.properties["best_epitope_netmhcpan4_WT"],
                epitope_mutation=self.properties["best_epitope_netmhcpan4"]), "Number_of_mismatches_mhcI")
            self.add_features(FeatureLiterature.number_of_mismatches(
                epitope_wild_type=self.properties["best_epitope_netmhcIIpan_WT"],
                epitope_mutation=self.properties["best_epitope_netmhcIIpan"]), "Number_of_mismatches_mhcII")
            if "mutation_found_in_proteome" not in self.properties:
                self.add_features(FeatureLiterature.match_in_proteome(
                    sequence=self.properties["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."], db=db),
                    "mutation_found_in_proteome")


            no_mismatch = self.properties["Number_of_mismatches_mhcI"]
            score_mut = self.properties["best%Rank_netmhcpan4"]
            score_wt = self.properties["best%Rank_netmhcpan4_WT"]
            mut_in_prot = self.properties["mutation_found_in_proteome"]

            self.add_features(FeatureLiterature.calc_priority_score(
                vaf_tumor=vaf_tumor, vaf_rna=vaf_rna, transcript_expr=transcript_expr, no_mismatch=no_mismatch,
                score_mut=score_mut, score_wt=score_wt, mut_in_prot=mut_in_prot), "Priority_score")

            # priority score using multiplexed representation score
            score_mut = self.properties["MB_score_top10_harmonic"]
            score_wt = self.properties["MB_score_WT_top10_harmonic"]
            self.add_features(FeatureLiterature.calc_priority_score(
                vaf_tumor=vaf_tumor, vaf_rna=vaf_rna, transcript_expr=transcript_expr, no_mismatch=no_mismatch,
                score_mut=score_mut, score_wt=score_wt, mut_in_prot=mut_in_prot), "Priority_score_MB")

            # neoag immunogenicity model
            sample_id = properties_manager.get_patient_id(self.properties)
            mut_peptide = self.properties["best_affinity_epitope_netmhcpan4"]
            score_mut = self.properties["best_affinity_netmhcpan4"]
            ref_peptide = self.properties["best_affinity_epitope_netmhcpan4_WT"]
            peptide_variant_position = self.properties["pos_MUT_MHCI_affinity_epi"]
            self.add_features(self.neoag_calculator.wrapper_neoag(
                sample_id=sample_id, mut_peptide=mut_peptide, score_mut=score_mut, ref_peptide=ref_peptide,
                peptide_variant_position=peptide_variant_position), "neoag_immunogencity")

            # IEDB immunogenicity only for epitopes with affinity < 500 nM (predicted with netMHCpan) --> in publications
            mhci_epitope = self.properties["best_affinity_epitope_netmhcpan4"]
            mhci_allele = self.properties["bestHLA_allele_affinity_netmhcpan4"]
            mhci_score = self.properties["best_affinity_netmhcpan4"]
            mhcii_epitope = self.properties["MHC_II_epitope_.best_prediction."]
            mhcii_allele = self.properties["MHC_II_allele_.best_prediction."]
            self.add_features(FeatureLiterature.calc_IEDB_immunogenicity(
                epitope=mhci_epitope, mhc_allele=mhci_allele, mhc_score=mhci_score), "IEDB_Immunogenicity_mhcI")
            self.add_features(FeatureLiterature.calc_IEDB_immunogenicity(
                epitope=mhcii_epitope, mhc_allele=mhcii_allele, mhc_score=None), "IEDB_Immunogenicity_mhcII")
            self.add_features(FeatureLiterature.calc_IEDB_immunogenicity(
                epitope=mhci_epitope, mhc_allele=mhci_allele, mhc_score=mhci_score, affin_filtering=True),
                              "IEDB_Immunogenicity_mhcI_affinity_filtered")

            self.add_mix_mhc_pred_features(xmer_wt=xmer_wt, xmer_mut=xmer_mut, patient_hlaI=patient_hlaI)

            self.add_mix_mhc2_pred_features(xmer_mut=xmer_mut, xmer_wt=xmer_wt, patient_hlaII=patient_hlaII)

            # dissimilarity to self-proteome

            # neoantigen fitness
            mhc_mutation = self.properties["best_affinity_epitope_netmhcpan4"]
            mhc_affinity = self.properties["best_affinity_netmhcpan4"]
            self.add_features(self.dissimilarity_calculator.calculate_dissimilarity(
                mhc_mutation=mhc_mutation, mhc_affinity=mhc_affinity, references=self.references),
                "dissimilarity")
            self.add_features(self.dissimilarity_calculator.calculate_dissimilarity(
                mhc_mutation=mhc_mutation, mhc_affinity=mhc_affinity, references=self.references,
                filter_binder=True), "dissimilarity_filter500")

            self.add_vax_rank_features()

            return self.properties

        def add_vax_rank_features(self):
            # vaxrank
            vaxrankscore = vaxrank.VaxRank()
            vaxrankscore.main(mutation_scores=self.properties["MB_affinities"],
                              expression_score=self.properties["Expression_Mutated_Transcript"])
            self.add_features(vaxrankscore.total_binding_score, "vaxrank_binding_score")
            self.add_features(vaxrankscore.ranking_score, "vaxrank_total_score")

        def add_mix_mhc2_pred_features(self, xmer_wt, xmer_mut, patient_hlaII):
            # MixMHC2pred
            #TODO:remove allele grep and pass as argument
            alleles = properties_manager.get_hla_allele(self.properties, patient_hlaII)
            self.predpresentation2.main(alleles=alleles, xmer_wt=xmer_wt, xmer_mut=xmer_mut)
            self.add_features(self.predpresentation2.all_peptides, "MixMHC2pred_all_peptides")
            self.add_features(self.predpresentation2.all_ranks, "MixMHC2pred_all_ranks")
            self.add_features(self.predpresentation2.all_alleles, "MixMHC2pred_all_alleles")
            self.add_features(self.predpresentation2.best_peptide, "MixMHC2pred_best_peptide")
            self.add_features(self.predpresentation2.best_rank, "MixMHC2pred_best_rank")
            self.add_features(self.predpresentation2.best_allele, "MixMHC2pred_best_allele")
            self.add_features(self.predpresentation2.best_peptide_wt, "MixMHC2pred_best_peptide_wt")
            self.add_features(self.predpresentation2.best_rank_wt, "MixMHC2pred_best_rank_wt")
            self.add_features(self.predpresentation2.difference_score_mut_wt, "MixMHC2pred_difference_rank_mut_wt")

        def add_mix_mhc_pred_features(self, xmer_wt, xmer_mut, patient_hlaI):
            # MixMHCpred
            #TODO:remove allele grep and pass as argument
            alleles = properties_manager.get_hla_allele(self.properties, patient_hlaI)
            self.predpresentation.main(xmer_wt=xmer_wt, xmer_mut=xmer_mut, alleles=alleles)
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

        def add_tcell_predictor_features(self, gene):
            # T cell predictor
            substitution = properties_manager.get_substitution(properties=self.properties)
            #TODO:remove --> IEDB
            epitope = self.properties["MHC_I_epitope_.best_prediction."]
            score = self.properties["MHC_I_score_.best_prediction."]
            self.add_features(self.tcell_predictor.calculate_tcell_predictor_score(
                gene=gene, substitution=substitution, epitope=epitope, score=score),
                "Tcell_predictor_score")
            epitope = self.properties["best_affinity_epitope_netmhcpan4_9mer"]
            score = self.properties["best_affinity_netmhcpan4_9mer"]
            self.add_features(self.tcell_predictor.calculate_tcell_predictor_score(
                gene=gene, substitution=substitution, epitope=epitope, score=score, threshold=500),
                "Tcell_predictor_score_9mersPredict")

        # TODO: remove this function
        def add_neoantigen_fitness_features(self, mutation_mhci, mutation_mhcii):

            # neoantigen fitness
            self.add_features(
                self.neoantigen_fitness_calculator.wrap_pathogen_similarity(
                    mutation=mutation_mhci, iedb=self.references.iedb),
                "Pathogensimiliarity_mhcI")
            self.add_features(
                self.neoantigen_fitness_calculator.wrap_pathogen_similarity(
                    mutation=mutation_mhcii, iedb=self.references.iedb),
                "Pathogensimiliarity_mhcII")

            score_mutation_mhci = self.properties["MHC_I_score_.best_prediction."]
            score_wild_type_mhci = self.properties["MHC_I_score_.WT."]
            self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=score_mutation_mhci, score_wild_type=score_wild_type_mhci), "Amplitude_mhcI")
            score_mutation_mhcii = self.properties["MHC_II_score_.best_prediction."]
            score_wild_type_mhcii = self.properties["MHC_II_score_.WT."]
            self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=score_mutation_mhcii, score_wild_type=score_wild_type_mhcii), "Amplitude_mhcII")

            # TODO: Franziska, please, review that this is the right value for mutation_in_anchor,
            #  it is possible that there was a bug here and that it should use "Mutation_in_anchor"
            self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
                amplitude=self.properties["Amplitude_mhcI"],
                pathogen_similarity=self.properties["Pathogensimiliarity_mhcI"],
                mutation_in_anchor=self.properties["Mutation_in_anchor_netmhcpan"]),
                "Recognition_Potential_mhcI")
            self.add_features(self.neoantigen_fitness_calculator.calculate_recognition_potential(
                amplitude=self.properties["Amplitude_mhcII"],
                pathogen_similarity=self.properties["Pathogensimiliarity_mhcII"],
                mutation_in_anchor=self.properties["Mutation_in_anchor_netmhcpan"]),
                "Recognition_Potential_mhcII")

        def add_self_similarity_features(self, mutation_mhci, mutation_mhcii, wild_type_mhci, wild_type_mhcii):
            # selfsimilarity
            self.add_features(self_similarity.get_self_similarity(mutation=mutation_mhci, wild_type=wild_type_mhci),
                              "Selfsimilarity_mhcI")
            self.add_features(self_similarity.get_self_similarity(
                wild_type=wild_type_mhcii, mutation=mutation_mhcii), "Selfsimilarity_mhcII")
            self.add_features(self_similarity.is_improved_binder(
                score_mutation=self.properties["best%Rank_netmhcpan4"],
                score_wild_type=self.properties["best%Rank_netmhcpan4_WT"]
            ), "ImprovedBinding_mhcI")
            self.add_features(self_similarity.is_improved_binder(
                # TODO: conversion from float representation needs to be changed
                score_mutation=self.properties["MHC_II_score_.best_prediction."],
                score_wild_type=self.properties["MHC_II_score_.WT."]
            ), "ImprovedBinding_mhcII")
            self.add_features(self_similarity.self_similarity_of_conserved_binder_only(
                has_conserved_binder=self.properties["ImprovedBinding_mhcI"],
                similarity=self.properties["Selfsimilarity_mhcI"]),
                "Selfsimilarity_mhcI_conserved_binder")

        def add_netmhcpan4_features(self):
            """
           returns netMHCpan affinity and rank scores of mutated epitope
            """
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

        def add_netmhcpan4_WT_features(self):
            """
            returns netMHCpan affinity and rank scores of WT epitope
            """
            # netmhcpan4 WT best affinity
            self.add_features(self.pred.best4_affinity_WT, "best_affinity_netmhcpan4_WT")
            self.add_features(self.pred.best4_affinity_epitope_WT, "best_affinity_epitope_netmhcpan4_WT")
            self.add_features(self.pred.best4_affinity_allele_WT, "bestHLA_allele_affinity_netmhcpan4_WT")
            # netmhcpan4 WT rank score
            self.add_features(self.pred.best4_mhc_score_WT, "best%Rank_netmhcpan4_WT")
            self.add_features(self.pred.best4_mhc_epitope_WT, "best_epitope_netmhcpan4_WT")
            self.add_features(self.pred.best4_mhc_allele_WT, "bestHLA_allele_netmhcpan4_WT")
            # netMHCpan WT best 9mer score
            self.add_features(self.pred.mhcI_score_9mer_WT, "best%Rank_netmhcpan4_9mer_WT")
            self.add_features(self.pred.mhcI_score_epitope_9mer_WT, "best_epitope_netmhcpan4_9mer_WT")
            self.add_features(self.pred.mhcI_score_allele_9mer_WT, "bestHLA_allele_netmhcpan4_9mer_Wt")
            # netmhcpan4 WT best 9mer affinity
            self.add_features(self.pred.mhcI_affinity_9mer_WT, "best_affinity_netmhcpan4_9mer_WT")
            self.add_features(self.pred.mhcI_affinity_allele_9mer_WT, "bestHLA_allele_affinity_netmhcpan4_9mer_WT")
            self.add_features(self.pred.mhcI_affinity_epitope_9mer_WT, "best_affinity_epitope_netmhcpan4_9mer_WT")

        def add_position_mutation(self, epi_wt, epi_mut, epi_wt_9mer, epi_mut_9mer):
            """
            returns position of mutation for best affinity epitope across all lengths and 9mer
            :return:
            """
            self.add_features(self_similarity.position_of_mutation_epitope(
                wild_type=epi_wt, mutation=epi_mut), "pos_MUT_MHCI_affinity_epi")
            self.add_features(self_similarity.position_of_mutation_epitope(
                wild_type=epi_wt_9mer, mutation=epi_mut_9mer),
                "pos_MUT_MHCI_affinity_epi_9mer")

        def add_mutation_in_anchor(self):
            """
            returns if mutation is in anchor position for best affinity epitope over all lengths and best 9mer affinity
            :return:
            """
            self.add_features(self_similarity.position_in_anchor_position(
                position_mhci=self.properties["pos_MUT_MHCI_affinity_epi"],
                peptide_length=self.properties["best_epitope_netmhcpan4"]),
                "Mutation_in_anchor_netmhcpan")
            self.add_features(self_similarity.position_in_anchor_position(
                position_mhci=self.properties["pos_MUT_MHCI_affinity_epi_9mer"],
                peptide_length=9),
                "Mutation_in_anchor_netmhcpan_9mer")

        def add_DAI_mhcI(self, aff_wt, aff_mut, sc_wt, sc_mut):
            """
            returns DAI based on affinity and based on rank score
            :return:
            """
            # DAI with affinity values
            self.add_features(
                FeatureLiterature.dai(score_mutation=aff_mut,
                                      score_wild_type=aff_wt, affin_filtering=True),
                "DAI_affinity_filtered")
            self.add_features(
                FeatureLiterature.dai(score_mutation=aff_mut,
                                      score_wild_type=aff_wt),
                "DAI_affinity")
            # DAI wiht rank scores by netmhcpan4
            self.add_features(
                FeatureLiterature.dai(score_mutation=sc_mut,
                                      score_wild_type=sc_wt),
                "DAI_rank_netmhcpan4")

        def add_multiple_binding_features(self):
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

        def add_multiple_binding_numdiff(self):
            """
            returns difference and ratio of # epitopes with rank scores < 1 or 2 for mutant and wt sequence
            """
            for threshold in [1, 2]:
                num_mutation = self.properties["MB_number_pep_MHCscore<{}".format(threshold)]
                num_wild_type = self.properties["MB_number_pep_WT_MHCscore<{}".format(threshold)]
                self.add_features(FeatureLiterature.diff_number_binders(
                    num_mutation=num_mutation,num_wild_type=num_wild_type), "Diff_numb_epis_<{}".format(threshold))
                self.add_features(FeatureLiterature.ratio_number_binders(
                    num_mutation=num_mutation, num_wild_type=num_wild_type), "Ratio_numb_epis_<{}".format(threshold))

        def add_multiple_binding_scorediff(self, mut_score, wt_score):
            """
            returns DAI and amplitude with multiple binding score
            """
            self.add_features(self.neoantigen_fitness_calculator.calculate_amplitude_mhc(
                score_mutation=mut_score, score_wild_type=wt_score),
                "Amplitude_mhcI_MB")
            self.add_features(
                FeatureLiterature.dai(score_mutation=mut_score, score_wild_type=wt_score),
                "DAI_mhcI_MB")


        def add_provean_score_features(self):
            # PROVEAN score
            ucsc_id = self.provean_annotator.build_ucsc_id_plus_position(
                substitution=self.properties["substitution"], ucsc_id=self.properties["UCSC_transcript"])
            self.add_features(ucsc_id, "UCSC_ID_position")
            self.add_features(self.provean_annotator.get_provean_annotation(
                mutated_aminoacid=self.properties['MUT_AA'], ucsc_id_position=ucsc_id),
                "PROVEAN_score")

        def add_aminoacid_index_features(self, aaindex1_dict, aaindex2_dict):
            # amino acid index
            mutation_aminoacid = self.properties["MUT_AA"]
            wild_type_aminoacid = self.properties["WT_AA"]
            for k in aaindex1_dict:
                self.add_features(aaindex1_dict[k].get(wild_type_aminoacid, "NA"), "{}_{}".format(k, "wt"))
                self.add_features(aaindex1_dict[k].get(mutation_aminoacid, "NA"), "{}_{}".format(k, "mut"))
            for k in aaindex2_dict:
                self.add_features(aaindex2_dict[k].get(wild_type_aminoacid, {}).get(mutation_aminoacid, "NA"), k)

        def add_aminoacid_frequency_features(self, aa_freq_dict, mutation_mhci, nmer_freq_dict):
            # amino acid frequency
            substitution = self.properties["substitution"]
            mutated_aminoacid = FeatureLiterature.wt_mut_aa(substitution=substitution, mut="mut")
            self.add_features(mutated_aminoacid, "MUT_AA")
            self.add_features(FeatureLiterature.wt_mut_aa(substitution=substitution, mut="wt"), "WT_AA")
            self.add_features(freq_score.freq_aa(mutated_aminoacid=mutated_aminoacid, dict_freq=aa_freq_dict),
                              "Frequency_mutated_AA")
            self.add_features(freq_score.freq_prod_4mer(mutation=mutation_mhci, dict_freq=aa_freq_dict),
                              "Product_Frequency_4mer")
            self.add_features(freq_score.freq_4mer(mutation=mutation_mhci, dict_freq=nmer_freq_dict), "Frequency_of_4mer")

        def add_expression_features(self, rna_avail, tumour_content):
            # expression
            transcript_expression = self.properties["transcript_expression"]
            patient_id = properties_manager.get_patient_id(self.properties)
            vaf_tumor = self.properties["VAF_in_tumor"]
            vaf_rna = vaf_tumor if rna_avail.get(patient_id, "False") == "False" else \
                self.properties.get("VAF_in_RNA", vaf_tumor)

            self.add_features(FeatureLiterature.rna_expression_mutation(
                transcript_expression=transcript_expression, vaf_rna=vaf_rna), "Expression_Mutated_Transcript")

            expression_mutated_transcript = self.properties.get("Expression_Mutated_Transcript")
            self.add_features(FeatureLiterature.expression_mutation_tc(
                transcript_expression=expression_mutated_transcript, patient_id=patient_id, tumour_content_dict=tumour_content),
                              "Expression_Mutated_Transcript_tumor_content")

        def add_differential_agretopicity_index_features(self, mutation_mhci, mutation_mhcii, wild_type_mhci,
                                                         wild_type_mhcii):
            # differential agretopicity index
            self.add_features(
                FeatureLiterature.dai(score_mutation=mutation_mhci, score_wild_type=wild_type_mhci), "DAI_mhcI")
            self.add_features(
                FeatureLiterature.dai(score_mutation=mutation_mhcii, score_wild_type=wild_type_mhcii), "DAI_mhcII")


        def add_differential_expression_features(self, gene, ref_dat):
            # differential expression
            expression_tumor = self.properties["transcript_expression"]
            expression_reference = differential_expression.add_rna_reference(gene, ref_dat, 0)
            expression_reference_sum = differential_expression.add_rna_reference(gene, ref_dat, 2)
            expression_reference_sd = differential_expression.add_rna_reference(gene, ref_dat, 1)
            self.add_features(expression_reference, "mean_ref_expression")
            self.add_features(expression_reference_sd, "sd_ref_expression")
            self.add_features(expression_reference_sum, "sum_ref_expression")
            self.add_features(differential_expression.fold_change(
                expression_tumor=expression_tumor, expression_reference=expression_reference), "log2_fc_tumour_ref")
            self.add_features(differential_expression.percentile_calc(
                expression_tumor=expression_tumor, expression_reference_sum=expression_reference_sum),
                "percentile_tumour_ref")
            self.add_features(differential_expression.pepper_calc(
                expression_tumor=expression_tumor, expression_reference=expression_reference,
                expression_reference_sd=expression_reference_sd), "DE_pepper")
