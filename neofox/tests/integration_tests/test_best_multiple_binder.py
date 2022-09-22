#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
from logzero import logger
from unittest import TestCase
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.factories import MhcFactory
from neofox.model.mhc_parser import MhcParser
import neofox.tests.integration_tests.integration_test_tools as integration_test_tools
from neofox.helpers.runner import Runner
from neofox.MHC_predictors.netmhcpan.combine_netmhcpan_pred_multiple_binders import (
    BestAndMultipleBinder,
)
from neofox.MHC_predictors.netmhcpan.netmhcpan_prediction import NetMhcPanPredictor
from neofox.MHC_predictors.netmhcpan.combine_netmhcIIpan_pred_multiple_binders import (
    BestAndMultipleBinderMhcII,
)
from neofox.MHC_predictors.netmhcpan.netmhcIIpan_prediction import NetMhcIIPanPredictor
from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.tests.tools import get_neoantigen


class TestBestMultipleBinder(TestCase):
    def setUp(self):
        references, self.configuration = integration_test_tools.load_references()
        self.runner = Runner()
        self.available_alleles_mhc1 = (
            references.get_available_alleles().get_available_mhc_i()
        )
        self.available_alleles_mhc2 = (
            references.get_available_alleles().get_available_mhc_ii()
        )
        self.hla_database = references.get_mhc_database()
        self.mhc_parser = MhcParser.get_mhc_parser(self.hla_database)
        self.test_mhc_one = integration_test_tools.get_hla_one_test(self.hla_database)
        self.test_mhc_two = integration_test_tools.get_hla_two_test(self.hla_database)
        self.uniprot = Uniprot(references.uniprot_pickle)
        self.proteome_blastp_runner = BlastpRunner(
            runner=self.runner, configuration=self.configuration,
            database=references.get_proteome_database())

    def test_best_multiple_run(self):
        best_multiple = BestAndMultipleBinder(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        # this is some valid example neoantigen candidate sequence
        mutation = get_neoantigen(
                mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
                wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
        )
        best_multiple.run(
            neoantigen=mutation,
            mhc1_alleles_patient=self.test_mhc_one,
            mhc1_alleles_available=self.available_alleles_mhc1,
            uniprot=self.uniprot,
        )
        self.assertEqual(602.12, best_multiple.best_epitope_by_affinity.affinity_mutated)
        self.assertEqual('HLA-A*02:01', best_multiple.best_epitope_by_affinity.allele_mhc_i.name)
        self.assertEqual(0.492, best_multiple.best_epitope_by_rank.rank_mutated)
        self.assertEqual('HLA-A*02:01', best_multiple.best_epitope_by_rank.allele_mhc_i.name)
        self.assertEqual("ILVTDQTRL", best_multiple.best_epitope_by_rank.mutated_peptide)

    def test_phbr1(self):
        best_multiple = BestAndMultipleBinder(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        netmhcpan = NetMhcPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        mutation = get_neoantigen(
                mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
                wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
        )
        # all alleles = heterozygous
        available_alleles = netmhcpan.get_only_available_alleles(self.test_mhc_one, self.available_alleles_mhc1)
        predictions = netmhcpan.mhc_prediction(available_alleles, mutation.mutated_xmer)

        predicted_neoepitopes = EpitopeHelper.remove_peptides_in_proteome(
            predictions=predictions, uniprot=self.uniprot
        )
        best_epitopes_per_allele = (
            BestAndMultipleBinder.extract_best_epitope_per_alelle(
                predicted_neoepitopes, self.test_mhc_one
            )
        )
        phbr_i = best_multiple.calculate_phbr_i(best_epitopes_per_allele, self.test_mhc_one)
        self.assertIsNotNone(phbr_i)
        self.assertAlmostEqual(1.359324592015038, phbr_i)
        # one homozygous allele present
        mhc_alleles = MhcFactory.build_mhc1_alleles(
            [
                "HLA-A*24:02",
                "HLA-A*02:01",
                "HLA-B*15:01",
                "HLA-B*44:02",
                "HLA-C*05:01",
                "HLA-C*05:01",
            ], self.hla_database
        )

        predictions = netmhcpan.mhc_prediction(available_alleles, mutation.mutated_xmer)

        predicted_neoepitopes = EpitopeHelper.remove_peptides_in_proteome(
            predictions=predictions,uniprot=self.uniprot
        )
        best_epitopes_per_allele = (
            BestAndMultipleBinder.extract_best_epitope_per_alelle(
                predicted_neoepitopes, mhc_alleles
            )
        )
        phbr_i = best_multiple.calculate_phbr_i(best_epitopes_per_allele, mhc_alleles)
        self.assertIsNotNone(phbr_i)
        self.assertAlmostEqual(1.0036998409510969, phbr_i)
        # mo info for one allele
        mhc_alleles = MhcFactory.build_mhc1_alleles(
            ["HLA-A*24:02", "HLA-A*02:01", "HLA-B*15:01", "HLA-B*44:02", "HLA-C*05:01"], self.hla_database
        )

        predictions = netmhcpan.mhc_prediction(available_alleles, mutation.mutated_xmer)
        predicted_neoepitopes = EpitopeHelper.remove_peptides_in_proteome(
            predictions=predictions, uniprot=self.uniprot
        )
        best_epitopes_per_allele = (
            BestAndMultipleBinder.extract_best_epitope_per_alelle(
                predicted_neoepitopes, mhc_alleles
            )
        )
        phbr_i = best_multiple.calculate_phbr_i(best_epitopes_per_allele, mhc_alleles)
        self.assertIsNone(phbr_i)

    def test_best_multiple_mhc2_run(self):
        best_multiple = BestAndMultipleBinderMhcII(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        # this is some valid example neoantigen candidate sequence
        mutation = get_neoantigen(
                mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
                wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
        )
        best_multiple.run(
            neoantigen=mutation,
            mhc2_alleles_patient=self.test_mhc_two,
            mhc2_alleles_available=self.available_alleles_mhc2,
            uniprot=self.uniprot
        )
        logger.info(best_multiple.best_predicted_epitope_rank)
        logger.info(best_multiple.best_predicted_epitope_affinity)
        logger.info(best_multiple.phbr_ii)
        self.assertEqual(3.26, best_multiple.best_predicted_epitope_rank.rank_mutated)
        self.assertEqual(
            1103.46, best_multiple.best_predicted_epitope_affinity.affinity_mutated
        )
        self.assertEqual(
            "SQDILVTDQTRLEAT", best_multiple.best_predicted_epitope_rank.mutated_peptide
        )

    def test_phbr2(self):
        best_multiple = BestAndMultipleBinderMhcII(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        netmhc2pan = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        mutation = get_neoantigen(
            mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
            wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
        )
        # all alleles = heterozygous
        allele_combinations = netmhc2pan.generate_mhc2_alelle_combinations(self.test_mhc_two)
        patient_mhc2_isoforms = best_multiple._get_only_available_combinations(
            allele_combinations, self.available_alleles_mhc2
        )
        predictions = netmhc2pan.mhc2_prediction(
            patient_mhc2_isoforms, mutation.mutated_xmer
        )
        filtered_predictions = EpitopeHelper.remove_peptides_in_proteome(
            predictions=predictions, uniprot=self.uniprot
        )
        logger.info(filtered_predictions)
        logger.info(self.test_mhc_two)
        best_predicted_epitopes_per_alelle = (
            best_multiple.extract_best_epitope_per_mhc2_alelle(predictions=filtered_predictions, mhc_isoforms=self.test_mhc_two
            )
        )
        phbr_ii = best_multiple.calculate_phbr_ii(best_predicted_epitopes_per_alelle)
        self.assertIsNotNone(phbr_ii)
        self.assertAlmostEqual(8.895757526065129, phbr_ii)
        # mo info for one allele
        mhc2_alleles = MhcFactory.build_mhc2_alleles(
            [
                "HLA-DRB1*04:02",
                "HLA-DRB1*08:01",
                "HLA-DQA1*03:01",
                "HLA-DQA1*04:01",
                "HLA-DQB1*03:02",
                "HLA-DQB1*04:02",
                "HLA-DPA1*01:03",
                "HLA-DPA1*02:01",
                "HLA-DPB1*13:01",
                "HLA-DPB1*13:01",
            ], self.hla_database
        )

        allele_combinations = netmhc2pan.generate_mhc2_alelle_combinations(mhc2_alleles)
        patient_mhc2_isoforms = best_multiple._get_only_available_combinations(
            allele_combinations, self.available_alleles_mhc2
        )
        predictions = netmhc2pan.mhc2_prediction(
            patient_mhc2_isoforms, mutation.mutated_xmer
        )
        filtered_predictions = EpitopeHelper.remove_peptides_in_proteome(
            predictions=predictions, uniprot=self.uniprot
        )
        best_predicted_epitopes_per_alelle = (
            best_multiple.extract_best_epitope_per_mhc2_alelle(
                filtered_predictions, mhc2_alleles
            )
        )
        logger.info(best_predicted_epitopes_per_alelle)
        logger.info(len(best_predicted_epitopes_per_alelle))
        phbr_ii = best_multiple.calculate_phbr_ii(best_predicted_epitopes_per_alelle)
        self.assertIsNone(phbr_ii)

        # one allele present
        mhc2_alleles = MhcFactory.build_mhc2_alleles(
            [
                "HLA-DRB1*04:02",
                "HLA-DRB1*08:01",
                "HLA-DQA1*03:01",
                "HLA-DQA1*04:01",
                "HLA-DQB1*03:02",
                "HLA-DQB1*04:02",
                "HLA-DPA1*01:03",
                "HLA-DPA1*02:01",
                "HLA-DPB1*13:01",
            ],
            self.hla_database
        )
        allele_combinations = netmhc2pan.generate_mhc2_alelle_combinations(mhc2_alleles)
        patient_mhc2_isoforms = best_multiple._get_only_available_combinations(
            allele_combinations, self.available_alleles_mhc2
        )
        predictions = netmhc2pan.mhc2_prediction(
            patient_mhc2_isoforms, mutation.mutated_xmer
        )
        filtered_predictions = EpitopeHelper.remove_peptides_in_proteome(
            predictions=predictions, uniprot=self.uniprot
        )
        best_predicted_epitopes_per_alelle = (
            best_multiple.extract_best_epitope_per_mhc2_alelle(
                filtered_predictions, mhc2_alleles
            )
        )
        logger.info(best_predicted_epitopes_per_alelle)
        logger.info(len(best_predicted_epitopes_per_alelle))
        phbr_ii = best_multiple.calculate_phbr_ii(best_predicted_epitopes_per_alelle)
        self.assertIsNone(phbr_ii)

    def test_generator_rate(self):
        best_multiple = BestAndMultipleBinder(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        netmhcpan = NetMhcPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        mutation = get_neoantigen(
            mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
            wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
        )
        # all alleles = heterozygous
        available_alleles = netmhcpan.get_only_available_alleles(self.test_mhc_one, self.available_alleles_mhc1)
        predictions = netmhcpan.mhc_prediction(available_alleles, mutation.mutated_xmer)
        predictions_wt = netmhcpan.mhc_prediction(available_alleles, mutation.wild_type_xmer)

        predicted_neoepitopes = EpitopeHelper.remove_peptides_in_proteome(predictions=predictions, uniprot=self.uniprot)
        filtered_predictions_wt = EpitopeHelper.filter_peptides_covering_snv(
            position_of_mutation=mutation.position, predictions=predictions_wt)
        paired_predictions = EpitopeHelper.pair_predictions(
            predictions=predicted_neoepitopes, predictions_wt=filtered_predictions_wt)

        generator_rate_ADN = best_multiple.determine_number_of_alternative_binders(predictions=paired_predictions)
        generator_rate_CDN = best_multiple.determine_number_of_binders(predictions=paired_predictions, threshold=50)
        self.assertEqual(generator_rate_ADN, 0)
        self.assertEqual(generator_rate_CDN, 0)


    def test_generator_rate_mhcII(self):
        best_multiple = BestAndMultipleBinderMhcII(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        netmhc2pan = NetMhcIIPanPredictor(
            runner=self.runner, configuration=self.configuration, mhc_parser=self.mhc_parser,
            blastp_runner=self.proteome_blastp_runner
        )
        mutation = get_neoantigen(
            mutated_xmer="RTNLLAALHRSVRWRAADQGHRSAFLV",
            wild_type_xmer="RTNLLAALHRSVRRRAADQGHRSAFLV",
        )
        # all alleles = heterozygous
        allele_combinations = netmhc2pan.generate_mhc2_alelle_combinations(self.test_mhc_two)
        patient_mhc2_isoforms = best_multiple._get_only_available_combinations(
            allele_combinations, self.available_alleles_mhc2
        )
        predictions = netmhc2pan.mhc2_prediction(
            patient_mhc2_isoforms, mutation.mutated_xmer
        )

        predictions_wt = netmhc2pan.mhc2_prediction(
            patient_mhc2_isoforms, mutation.wild_type_xmer
        )

        predicted_neoepitopes = EpitopeHelper.remove_peptides_in_proteome(
            predictions=predictions, uniprot=self.uniprot
        )
        filtered_predictions_wt = EpitopeHelper.filter_peptides_covering_snv(
            position_of_mutation=mutation.position, predictions=predictions_wt
        )

        paired_predictions = EpitopeHelper.pair_predictions(
            predictions=predicted_neoepitopes, predictions_wt=filtered_predictions_wt)

        generator_rate_ADN = best_multiple.determine_number_of_alternative_binders(predictions=paired_predictions)
        generator_rate_CDN = best_multiple.determine_number_of_binders(predictions=paired_predictions)
        self.assertEqual(generator_rate_ADN, 6)
        self.assertEqual(generator_rate_CDN, 0)
