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
import unittest
from unittest import TestCase
from datetime import datetime
import pkg_resources
import orjson as json
from neofox.exceptions import NeofoxConfigurationException

from neofox import NEOFOX_MIXMHCPRED_ENV, NEOFOX_MIXMHC2PRED_ENV, NEOFOX_PRIME_ENV

import neofox.tests
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import Neoantigen, Patient
from neofox.model.factories import NOT_AVAILABLE_VALUE, PatientFactory, MhcFactory
from neofox.neofox import NeoFox
from neofox.references.references import ORGANISM_MUS_MUSCULUS
from neofox.tests.integration_tests import integration_test_tools
import pandas as pd
from logzero import logger
import os
import shutil
import math
import numpy as np
import timeit


class TestNeofox(TestCase):

    def setUp(self):
        self.references, self.configuration = integration_test_tools.load_references()
        self.references_mouse, self.configuration_mouse = integration_test_tools.load_references(
            organism=ORGANISM_MUS_MUSCULUS)
        # self.fastafile = integration_test_tools.create_temp_aminoacid_fasta_file()
        # self.runner = Runner()
        self.patient_id = "Pt29"
        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_model_realistic.txt"
        )
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_file.txt"
        )
        patients_file_mouse = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_file_mouse.txt"
        )
        self.hla_database = self.references.get_mhc_database()
        self.h2_database = self.references_mouse.get_mhc_database()
        self.patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        self.patients_mouse = ModelConverter.parse_patients_file(patients_file_mouse, self.h2_database)
        self.neoantigens = ModelConverter.parse_candidate_file(input_file)
        self.neoantigens_mouse = ModelConverter.parse_candidate_file(input_file)

    def test_neoantigens_without_gene(self):
        """"""
        neoantigens, patients = self._get_test_data()
        for n in neoantigens:
            n.gene = ""
        annotations = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        ).get_annotations()
        self.assertEqual(5, len(annotations))
        self.assertIsInstance(annotations[0], Neoantigen)
        self.assertTrue(len(annotations[0].neofox_annotations.annotations) > 10)

    def _get_test_data(self):
        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_model.txt"
        )
        data = pd.read_csv(input_file, sep="\t")
        data = data.replace({np.nan: None})
        neoantigens = ModelConverter._neoantigens_csv2objects(data)
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_file.txt"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        return neoantigens, patients

    def test_neofox(self):
        """
        This test is equivalent to the command line call of neofox

        NOTE: we will need to check the output when the calculation of resuls and printing to stdout have been decoupled
        """
        output_file = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/output_{:%Y%m%d%H%M%S}_neoantigen_candidates_annotated.tsv".format(
                datetime.now()
            ),
        )
        output_json_neoantigens = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/output_{:%Y%m%d%H%M%S}.neoantigen_candidates.json".format(
                datetime.now()
            ),
        )
        annotations = NeoFox(
            neoantigens=self.neoantigens,
            patients=self.patients,
            num_cpus=4,
        ).get_annotations()
        annotation_names = [a.name for n in annotations for a in n.neofox_annotations.annotations]

        # check it does contain any of the MixMHCpred annotations
        self.assertIn("MixMHC2pred_bestRank_peptide", annotation_names)
        self.assertIn("MixMHC2pred_bestRank_rank", annotation_names)
        self.assertIn("MixMHC2pred_bestRank_allele", annotation_names)
        self.assertIn("MixMHCpred_bestScore_peptide", annotation_names)
        self.assertIn("MixMHCpred_bestScore_score", annotation_names)
        self.assertIn("MixMHCpred_bestScore_rank", annotation_names)
        self.assertIn("MixMHCpred_bestScore_allele", annotation_names)
        # checks it does have some of the NetMHCpan annotations
        self.assertIn("NetMHCpan_bestAffinity9mer_positionMutation", annotation_names)
        self.assertIn("NetMHCIIpan_bestRank_rank", annotation_names)

        # writes output
        ModelConverter.annotations2neoantigens_table(neoantigens=annotations).to_csv(
            output_file, sep="\t", index=False)

        with open(output_json_neoantigens, "wb") as f:
            f.write(json.dumps(ModelConverter.objects2json(annotations)))

        # regression test
        self._regression_test_on_output_file(new_file=output_file)

    def test_neomouse(self):
        output_file = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/output_mouse_{:%Y%m%d%H%M%S}_neoantigen_candidates_annotated.tsv".format(
                datetime.now()
            ),
        )

        output_json_neoantigens = pkg_resources.resource_filename(
            neofox.tests.__name__,
            "resources/output_mouse_{:%Y%m%d%H%M%S}.neoantigen_candidates.json".format(
                datetime.now()
            ),
        )
        annotations = NeoFox(
            neoantigens=self.neoantigens_mouse,
            patients=self.patients_mouse,
            num_cpus=4,
            reference_folder=self.references_mouse
        ).get_annotations()
        annotation_names = [a.name for n in annotations for a in n.neofox_annotations.annotations]

        # checks it does have some of the NetMHCpan annotations
        self.assertIn("NetMHCpan_bestAffinity9mer_positionMutation", annotation_names)
        self.assertIn("NetMHCIIpan_bestRank_rank", annotation_names)

        # writes output
        ModelConverter.annotations2neoantigens_table(neoantigens=annotations).to_csv(
            output_file, sep="\t", index=False)

        with open(output_json_neoantigens, "wb") as f:
            f.write(json.dumps(ModelConverter.objects2json(annotations)))

        # regression test
        self._regression_test_on_output_file(new_file=output_file, previous_filename="resources/output_previous_mouse.txt")

    def test_neofox_only_one_neoantigen(self):
        """"""
        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_model_only_one.txt"
        )
        neoantigens = ModelConverter.parse_candidate_file(
            input_file
        )
        annotations = NeoFox(
            neoantigens=neoantigens,
            patients=self.patients,
            num_cpus=4,
        ).get_annotations()
        self.assertEqual(1, len(annotations))
        self.assertIsInstance(annotations[0], Neoantigen)
        self.assertTrue(len(annotations[0].neofox_annotations.annotations) > 10)

    def test_neofox_model_input(self):
        """"""
        neoantigens, patients = self._get_test_data()
        annotations = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        ).get_annotations()
        self.assertEqual(5, len(annotations))
        self.assertIsInstance(annotations[0], Neoantigen)
        self.assertTrue(len(annotations[0].neofox_annotations.annotations) == 82)

    def test_neofox_without_mixmhcpreds(self):
        """
        This test aims at testing neofox when MixMHCpred and MixMHC2pred are not configured. As these are optional it
        shoudl just run, but without these annotations in the output
        """
        del os.environ[NEOFOX_MIXMHCPRED_ENV]
        del os.environ[NEOFOX_MIXMHC2PRED_ENV]
        annotations = NeoFox(
            neoantigens=self.neoantigens,
            patients=self.patients,
            num_cpus=4,
        ).get_annotations()
        annotation_names = [a.name for n in annotations for a in n.neofox_annotations.annotations]
        # check it does not contain any of the MixMHCpred annotations
        self.assertNotIn("MixMHC2pred_bestRank_peptide", annotation_names)
        self.assertNotIn("MixMHC2pred_bestRank_rank", annotation_names)
        self.assertNotIn("MixMHC2pred_bestRank_allele", annotation_names)
        self.assertNotIn("MixMHCpred_bestScore_peptide", annotation_names)
        self.assertNotIn("MixMHCpred_bestScore_score", annotation_names)
        self.assertNotIn("MixMHCpred_bestScore_rank", annotation_names)
        self.assertNotIn("MixMHCpred_bestScore_allele", annotation_names)
        # checks it does have some of the NetMHCpan annotations
        self.assertIn("NetMHCpan_bestAffinity9mer_positionMutation", annotation_names)
        self.assertIn("NetMHCIIpan_bestRank_rank", annotation_names)

    def test_neofox_without_prime(self):
        """
        This test aims at testing neofox when Prime is not configured. As these are optional it
        shoudl just run, but without these annotations in the output
        """
        del os.environ[NEOFOX_PRIME_ENV]
        annotations = NeoFox(
            neoantigens=self.neoantigens[0:1],
            patients=self.patients,
            num_cpus=4,
        ).get_annotations()
        annotation_names = [a.name for n in annotations for a in n.neofox_annotations.annotations]
        # check it does not contain any of the MixMHCpred annotations
        self.assertNotIn("PRIME_best_peptide", annotation_names)
        self.assertNotIn("PRIME_best_rank", annotation_names)
        self.assertNotIn("PRIME_best_allele", annotation_names)
        # checks it does have some of the NetMHCpan annotations
        self.assertIn("NetMHCpan_bestAffinity9mer_positionMutation", annotation_names)
        self.assertIn("NetMHCIIpan_bestRank_rank", annotation_names)

    def test_neofox_with_prime_and_without_mixmhcpred(self):
        """
        This test aims at testing neofox when Prime is configured, but not MixMHCpred. As PRIME depends on
        MixMHCpred no PRIME annotations should be provided
        """
        del os.environ[NEOFOX_MIXMHCPRED_ENV]
        annotations = NeoFox(
            neoantigens=self.neoantigens[0:1],
            patients=self.patients,
            num_cpus=4,
        ).get_annotations()
        annotation_names = [a.name for n in annotations for a in n.neofox_annotations.annotations]
        # check it does not contain any of the MixMHCpred annotations
        self.assertNotIn("PRIME_best_peptide", annotation_names)
        self.assertNotIn("PRIME_best_rank", annotation_names)
        self.assertNotIn("PRIME_best_allele", annotation_names)
        # checks it does have some of the NetMHCpan annotations
        self.assertIn("NetMHCpan_bestAffinity9mer_positionMutation", annotation_names)
        self.assertIn("NetMHCIIpan_bestRank_rank", annotation_names)

    @unittest.skip
    def test_neofox_performance(self):
        def compute_annotations():
            return NeoFox(
                neoantigens=self.neoantigens,
                patients=self.patients,
                num_cpus=4,
            ).get_annotations()

        print("Average time: {}".format(timeit.timeit(compute_annotations, number=5)))

    @unittest.skip
    def test_neofox_performance_single_neoantigen(self):

        input_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_only_one.txt"
        )
        neoantigens, _ = ModelConverter.parse_candidate_file(input_file)

        def compute_annotations():
            return NeoFox(
                neoantigens=neoantigens,
                patients=self.patients,
                num_cpus=4,
            ).get_annotations()

        print("Average time: {}".format(timeit.timeit(compute_annotations, number=10)))

    def test_neofox_with_config(self):
        neoantigens, patients = self._get_test_data()
        config_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/neofox_config.txt"
        )
        try:
            NeoFox(
                neoantigens=neoantigens,
                patients=patients,
                num_cpus=4,
                configuration_file=config_file,
            )
        except NeofoxConfigurationException as e:
            assert "/neofox/testing/reference_data" in str(e)
            return
        assert False

    def test_neofox_without_mhc2(self):
        """"""
        neoantigens, patients = self._get_test_data()
        for p in patients:
            p.mhc2 = []
        annotations = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        ).get_annotations()
        self.assertEqual(5, len(annotations))
        self.assertIsInstance(annotations[0], Neoantigen)
        self.assertEqual(len(annotations[0].neofox_annotations.annotations), 63)

    def test_neofox_without_mhc1(self):
        neoantigens, patients = self._get_test_data()
        for p in patients:
            p.mhc1 = []
        annotations = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        ).get_annotations()
        self.assertEqual(5, len(annotations))
        self.assertIsInstance(annotations[0], Neoantigen)
        self.assertEqual(len(annotations[0].neofox_annotations.annotations), 39)

    def test_gene_expression_imputation(self):
        neoantigens, patients = self._get_test_data()
        for p in patients:
            p.is_rna_available = False
        neofox = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        )
        for n in neofox.neoantigens:
            self.assertIsNotNone(n.imputed_gene_expression)
            self.assertGreater(n.imputed_gene_expression, 0)

    def test_neoantigens_with_non_existing_gene(self):
        """"""
        neoantigens, patients = self._get_test_data()
        for n in neoantigens:
            n.gene = "IDONTEXIST"
        neofox = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        )
        for n in neofox.neoantigens:
            self.assertIsNone(n.imputed_gene_expression)

    def test_neoantigens_with_empty_gene(self):
        """"""
        neoantigens, patients = self._get_test_data()
        for n in neoantigens:
            n.gene = ""
        neofox = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        )
        for n in neofox.neoantigens:
            self.assertIsNone(n.imputed_gene_expression)

    def test_neoantigens_with_empty_rna_expression(self):
        """"""
        neoantigens, patients = self._get_test_data()
        for n in neoantigens:
            n.rna_expression = None
            n.patient_identifier = "patient_without_tumor_type"
        neofox = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        )
        for n in neofox.neoantigens:
            self.assertIsNone(n.rna_expression)

    def test_neoantigens_with_rna_expression(self):
        """"""
        neoantigens, patients = self._get_test_data()
        for n in neoantigens:
            n.rna_expression = 1.2
        neofox = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        )

        for n in neofox.neoantigens:
            self.assertEqual(n.rna_expression, 1.2)


    def test_patient_with_non_existing_allele_does_not_crash(self):
        """"""
        neoantigens, patients = self._get_test_data()
        for p in patients:
            # sets one MHC I allele to a non existing allele
            p.mhc1[0].alleles[0] = MhcFactory.build_mhc1_alleles(["HLA-A*99:99"], mhc_database=self.hla_database)[0].alleles[0]
        neofox = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        )
        neofox.get_annotations()

    def test_neoantigens_with_rare_aminoacids(self):
        """"""
        neoantigens, patients = self._get_test_data()
        for n in neoantigens:
            position_to_replace = int(len(n.mutated_xmer)/2)
            n.mutated_xmer = n.mutated_xmer[:position_to_replace] + "U" + \
                                      n.mutated_xmer[position_to_replace+1:]
        annotations = NeoFox(
            neoantigens=neoantigens,
            patients=patients,
            num_cpus=4,
        ).get_annotations()
        self.assertEqual(5, len(annotations))
        self.assertIsInstance(annotations[0], Neoantigen)
        self.assertTrue(len(annotations[0].neofox_annotations.annotations) > 10)
        for na in annotations:
            for a in na.neofox_annotations.annotations:
                if a.name in ["Selfsimilarity_MHCI_conserved_binder", "Tcell_predictor_score_cutoff"]:
                    self.assertEqual(a.value, NOT_AVAILABLE_VALUE)

    def test_neoantigen_without_9mer_netmhcpan_results(self):
        patient_identifier = "12345"
        neoantigen = Neoantigen(
            wild_type_xmer="HLAQHQRVHTGEKPYKCNECGKTFRQT",
            mutated_xmer="HLAQHQRVHTGEKAYKCNECGKTFRQT",
            patient_identifier=patient_identifier
        )
        patient = PatientFactory.build_patient(
            identifier=patient_identifier,
            mhc_alleles=["HLA-A*24:106", "HLA-A*02:200", "HLA-B*08:33", "HLA-B*40:94", "HLA-C*02:20", "HLA-C*07:86"],
            mhc2_alleles=["HLA-DRB1*07:14", "HLA-DRB1*04:18", "HLA-DPA1*01:05", "HLA-DPA1*03:01", "HLA-DPB1*17:01",
                          "HLA-DPB1*112:01", "HLA-DQA1*01:06", "HLA-DQA1*01:09", "HLA-DQB1*03:08", "HLA-DQB1*06:01"],
            mhc_database=self.references.get_mhc_database()
        )

        annotations = NeoFox(
            neoantigens=[neoantigen],
            patients=[patient],
            num_cpus=4,
        ).get_annotations()
        # it does not crash even though there are no best 9mers
        self.assertIsNotNone(annotations)

    def test_neoantigen_in_proteome(self):
        patient_identifier = "12345"
        neoantigen = Neoantigen(
            mutated_xmer="PKLLENLLSKGETISFLECF",
            patient_identifier=patient_identifier
        )
        patient = PatientFactory.build_patient(
            identifier=patient_identifier,
            mhc_alleles=["HLA-A*24:106", "HLA-A*02:200", "HLA-B*08:33", "HLA-B*40:94", "HLA-C*02:20", "HLA-C*07:86"],
            mhc2_alleles=["HLA-DRB1*07:14", "HLA-DRB1*04:18", "HLA-DPA1*01:05", "HLA-DPA1*03:01", "HLA-DPB1*17:01",
                          "HLA-DPB1*112:01", "HLA-DQA1*01:06", "HLA-DQA1*01:09", "HLA-DQB1*03:08", "HLA-DQB1*06:01"],
            mhc_database=self.references.get_mhc_database()
        )

        annotations = NeoFox(
            neoantigens=[neoantigen],
            patients=[patient],
            num_cpus=4,
        ).get_annotations()
        # it does not crash even though there are no best 9mers
        self.assertIsNotNone(annotations)

    def test_neoantigen_failing(self):
        patient_identifier = "12345"
        neoantigen = Neoantigen(
            wild_type_xmer="ARPDMFCLFHGKRYFPGESWHPYLEPQ",
            mutated_xmer="ARPDMFCLFHGKRHFPGESWHPYLEPQ",
            patient_identifier=patient_identifier
        )
        patient = Patient(
            identifier=patient_identifier,
            mhc1=MhcFactory.build_mhc1_alleles([
                "HLA-A*03:01", "HLA-A*29:02", "HLA-B*07:02", "HLA-B*44:03", "HLA-C*07:02", "HLA-C*16:01"],
                mhc_database=self.references.get_mhc_database()),
        )

        annotations = NeoFox(
            neoantigens=[neoantigen],
            patients=[patient],
            num_cpus=4,
        ).get_annotations()
        # it does not crash even though there are no best 9mers
        self.assertIsNotNone(annotations)

    def test_neoantigen_no_wt_failing(self):
        patient_identifier = "12345"
        neoantigen = Neoantigen(
            mutated_xmer="SPSFPLEPDDEVFTAIAKAMEEMVEDS",
            patient_identifier=patient_identifier
        )
        patient = Patient(
            identifier=patient_identifier,
            mhc1=MhcFactory.build_mhc1_alleles([
                "HLA-A*02:24", "HLA-A*36:04", "HLA-B*58:25", "HLA-B*35:102", "HLA-C*02:30", "HLA-C*07:139"],
                mhc_database=self.references.get_mhc_database()),
        )

        annotations = NeoFox(
            neoantigens=[neoantigen],
            patients=[patient],
            num_cpus=4,
        ).get_annotations()
        # it does not crash even though there are no best 9mers
        self.assertIsNotNone(annotations)

    @unittest.skip
    def test_neofox_synthetic_data(self):
        """
        this test just ensures that NeoFox does not crash with the synthetic data
        """
        data = [
            ("resources/synthetic_data/neoantigens_1patients_10neoantigens.2.txt",
             "resources/synthetic_data/patients_1patients_10neoantigens.2.txt"),
            ("resources/synthetic_data/neoantigens_10patients_10neoantigens.0.txt",
             "resources/synthetic_data/patients_10patients_10neoantigens.0.txt"),
            #("resources/synthetic_data/neoantigens_100patients_10neoantigens.2.txt",
            # "resources/synthetic_data/patients_100patients_10neoantigens.2.txt"),

            #("resources/synthetic_data/neoantigens_no_wt_1patients_10neoantigens.3.txt",
            # "resources/synthetic_data/patients_no_wt_1patients_10neoantigens.3.txt"),
            #("resources/synthetic_data/poltergeist_neoantigens.txt",
            # "resources/synthetic_data/poltergeist_patients.txt")
            ("resources/synthetic_data/neoantigens_no_wt_10patients_10neoantigens.4.txt",
             "resources/synthetic_data/patients_no_wt_10patients_10neoantigens.4.txt"),
            #("resources/synthetic_data/neoantigens_100patients_10neoantigens.4.txt",
            # "resources/synthetic_data/patients_100patients_10neoantigens.4.txt"),
        ]

        for n, p, in data:
            input_file = pkg_resources.resource_filename(
                neofox.tests.__name__, n)
            data = pd.read_csv(input_file, sep="\t")
            data = data.replace({np.nan: None})
            neoantigens = ModelConverter._neoantigens_csv2objects(data)
            patients_file = pkg_resources.resource_filename(
                neofox.tests.__name__, p)
            patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
            annotations = NeoFox(
                neoantigens=neoantigens,
                patients=patients,
                num_cpus=4,
            ).get_annotations()
            self.assertIsNotNone(annotations)

    def test_with_all_neoepitopes(self):
        """
        This test aims at testing neofox when MixMHCpred and MixMHC2pred are not configured. As these are optional it
        shoudl just run, but without these annotations in the output
        """
        annotations = NeoFox(
            neoantigens=self.neoantigens[2:3],  # only one as this is quite slow, the first one is an indel!
            patients=self.patients,
            num_cpus=4,
            with_all_neoepitopes=True
        ).get_annotations()

        found_recognition_potential = False
        self.assertIsNotNone(annotations)
        for n in annotations:
            self.assertIsNotNone(n.neoepitopes_mhc_i)
            self.assertIsNotNone(n.neoepitopes_mhc_i_i)
            for e in n.neoepitopes_mhc_i:
                self.assertIsNotNone(e.mutated_peptide)
                self.assertIsNotNone(e.wild_type_peptide)
                self.assertGreater(e.affinity_mutated, 0)
                self.assertGreater(e.affinity_wild_type, 0)
                self.assertGreater(e.rank_mutated, 0)
                self.assertLess(e.rank_mutated, neofox.RANK_MHCI_THRESHOLD_DEFAULT)
                self.assertGreater(e.rank_wild_type, 0)
                self.assertIsNotNone(e.allele_mhc_i)
                self.assertNotEqual(e.allele_mhc_i.name, '')
                self.assertIsNotNone(e.isoform_mhc_i_i)
                self.assertEqual(e.isoform_mhc_i_i.name, '')
                recognition_potential = EpitopeHelper.get_annotation_by_name(
                    e.neofox_annotations.annotations, name='recognition_potential')
                found_recognition_potential = found_recognition_potential or recognition_potential != 'NA'
            for e in n.neoepitopes_mhc_i_i:
                self.assertIsNotNone(e.mutated_peptide)
                self.assertIsNotNone(e.wild_type_peptide)
                self.assertGreater(e.affinity_mutated, 0)
                self.assertGreater(e.affinity_wild_type, 0)
                self.assertGreater(e.rank_mutated, 0)
                self.assertLess(e.rank_mutated, neofox.RANK_MHCII_THRESHOLD_DEFAULT)
                self.assertGreater(e.rank_wild_type, 0)
                self.assertIsNotNone(e.allele_mhc_i)
                self.assertEqual(e.allele_mhc_i.name, '')
                self.assertIsNotNone(e.isoform_mhc_i_i)
                self.assertNotEqual(e.isoform_mhc_i_i.name, '')

        self.assertTrue(found_recognition_potential)

        df_epitopes_mhci = ModelConverter.annotations2epitopes_table(annotations, mhc=neofox.MHC_I)
        self.assertFalse(any(c.startswith('isoformMhcII') for c in df_epitopes_mhci.columns))

        df_epitopes_mhcii = ModelConverter.annotations2epitopes_table(annotations, mhc=neofox.MHC_II)
        self.assertFalse(any(c.startswith('alleleMhcI') for c in df_epitopes_mhcii.columns))

    def _regression_test_on_output_file(self, new_file, previous_filename="resources/output_previous.txt"):
        previous_file = pkg_resources.resource_filename(neofox.tests.__name__, previous_filename)
        if os.path.exists(previous_file):
            NeofoxChecker(previous_file, new_file)
        else:
            logger.warning("No previous file to compare output with")
        # copies the new file to the previous file for the next execution only if no values were different
        shutil.copyfile(new_file, previous_file)


class NeofoxChecker:
    def __init__(self, previous_file, new_file):
        previous_df = pd.read_csv(previous_file, sep="\t")
        new_df = pd.read_csv(new_file, sep="\t")
        self._check_rows(new_df, previous_df)
        shared_columns = self._check_columns(new_df, previous_df)
        differing_columns = []
        for c in shared_columns:
            if self._check_values(c, new_df, previous_df):
                differing_columns.append(c)
        if len(differing_columns) > 0:
            differing_columns.sort()
            logger.error(
                "There are {} columns with differing values".format(
                    len(differing_columns)
                )
            )
            logger.error("Differing columns {}".format(differing_columns))
            logger.error("The regression test contains errors")
        else:
            logger.info("The regression test was successful!")

    def _check_values(self, column_name, new_df, previous_df):
        error = False
        ok_values_count = 0
        ko_values_count = 0
        for s1, s2 in zip(previous_df[column_name], new_df[column_name]):
            if self._check_single_value(s1, s2):
                ok_values_count += 1
            else:
                ko_values_count += 1

        if ok_values_count == 0:
            logger.error(
                "There are no equal values at all for column {}".format(column_name)
            )

        if ko_values_count > 0:
            logger.error(
                "There are {} different values for column {}".format(
                    ko_values_count, column_name
                )
            )
            logger.error(
                "Previous version: {}".format(
                    previous_df[column_name]
                    .transform(lambda x: x[0:20] + "..." if len(x) > 20 else x)
                    .values
                )
            )
            logger.error(
                "New version: {}".format(
                    new_df[column_name]
                    .transform(lambda x: x[0:20] + "..." if len(x) > 20 else x)
                    .values
                )
            )
            error = True

        return error

    def _check_single_value(self, s1, s2):
        is_float_s1 = self._is_float(s1)
        is_float_s2 = self._is_float(s2)
        if is_float_s1 and is_float_s2:
            # equality of NaN is never true so we force it
            # relative tolerance set to consider equal very close floats
            is_equal = (
                True
                if np.isnan(s1) and np.isnan(s2)
                else math.isclose(s1, s2, rel_tol=0.0001)
            )
        else:
            is_equal = s1 == s2
        return is_equal

    def _is_float(self, s1):
        return isinstance(s1, float) or isinstance(s1, np.float)

    def _check_columns(self, new_df, previous_df):
        shared_columns = set(previous_df.columns).intersection(set(new_df.columns))
        lost_columns = set(previous_df.columns).difference(set(new_df.columns))
        if len(lost_columns) > 0:
            logger.warning(
                "There are {} lost columns: {}".format(len(lost_columns), lost_columns)
            )
        gained_columns = set(new_df.columns).difference(set(previous_df.columns))
        if len(gained_columns) > 0:
            logger.warning(
                "There are {} gained columns: {}".format(
                    len(gained_columns), gained_columns
                )
            )
        # fails the test if there are no shared columns
        assert len(shared_columns) > 0, "No shared columns"
        return shared_columns

    def _check_rows(self, new_df, previous_df):
        # fails the test if the number of rows differ
        assert previous_df.shape[0] == new_df.shape[0], "Mismatching number of rows"


