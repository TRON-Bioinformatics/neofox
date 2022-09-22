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
import random
import tempfile
from unittest import TestCase

import dotenv
from Bio.Alphabet.IUPAC import IUPACData

from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.factories import MhcFactory
from neofox.model.neoantigen import PredictedEpitope, MhcAllele, Mhc2Isoform
from neofox.references.references import ReferenceFolder, DependenciesConfiguration, ORGANISM_HOMO_SAPIENS, \
    ORGANISM_MUS_MUSCULUS


def load_references(organism=ORGANISM_HOMO_SAPIENS):
    dotenv.load_dotenv(override=True)
    return ReferenceFolder(organism=organism), DependenciesConfiguration()


def create_temp_aminoacid_fasta_file():
    fastafile = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with fastafile as f:
        f.write(get_random_kmer())
    return fastafile


def get_random_kmer(k=25):
    return "".join(random.choices(list(IUPACData.protein_letters), k=k))


def get_hla_one_test(hla_database):
    return MhcFactory.build_mhc1_alleles(
        [
            "HLA-A*24:02",
            "HLA-A*02:01",
            "HLA-B*15:01",
            "HLA-B*44:02",
            "HLA-C*07:02",
            "HLA-C*05:01",
        ], hla_database
    )


def get_h2_one_test(h2_database):
    return MhcFactory.build_mhc1_alleles(
        [
            "H2Kd",
            "H2Kd",
            "H2Dd",
            "H2Dd",
            "H2Ld",
            "H2Ld",
        ], h2_database
    )


def get_hla_two_test(hla_database):
    return MhcFactory.build_mhc2_alleles(
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
                "HLA-DPB1*04:01",
            ], hla_database
        )


def get_h2_two_test(h2_database):
    return MhcFactory.build_mhc2_alleles(
        [
            "H2Ad",
            "H2Ad",
            "H2Ed",
            "H2Ed"
        ], h2_database
    )


mutations_with_rare_aminoacids = [
            ("UTTDSDGKF", "UTTDSWGKF"),  # this is an epitope from IEDB of length 9
            ("XTTDSDGKF", "XTTDSWGKF"),
            ("BTTDSDGKF", "BTTDSWGKF"),
            ("JTTDSDGKF", "JTTDSWGKF"),
            ("OTTDSDGKF", "OTTDSWGKF"),
            ("ZTTDSDGKF", "ZTTDSWGKF"),
            # only present in the wild type
            ("UTTDSDGKF", "TTTDSWGKF"),  # this is an epitope from IEDB of length 9
            ("XTTDSDGKF", "TTTDSWGKF"),
            ("BTTDSDGKF", "TTTDSWGKF"),
            ("JTTDSDGKF", "TTTDSWGKF"),
            ("OTTDSDGKF", "TTTDSWGKF"),
            ("ZTTDSDGKF", "TTTDSWGKF"),
            # only present in the mutation
            ("TTTDSDGKF", "UTTDSWGKF"),  # this is an epitope from IEDB of length 9
            ("TTTDSDGKF", "XTTDSWGKF"),
            ("TTTDSDGKF", "BTTDSWGKF"),
            ("TTTDSDGKF", "JTTDSWGKF"),
            ("TTTDSDGKF", "OTTDSWGKF"),
            ("TTTDSDGKF", "ZTTDSWGKF")
        ]


class BaseIntegrationTest(TestCase):

    def setUp(self):
        self.references, self.configuration = load_references()
        self.references_mouse, self.configuration_mouse = load_references(
            organism=ORGANISM_MUS_MUSCULUS)

    def _get_test_mhci_allele(self, allele) -> MhcAllele:
        mhci = MhcFactory.build_mhc1_alleles([allele], mhc_database=self.references.get_mhc_database())
        return mhci[0].alleles[0]

    def _get_test_mhcii_isoform(self, isoform) -> Mhc2Isoform:
        mhcii = MhcFactory.build_mhc2_alleles([isoform], mhc_database=self.references.get_mhc_database())
        return mhcii[0].isoforms[0]

    def assert_float_annotation(self, annotated_neoepitope, annotation_name):
        annotation_value = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, annotation_name)
        self.assertIsInstance(annotation_value, str)
        self.assertIsInstance(float(annotation_value), float)

    def assert_annotation(self, annotated_neoepitope, annotation_name):
        annotation_value = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, annotation_name)
        self.assertIsInstance(annotation_value, str)

    def assert_neoepitope_mhci(self, original_neoepitope: PredictedEpitope, annotated_neoepitope: PredictedEpitope):
        self.assertIsNotNone(annotated_neoepitope)
        # input data remains the same
        self.assertEqual(annotated_neoepitope.mutated_peptide, original_neoepitope.mutated_peptide)
        if original_neoepitope.wild_type_peptide is not None and original_neoepitope.wild_type_peptide != '':
            self.assertEqual(annotated_neoepitope.wild_type_peptide, original_neoepitope.wild_type_peptide)
        else:
            self.assertIsNotNone(annotated_neoepitope.wild_type_peptide)
            self.assertNotEqual(annotated_neoepitope.wild_type_peptide, '')
        self.assertEqual(annotated_neoepitope.allele_mhc_i.name, original_neoepitope.allele_mhc_i.name)

        # netMHCpan annotations
        self.assertIsInstance(annotated_neoepitope.rank_mutated, float)
        self.assertIsInstance(annotated_neoepitope.rank_wild_type, float)
        self.assertIsInstance(annotated_neoepitope.affinity_mutated, float)
        self.assertIsInstance(annotated_neoepitope.affinity_wild_type, float)

        # MixMHCpred annotations
        self.assert_float_annotation(annotated_neoepitope, annotation_name="MixMHCpred_score")
        self.assert_float_annotation(annotated_neoepitope, annotation_name="MixMHCpred_rank")
        self.assert_float_annotation(annotated_neoepitope, annotation_name="MixMHCpred_WT_score")
        self.assert_float_annotation(annotated_neoepitope, annotation_name="MixMHCpred_WT_rank")

        # PRIME annotations
        self.assert_float_annotation(annotated_neoepitope, annotation_name="PRIME_score")
        self.assert_float_annotation(annotated_neoepitope, annotation_name="PRIME_rank")
        self.assert_float_annotation(annotated_neoepitope, annotation_name="PRIME_WT_score")
        self.assert_float_annotation(annotated_neoepitope, annotation_name="PRIME_WT_rank")

        # additional annotations
        self.assert_annotation(annotated_neoepitope, annotation_name="position_mutation")
        self.assert_annotation(annotated_neoepitope, annotation_name="anchor_mutated")
        self.assert_annotation(annotated_neoepitope, annotation_name="amplitude")
        self.assert_annotation(annotated_neoepitope, annotation_name="pathogen_similarity")
        self.assert_annotation(annotated_neoepitope, annotation_name="recognition_potential")
        self.assert_annotation(annotated_neoepitope, annotation_name="DAI")
        self.assert_annotation(annotated_neoepitope, annotation_name="Improved_Binder_MHCI")
        self.assert_annotation(annotated_neoepitope, annotation_name="Selfsimilarity")
        self.assert_annotation(annotated_neoepitope, annotation_name="Selfsimilarity_conserved_binder")
        self.assert_annotation(annotated_neoepitope, annotation_name="mutation_not_found_in_proteome")
        self.assert_annotation(annotated_neoepitope, annotation_name="dissimilarity_score")
        self.assert_annotation(annotated_neoepitope, annotation_name="number_of_mismatches")
        self.assert_annotation(annotated_neoepitope, annotation_name="IEDB_Immunogenicity")
        self.assert_annotation(annotated_neoepitope, annotation_name="hex_alignment_score")

        # others to comes
        self.assert_annotation(annotated_neoepitope, annotation_name="Priority_score")
        self.assert_annotation(annotated_neoepitope, annotation_name="Tcell_predictor")

    def assert_neoepitope_mhcii(self, original_neoepitope: PredictedEpitope, annotated_neoepitope: PredictedEpitope):
        self.assertIsNotNone(annotated_neoepitope)
        # input data remains the same
        self.assertEqual(annotated_neoepitope.mutated_peptide, original_neoepitope.mutated_peptide)
        if original_neoepitope.wild_type_peptide is not None and original_neoepitope.wild_type_peptide != '':
            self.assertEqual(annotated_neoepitope.wild_type_peptide, original_neoepitope.wild_type_peptide)
        else:
            self.assertIsNotNone(annotated_neoepitope.wild_type_peptide)
            self.assertNotEqual(annotated_neoepitope.wild_type_peptide, '')
        self.assertEqual(annotated_neoepitope.isoform_mhc_i_i.name, original_neoepitope.isoform_mhc_i_i.name)

        # netMHCpan annotations
        self.assertIsInstance(annotated_neoepitope.rank_mutated, float)
        self.assertIsInstance(annotated_neoepitope.rank_wild_type, float)
        self.assertIsInstance(annotated_neoepitope.affinity_mutated, float)
        self.assertIsInstance(annotated_neoepitope.affinity_wild_type, float)

        # MixMHCpred annotations
        self.assert_float_annotation(annotated_neoepitope, annotation_name="MixMHC2pred_rank")
        self.assert_float_annotation(annotated_neoepitope, annotation_name="MixMHC2pred_WT_rank")

        # additional annotations
        self.assert_annotation(annotated_neoepitope, annotation_name="amplitude")
        self.assert_annotation(annotated_neoepitope, annotation_name="pathogen_similarity")
        self.assert_annotation(annotated_neoepitope, annotation_name="Selfsimilarity")
        self.assert_annotation(annotated_neoepitope, annotation_name="mutation_not_found_in_proteome")
        self.assert_annotation(annotated_neoepitope, annotation_name="dissimilarity_score")
        self.assert_annotation(annotated_neoepitope, annotation_name="IEDB_Immunogenicity")
        self.assert_annotation(annotated_neoepitope, annotation_name="hex_alignment_score")
