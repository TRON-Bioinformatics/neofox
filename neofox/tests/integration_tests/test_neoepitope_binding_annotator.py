from unittest import TestCase

from neofox.annotation_resources.uniprot.uniprot import Uniprot
from neofox.annotator.annotator import NeoantigenAnnotator
from neofox.annotator.neoepitope_mhc_binding_annotator import NeoepitopeMhcBindingAnnotator
from neofox.helpers.blastp_runner import BlastpRunner
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.helpers.runner import Runner
from neofox.model.factories import MhcFactory, NeoantigenFactory
from neofox.model.neoantigen import PredictedEpitope, MhcAllele, Neoantigen, Mhc2Isoform, Patient
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.tests.integration_tests import integration_test_tools
from neofox.tests.integration_tests.integration_test_tools import get_hla_one_test, get_hla_two_test


class NeoepitopeMhcBindingAnnotatorTest(TestCase):

    def setUp(self) -> None:
        self.references, self.configuration = integration_test_tools.load_references()
        proteome_blastp_runner = BlastpRunner(
            runner=Runner(), configuration=self.configuration,
            database=self.references.get_proteome_database())
        self.annotator = NeoepitopeMhcBindingAnnotator(
            references=self.references,
            configuration=self.configuration,
            proteome_blastp_runner=proteome_blastp_runner,
            uniprot=Uniprot(self.references.uniprot_pickle))

    def test_neoepitope_mhc1(self):

        allele = MhcFactory.build_mhc1_alleles(["HLA-A*03:01"], self.references.get_mhc_database())[0].alleles[0]
        neoepitope = PredictedEpitope(
            mutated_peptide="PSQDILVID",
            wild_type_peptide="PSQDILVTD",
            allele_mhc_i=allele
        )
        annotated_neoepitope = self.annotator.get_mhc_binding_annotations(neoepitope=neoepitope)

        self.assertIsNotNone(annotated_neoepitope)
        # input data remains the same
        self.assertEqual(annotated_neoepitope.mutated_peptide, neoepitope.mutated_peptide)
        self.assertEqual(annotated_neoepitope.wild_type_peptide, neoepitope.wild_type_peptide)
        self.assertEqual(annotated_neoepitope.allele_mhc_i.name, neoepitope.allele_mhc_i.name)
        # netMHCpan annotations
        self.assertIsInstance(annotated_neoepitope.rank_mutated, float)
        self.assertIsInstance(annotated_neoepitope.rank_wild_type, float)
        self.assertIsInstance(annotated_neoepitope.affinity_mutated, float)
        self.assertIsInstance(annotated_neoepitope.affinity_wild_type, float)
        # MixMHCpred annotations
        mixmhcpred_affinity = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, "MixMHCpred_affinity_score")
        self.assertIsInstance(mixmhcpred_affinity, str)
        self.assertIsInstance(float(mixmhcpred_affinity), float)
        mixmhcpred_rank = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, "MixMHCpred_rank")
        self.assertIsInstance(mixmhcpred_rank, str)
        self.assertIsInstance(float(mixmhcpred_rank), float)
        # PRIME annotations
        prime_affinity = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, "PRIME_affinity_score")
        self.assertIsInstance(prime_affinity, str)
        self.assertIsInstance(float(prime_affinity), float)
        prime_rank = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, "PRIME_rank")
        self.assertIsInstance(prime_rank, str)
        self.assertIsInstance(float(prime_rank), float)

    def test_neoepitope_mhc2(self):

        isoform = MhcFactory.build_mhc2_alleles(["HLA-DRB1*01:01"], self.references.get_mhc_database())[0].isoforms[0]
        neoepitope = PredictedEpitope(
            mutated_peptide="DEVLGEPSQDILVTDQTR",
            wild_type_peptide="DEVLGEPSQDILVIDQTR",
            isoform_mhc_i_i=isoform
        )
        annotated_neoepitope = self.annotator.get_mhc_binding_annotations(neoepitope=neoepitope)

        self.assertIsNotNone(annotated_neoepitope)
        # input data remains the same
        self.assertEqual(annotated_neoepitope.mutated_peptide, neoepitope.mutated_peptide)
        self.assertEqual(annotated_neoepitope.wild_type_peptide, neoepitope.wild_type_peptide)
        self.assertEqual(annotated_neoepitope.isoform_mhc_i_i.name, neoepitope.isoform_mhc_i_i.name)
        # netMHC2pan annotations
        self.assertIsInstance(annotated_neoepitope.rank_mutated, float)
        self.assertIsInstance(annotated_neoepitope.rank_wild_type, float)
        self.assertIsInstance(annotated_neoepitope.affinity_mutated, float)
        self.assertIsInstance(annotated_neoepitope.affinity_wild_type, float)
        # MixMHC2pred annotations
        mixmhcpred_rank = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, "MixMHC2pred_rank")
        self.assertIsInstance(mixmhcpred_rank, str)
        self.assertIsInstance(float(mixmhcpred_rank), float)
