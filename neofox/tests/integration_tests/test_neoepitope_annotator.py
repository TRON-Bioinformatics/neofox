from unittest import TestCase

from neofox.annotator.neoepitope_annotator import NeoepitopeAnnotator
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.factories import MhcFactory
from neofox.model.neoantigen import PredictedEpitope, MhcAllele, Mhc2Isoform, Patient
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.tests.integration_tests import integration_test_tools
from neofox.tests.integration_tests.integration_test_tools import get_hla_one_test, get_hla_two_test


class NeoepitopeAnnotatorTest(TestCase):

    def setUp(self) -> None:
        self.references, self.configuration = integration_test_tools.load_references()
        self.annotator = NeoepitopeAnnotator(
            references=self.references,
            configuration=self.configuration,
            tcell_predictor=TcellPrediction(),
            self_similarity=SelfSimilarityCalculator()
        )
        self.patient = Patient(
            identifier="123",
            mhc1=get_hla_one_test(self.references.get_mhc_database()),
            mhc2=get_hla_two_test(self.references.get_mhc_database())
        )

    def test_neoepitope_mhci_annotation(self):

        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            wild_type_peptide="DILVIDQTR",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01')
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope)
        self._assert_neoepitope_mhci(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)

    def test_neoepitope_mhci_without_wild_type(self):

        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01')
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope)
        self._assert_neoepitope_mhci(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)

    def test_neoepitope_mhci_9mer_with_frequencies_and_gene(self):
        """
        this checks fields that are only annotated when expression, vaf and/or gene are provided
        """
        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            wild_type_peptide="DILVIDQTR",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            rna_variant_allele_frequency=0.5,
            dna_variant_allele_frequency=1.0,
            rna_expression=125,
            gene='BRCA2'
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope)
        self._assert_neoepitope_mhci(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)
        self._assert_float_annotation(annotated_neoepitope, annotation_name="Priority_score")
        self._assert_float_annotation(annotated_neoepitope, annotation_name="Tcell_predictor_score")

    def test_neoepitope_mhci_10mer_no_tcell_predictor(self):
        """
        this checks fields that are only annotated when expression, vaf and/or gene are provided
        """
        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTRA",
            wild_type_peptide="DILVIDQTRA",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            gene='BRCA2'
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope)
        self._assert_neoepitope_mhci(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)
        annotation_value = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, "Tcell_predictor_score")
        self.assertEqual(annotation_value, "NA")

    def test_neoepitope_mhci_without_dna_vaf(self):
        """
        this checks fields that are only annotated when expression, vaf and/or gene are provided
        """
        neoepitope_with_dna_vaf = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            wild_type_peptide="DILVIDQTR",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            dna_variant_allele_frequency=1.0,
            rna_variant_allele_frequency=0.1,
            rna_expression=125,
            gene='BRCA2'
        )
        neoepitope_without_dna_vaf = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            wild_type_peptide="DILVIDQTR",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            rna_variant_allele_frequency=0.1,
            rna_expression=125,
            gene='BRCA2'
        )

        annotated_neoepitope1 = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope_with_dna_vaf)
        self._assert_neoepitope_mhci(original_neoepitope=neoepitope_with_dna_vaf,
                                     annotated_neoepitope=annotated_neoepitope1)

        annotated_neoepitope2 = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope_without_dna_vaf)
        self._assert_neoepitope_mhci(original_neoepitope=neoepitope_without_dna_vaf,
                                     annotated_neoepitope=annotated_neoepitope2)

        self.assertNotEqual(
            EpitopeHelper.get_annotation_by_name(
                annotated_neoepitope1.neofox_annotations.annotations, "Priority_score"),
            EpitopeHelper.get_annotation_by_name(
                annotated_neoepitope2.neofox_annotations.annotations, "Priority_score")
        )

    def test_neoepitope_mhci_without_vaf(self):
        """
        this checks fields that are only annotated when expression, vaf and/or gene are provided
        """
        neoepitope_without_vaf = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            wild_type_peptide="DILVIDQTR",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            dna_variant_allele_frequency=None,
            rna_variant_allele_frequency=None,
            gene='BRCA2'
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope_without_vaf)
        self._assert_neoepitope_mhci(original_neoepitope=neoepitope_without_vaf,
                                     annotated_neoepitope=annotated_neoepitope)

        self.assertEqual(
            EpitopeHelper.get_annotation_by_name(
                annotated_neoepitope.neofox_annotations.annotations, "Priority_score"), "NA")

    def test_neoepitope_mhcii_annotation(self):

        isoform = MhcFactory.build_mhc2_alleles(["HLA-DRB1*01:01"], self.references.get_mhc_database())[0].isoforms[0]
        neoepitope = PredictedEpitope(
            mutated_peptide="DEVLGEPSQDILVTDQTR",
            wild_type_peptide="DEVLGEPSQDILVIDQTR",
            isoform_mhc_i_i=isoform
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope)
        self._assert_neoepitope_mhcii(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)

    def test_neoepitope_mhcii_without_wild_type(self):

        isoform = MhcFactory.build_mhc2_alleles(["HLA-DRB1*01:01"], self.references.get_mhc_database())[0].isoforms[0]
        neoepitope = PredictedEpitope(
            mutated_peptide="DEVLGEPSQDILVTDQTR",
            isoform_mhc_i_i=isoform
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope)
        self._assert_neoepitope_mhcii(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)

    def _get_test_mhci_allele(self, allele) -> MhcAllele:
        mhci = MhcFactory.build_mhc1_alleles([allele], mhc_database=self.references.get_mhc_database())
        return mhci[0].alleles[0]

    def _get_test_mhcii_isoform(self, isoform) -> Mhc2Isoform:
        mhcii = MhcFactory.build_mhc2_alleles([isoform], mhc_database=self.references.get_mhc_database())
        return mhcii[0].isoforms[0]

    def _assert_neoepitope_mhci(self, original_neoepitope: PredictedEpitope, annotated_neoepitope: PredictedEpitope):
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
        self._assert_float_annotation(annotated_neoepitope, annotation_name="MixMHCpred_affinity_score")
        self._assert_float_annotation(annotated_neoepitope, annotation_name="MixMHCpred_rank")
        self._assert_float_annotation(annotated_neoepitope, annotation_name="MixMHCpred_WT_affinity_score")
        self._assert_float_annotation(annotated_neoepitope, annotation_name="MixMHCpred_WT_rank")

        # PRIME annotations
        self._assert_float_annotation(annotated_neoepitope, annotation_name="PRIME_affinity_score")
        self._assert_float_annotation(annotated_neoepitope, annotation_name="PRIME_rank")
        self._assert_float_annotation(annotated_neoepitope, annotation_name="PRIME_WT_affinity_score")
        self._assert_float_annotation(annotated_neoepitope, annotation_name="PRIME_WT_rank")

        # additional annotations
        self._assert_annotation(annotated_neoepitope, annotation_name="position_mutation")
        self._assert_annotation(annotated_neoepitope, annotation_name="anchor_mutated")
        self._assert_annotation(annotated_neoepitope, annotation_name="amplitude")
        self._assert_annotation(annotated_neoepitope, annotation_name="pathogen_similarity")
        self._assert_annotation(annotated_neoepitope, annotation_name="recognition_potential")
        self._assert_annotation(annotated_neoepitope, annotation_name="DAI")
        self._assert_annotation(annotated_neoepitope, annotation_name="Improved_Binder_MHCI")
        self._assert_annotation(annotated_neoepitope, annotation_name="Selfsimilarity")
        self._assert_annotation(annotated_neoepitope, annotation_name="Selfsimilarity_conserved_binder")
        self._assert_annotation(annotated_neoepitope, annotation_name="mutation_not_found_in_proteome")
        self._assert_annotation(annotated_neoepitope, annotation_name="dissimilarity_score")
        self._assert_annotation(annotated_neoepitope, annotation_name="number_of_mismatches")
        self._assert_annotation(annotated_neoepitope, annotation_name="IEDB_Immunogenicity")
        self._assert_annotation(annotated_neoepitope, annotation_name="hex_alignment_score")

        # others to comes
        self._assert_annotation(annotated_neoepitope, annotation_name="Priority_score")
        self._assert_annotation(annotated_neoepitope, annotation_name="Tcell_predictor_score")

    def _assert_neoepitope_mhcii(self, original_neoepitope: PredictedEpitope, annotated_neoepitope: PredictedEpitope):
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
        self._assert_float_annotation(annotated_neoepitope, annotation_name="MixMHC2pred_rank")
        self._assert_float_annotation(annotated_neoepitope, annotation_name="MixMHC2pred_WT_rank")

        # additional annotations
        self._assert_annotation(annotated_neoepitope, annotation_name="amplitude")
        self._assert_annotation(annotated_neoepitope, annotation_name="pathogen_similarity")
        self._assert_annotation(annotated_neoepitope, annotation_name="Selfsimilarity")
        self._assert_annotation(annotated_neoepitope, annotation_name="mutation_not_found_in_proteome")
        self._assert_annotation(annotated_neoepitope, annotation_name="dissimilarity_score")
        self._assert_annotation(annotated_neoepitope, annotation_name="IEDB_Immunogenicity")
        self._assert_annotation(annotated_neoepitope, annotation_name="hex_alignment_score")

    def _assert_float_annotation(self, annotated_neoepitope, annotation_name):
        annotation_value = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, annotation_name)
        self.assertIsInstance(annotation_value, str)
        self.assertIsInstance(float(annotation_value), float)

    def _assert_annotation(self, annotated_neoepitope, annotation_name):
        annotation_value = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, annotation_name)
        self.assertIsInstance(annotation_value, str)
