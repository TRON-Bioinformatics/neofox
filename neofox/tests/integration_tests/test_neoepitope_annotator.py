from neofox.annotator.neoepitope_annotator import NeoepitopeAnnotator
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.factories import MhcFactory
from neofox.model.neoantigen import PredictedEpitope, MhcAllele, Mhc2Isoform, Patient
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.tests.integration_tests.integration_test_tools import get_hla_one_test, get_hla_two_test, \
    BaseIntegrationTest


class NeoepitopeAnnotatorTest(BaseIntegrationTest):

    def setUp(self) -> None:
        super().setUp()
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
        self.assert_neoepitope_mhci(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)

    def test_neoepitope_mhci_without_wild_type(self):

        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01')
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope)
        self.assert_neoepitope_mhci(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)

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
        self.assert_neoepitope_mhci(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)
        self.assert_float_annotation(annotated_neoepitope, annotation_name="Priority_score")
        self.assert_float_annotation(annotated_neoepitope, annotation_name="Tcell_predictor")

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
        self.assert_neoepitope_mhci(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)
        annotation_value = EpitopeHelper.get_annotation_by_name(
            annotated_neoepitope.neofox_annotations.annotations, "Tcell_predictor")
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
        self.assert_neoepitope_mhci(original_neoepitope=neoepitope_with_dna_vaf,
                                     annotated_neoepitope=annotated_neoepitope1)

        annotated_neoepitope2 = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope_without_dna_vaf)
        self.assert_neoepitope_mhci(original_neoepitope=neoepitope_without_dna_vaf,
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
        self.assert_neoepitope_mhci(original_neoepitope=neoepitope_without_vaf,
                                     annotated_neoepitope=annotated_neoepitope)

        self.assertEqual(
            EpitopeHelper.get_annotation_by_name(
                annotated_neoepitope.neofox_annotations.annotations, "Priority_score"), "NA")

    def test_neoepitope_mhcii_annotation(self):

        neoepitope = PredictedEpitope(
            mutated_peptide="DEVLGEPSQDILVTDQTR",
            wild_type_peptide="DEVLGEPSQDILVIDQTR",
            isoform_mhc_i_i=self._get_test_mhcii_isoform("HLA-DRB1*01:01")
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope)
        self.assert_neoepitope_mhcii(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)

    def test_neoepitope_mhcii_without_wild_type(self):

        neoepitope = PredictedEpitope(
            mutated_peptide="DEVLGEPSQDILVTDQTR",
            isoform_mhc_i_i=self._get_test_mhcii_isoform("HLA-DRB1*01:01")
        )

        annotated_neoepitope = self.annotator.get_annotated_neoepitope(neoepitope=neoepitope)
        self.assert_neoepitope_mhcii(original_neoepitope=neoepitope, annotated_neoepitope=annotated_neoepitope)
