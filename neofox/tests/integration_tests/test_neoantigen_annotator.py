from unittest import TestCase

from neofox.annotator import NeoantigenAnnotator
from neofox.model.factories import MhcFactory
from neofox.model.neoantigen import PredictedEpitope, MhcAllele, Neoantigen, Mhc2Isoform
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.tests.integration_tests import integration_test_tools


class NeoantigenAnnotatorTest(TestCase):

    def setUp(self) -> None:
        self.references, self.configuration = integration_test_tools.load_references()
        self.annotator = NeoantigenAnnotator(
            references=self.references,
            configuration=self.configuration,
            tcell_predictor=TcellPrediction(),
            self_similarity=SelfSimilarityCalculator()
        )

    def test_neoepitope_annotation_mhci(self):
        epitope = PredictedEpitope(
            mutated_peptide="AAAAAADAAAAA",
            wild_type_peptide="AAAAAAAAAAAA",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01')
        )
        neoantigen = Neoantigen(gene='BRCA2', dna_variant_allele_frequency=0.5, rna_variant_allele_frequency=0.5,
                                rna_expression=0.5)
        annotated_epitope = self.annotator.get_additional_annotations_neoepitope_mhci(
            epitope=epitope,
            neoantigen=neoantigen,
            vaf_rna=0.5
        )
        self.assertIsNotNone(annotated_epitope)
        self.assertEqual(annotated_epitope.mutated_peptide, epitope.mutated_peptide)
        self.assertEqual(annotated_epitope.wild_type_peptide, epitope.wild_type_peptide)
        self.assertEqual(annotated_epitope.allele_mhc_i.name, epitope.allele_mhc_i.name)
        self.assertGreater(len(annotated_epitope.neofox_annotations.annotations), 0)

    def test_neoepitope_annotation_mhci_without_wild_type(self):
        epitope = PredictedEpitope(
            mutated_peptide="AAAAAADAAAAA",
            #wild_type_peptide="AAAAAAAAAAAA",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01')
        )
        neoantigen = Neoantigen(gene='BRCA2', dna_variant_allele_frequency=0.5, rna_variant_allele_frequency=0.5,
                                rna_expression=0.5)
        annotated_epitope = self.annotator.get_additional_annotations_neoepitope_mhci(
            epitope=epitope,
            neoantigen=neoantigen,
            vaf_rna=0.5
        )
        self.assertIsNotNone(annotated_epitope)
        self.assertEqual(annotated_epitope.mutated_peptide, epitope.mutated_peptide)
        self.assertEqual(annotated_epitope.wild_type_peptide, epitope.wild_type_peptide)
        self.assertEqual(annotated_epitope.allele_mhc_i.name, epitope.allele_mhc_i.name)
        self.assertGreater(len(annotated_epitope.neofox_annotations.annotations), 0)

    def test_neoepitope_annotation_mhcii(self):
        epitope = PredictedEpitope(
            mutated_peptide="AAAAAADAAAAA",
            wild_type_peptide="AAAAAAAAAAAA",
            allele_mhc_i=MhcAllele(name='HLA-A*01:01'),
            isoform_mhc_i_i=self._get_test_mhcii_isoform('HLA-DRB1*04:02')
        )
        annotated_epitope = self.annotator.get_additional_annotations_neoepitope_mhcii(epitope=epitope)
        self.assertIsNotNone(annotated_epitope)
        self.assertEqual(annotated_epitope.mutated_peptide, epitope.mutated_peptide)
        self.assertEqual(annotated_epitope.wild_type_peptide, epitope.wild_type_peptide)
        self.assertEqual(annotated_epitope.allele_mhc_i.name, epitope.allele_mhc_i.name)
        self.assertGreater(len(annotated_epitope.neofox_annotations.annotations), 0)

    def _get_test_mhci_allele(self, allele) -> MhcAllele:
        mhci = MhcFactory.build_mhc1_alleles([allele], mhc_database=self.references.get_mhc_database())
        return mhci[0].alleles[0]

    def _get_test_mhcii_isoform(self, isoform) -> Mhc2Isoform:
        mhcii = MhcFactory.build_mhc2_alleles([isoform], mhc_database=self.references.get_mhc_database())
        return mhcii[0].isoforms[0]
