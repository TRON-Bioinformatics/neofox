from unittest import TestCase

from neofox.annotator.neoantigen_annotator import NeoantigenAnnotator
from neofox.model.factories import MhcFactory, NeoantigenFactory
from neofox.model.neoantigen import PredictedEpitope, MhcAllele, Neoantigen, Mhc2Isoform, Patient
from neofox.published_features.Tcell_predictor.tcellpredictor_wrapper import TcellPrediction
from neofox.published_features.self_similarity.self_similarity import SelfSimilarityCalculator
from neofox.tests.integration_tests import integration_test_tools
from neofox.tests.integration_tests.integration_test_tools import get_hla_one_test, get_hla_two_test


class NeoantigenAnnotatorTest(TestCase):

    def setUp(self) -> None:
        self.references, self.configuration = integration_test_tools.load_references()
        self.annotator = NeoantigenAnnotator(
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

    def test_neoantigen_annotation(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
            mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
            patient_identifier="123",
            gene="BRCA2"
        )
        annotated_neoantigen = self.annotator.get_annotated_neoantigen(neoantigen=neoantigen, patient=self.patient)
        self._assert_neoantigen(annotated_neoantigen, neoantigen)
        self._assert_epitopes(annotated_neoantigen=annotated_neoantigen)

    def test_neoantigen_annotation_with_all_epitopes(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
            mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
            patient_identifier="123",
            gene="BRCA2"
        )
        annotated_neoantigen = self.annotator.get_annotated_neoantigen(
            neoantigen=neoantigen, patient=self.patient, with_all_neoepitopes=True)
        self._assert_neoantigen(annotated_neoantigen, neoantigen)
        self._assert_epitopes(annotated_neoantigen=annotated_neoantigen, with_all_epitopes=True)

    def test_neoantigen_annotation_without_wild_type(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
            patient_identifier="123",
            gene="BRCA2"
        )
        annotated_neoantigen = self.annotator.get_annotated_neoantigen(neoantigen=neoantigen, patient=self.patient)
        self._assert_neoantigen(annotated_neoantigen, neoantigen)
        # wild type xmer is still empty!
        self.assertIsNone(annotated_neoantigen.wild_type_xmer)
        self._assert_epitopes(annotated_neoantigen=annotated_neoantigen)

    def test_neoantigen_annotation_without_wild_type_and_with_all_epitopes(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
            patient_identifier="123",
            gene="BRCA2"
        )
        annotated_neoantigen = self.annotator.get_annotated_neoantigen(
            neoantigen=neoantigen, patient=self.patient, with_all_neoepitopes=True)
        self._assert_neoantigen(annotated_neoantigen, neoantigen)
        # wild type xmer is still empty!
        self.assertIsNone(annotated_neoantigen.wild_type_xmer)
        self._assert_epitopes(annotated_neoantigen, with_all_epitopes=True)

    def test_neoantigen_annotation_with_vaf_and_without_tx_expression(self):
        neoantigen = NeoantigenFactory.build_neoantigen(
            wild_type_xmer="DEVLGEPSQDILVIDQTRLEATISPET",
            mutated_xmer="DEVLGEPSQDILVTDQTRLEATISPET",
            patient_identifier="123",
            gene="BRCA2",
            dna_variant_allele_frequency=0.5
        )
        annotated_neoantigen = self.annotator.get_annotated_neoantigen(
            neoantigen=neoantigen, patient=self.patient, with_all_neoepitopes=True)
        self._assert_neoantigen(annotated_neoantigen, neoantigen)

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
            neoantigen=neoantigen
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
            neoantigen=neoantigen
        )
        self.assertIsNotNone(annotated_epitope)
        self.assertEqual(annotated_epitope.mutated_peptide, epitope.mutated_peptide)
        self.assertEqual(annotated_epitope.wild_type_peptide, epitope.wild_type_peptide)
        self.assertEqual(annotated_epitope.allele_mhc_i.name, epitope.allele_mhc_i.name)
        self.assertGreater(len(annotated_epitope.neofox_annotations.annotations), 0)

    def test_neoepitope_annotation_mhcii(self):
        epitope = PredictedEpitope(
            mutated_peptide="AAAAAADAAAAAAA",
            wild_type_peptide="AAAAAAAAAAAAAA",
            isoform_mhc_i_i=self._get_test_mhcii_isoform('HLA-DRB1*04:02')
        )
        annotated_epitope = self.annotator.get_additional_annotations_neoepitope_mhcii(epitope=epitope)
        self.assertIsNotNone(annotated_epitope)
        self.assertEqual(annotated_epitope.mutated_peptide, epitope.mutated_peptide)
        self.assertEqual(annotated_epitope.wild_type_peptide, epitope.wild_type_peptide)
        self.assertEqual(annotated_epitope.isoform_mhc_i_i.name, epitope.isoform_mhc_i_i.name)
        self.assertGreater(len(annotated_epitope.neofox_annotations.annotations), 0)

    def _get_test_mhci_allele(self, allele) -> MhcAllele:
        mhci = MhcFactory.build_mhc1_alleles([allele], mhc_database=self.references.get_mhc_database())
        return mhci[0].alleles[0]

    def _get_test_mhcii_isoform(self, isoform) -> Mhc2Isoform:
        mhcii = MhcFactory.build_mhc2_alleles([isoform], mhc_database=self.references.get_mhc_database())
        return mhcii[0].isoforms[0]

    def _assert_neoantigen(self, annotated_neoantigen: Neoantigen, neoantigen: Neoantigen):
        self.assertIsNotNone(annotated_neoantigen)
        self.assertEqual(annotated_neoantigen.mutated_xmer, neoantigen.mutated_xmer)
        self.assertEqual(annotated_neoantigen.wild_type_xmer, neoantigen.wild_type_xmer)
        self.assertEqual(annotated_neoantigen.position, neoantigen.position)
        self.assertGreater(len(annotated_neoantigen.neofox_annotations.annotations), 0)
        annotation_names = [a.name for a in annotated_neoantigen.neofox_annotations.annotations]
        self.assertTrue("NetMHCpan_bestRank_peptide" in annotation_names)
        self.assertTrue("NetMHCpan_bestRank_allele" in annotation_names)

    def _assert_epitopes(self, annotated_neoantigen, with_all_epitopes=False):
        # neoepitopes for both MHC I and MHC II are not empty
        self.assertGreater(len(annotated_neoantigen.neoepitopes_mhc_i), 0)
        self.assertGreater(len(annotated_neoantigen.neoepitopes_mhc_i_i), 0)
        observed_mixmhcpred_annotations = False
        for e in annotated_neoantigen.neoepitopes_mhc_i + annotated_neoantigen.neoepitopes_mhc_i_i:
            # WT peptides are not empty!
            self.assertIsNotNone(e.wild_type_peptide)
            self.assertIsNotNone(e.mutated_peptide)
            annotation_names = [a.name for a in e.neofox_annotations.annotations]
            observed_mixmhcpred_annotations = \
                observed_mixmhcpred_annotations or \
                "MixMHCpred_score" in annotation_names or "MixMHC2pred_score" in annotation_names
            if with_all_epitopes:
                # they do have extra annotations
                self.assertTrue(
                    "dissimilarity_score" in annotation_names, msg="Annotations: {}".format(annotation_names))
                self.assertTrue(
                    "IEDB_Immunogenicity" in annotation_names, msg="Annotations: {}".format(annotation_names))
            else:
                # they do not have extra annotations
                self.assertFalse(
                    "dissimilarity_score" in annotation_names, msg="Annotations: {}".format(annotation_names))
                self.assertFalse(
                    "IEDB_Immunogenicity" in annotation_names, msg="Annotations: {}".format(annotation_names))
        # not all epitopes may have results for MixMHCpred
        self.assertTrue(observed_mixmhcpred_annotations)
