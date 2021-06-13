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
from unittest import TestCase
import pkg_resources
import pandas as pd
from neofox.exceptions import NeofoxDataValidationException

import neofox.tests
from neofox.model.conversion import ModelConverter, EXTERNAL_ANNOTATIONS_NAME
from neofox.model.neoantigen import (
    Neoantigen,
    Mutation,
    Patient,
    Annotation,
    NeoantigenAnnotations,
    Zygosity,
    Mhc2Name,
)
from neofox.model.conversion import ModelValidator
from neofox.tests.fake_classes import FakeHlaDatabase
from neofox.tests.unit_tests.tools import get_random_neoantigen


class ModelConverterTest(TestCase):

    def setUp(self) -> None:
        self.hla_database = FakeHlaDatabase()

    def test_model2json(self):
        neoantigens = [get_random_neoantigen() for _ in range(5)]
        json_data = [n.to_json() for n in neoantigens]
        self.assertIsInstance(json_data, list)
        self.assertEqual(5, len(json_data))
        neoantigens2 = [Neoantigen().from_json(j) for j in json_data]
        self._assert_lists_equal(neoantigens, neoantigens2)

    def test_model2dict(self):
        neoantigens = [get_random_neoantigen() for _ in range(5)]
        json_data = [n.to_dict() for n in neoantigens]
        self.assertIsInstance(json_data, list)
        self.assertEqual(5, len(json_data))
        neoantigens2 = [Neoantigen().from_dict(j) for j in json_data]
        self._assert_lists_equal(neoantigens, neoantigens2)

    def test_model2csv(self):
        neoantigen = get_random_neoantigen()
        csv_data = ModelConverter.object2series(neoantigen)
        self.assertIsNotNone(csv_data)
        self.assertIsInstance(csv_data, pd.Series)
        self.assertEqual(
            neoantigen.dna_variant_allele_frequency,
            csv_data.dna_variant_allele_frequency,
        )

    def test_model2flat_dict(self):
        neoantigen = get_random_neoantigen()
        flat_dict = ModelConverter.object2flat_dict(neoantigen)
        self.assertIsNotNone(flat_dict)
        self.assertEqual(
            neoantigen.dna_variant_allele_frequency,
            flat_dict["dna_variant_allele_frequency"],
        )
        self.assertEqual(
            neoantigen.mutation.mutated_xmer, flat_dict["mutation.mutated_xmer"]
        )

    def test_model2csv2model(self):
        neoantigen = get_random_neoantigen()
        csv_data = ModelConverter.object2series(neoantigen)
        neoantigen2 = ModelConverter.neoantigens_csv2object(csv_data)
        self.assertEqual(neoantigen, neoantigen2)

    def test_neoantigen_annotations(self):
        neoantigen = get_random_neoantigen()
        annotations = NeoantigenAnnotations()
        annotations.neoantigen_identifier = (
            ModelValidator.generate_neoantigen_identifier(neoantigen)
        )
        annotations.annotations = [
            Annotation(name="string_annotation", value="blabla"),
            Annotation(name="integer_annotation", value=1),
            Annotation(name="float_annotation", value=1.1),
        ]
        annotations_dict = annotations.to_dict()
        self.assertTrue(len(annotations_dict.get("annotations")) == 3)
        self.assertEqual(annotations_dict.get("annotations")[0].get("value"), "blabla")
        # this does not fail, but it will fail validation
        self.assertEqual(annotations_dict.get("annotations")[1].get("value"), 1)
        self.assertEqual(annotations_dict.get("annotations")[2].get("value"), 1.1)
        self.assertTrue(
            annotations_dict.get("neoantigenIdentifier"),
            ModelValidator.generate_neoantigen_identifier(neoantigen),
        )

    def test_candidate_neoantigens2model(self):
        canidate_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data.txt"
        )
        with open(canidate_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens, external_annotations = ModelConverter().parse_candidate_file(
            canidate_file
        )
        self.assertIsNotNone(neoantigens)
        self.assertEqual(self.count_lines -1, len(neoantigens))
        for n in neoantigens:
            self.assertIsInstance(n, Neoantigen)
            self.assertIsInstance(n.mutation, Mutation)
            self.assertTrue(n.gene is not None and len(n.gene) > 0)
            self.assertTrue(
                n.mutation.mutated_xmer is not None and len(n.mutation.mutated_xmer) > 1
            )
            self.assertTrue(
                n.mutation.wild_type_xmer is not None
                and len(n.mutation.wild_type_xmer) > 1
            )
            self.assertTrue(
                n.mutation.position is not None and len(n.mutation.position) >= 1
            )
            self.assertTrue(
                n.rna_variant_allele_frequency is None
                or n.rna_variant_allele_frequency == -1
                or (0 <= n.rna_variant_allele_frequency <= 1)
            )
            self.assertTrue(n.rna_expression is None or n.rna_expression >= 0)
            self.assertTrue(0 <= n.dna_variant_allele_frequency <= 1)

        # test external annotations
        self._assert_external_annotations(
            expected_number_external_annotations=44,
            external_annotations=external_annotations,
        )

    def _assert_external_annotations(
        self, expected_number_external_annotations, external_annotations
    ):
        for neoantigen_annotation in external_annotations:
            self.assertIsInstance(neoantigen_annotation, NeoantigenAnnotations)
            self.assertNotEmpty(neoantigen_annotation.neoantigen_identifier)
            self.assertEqual(neoantigen_annotation.annotator, EXTERNAL_ANNOTATIONS_NAME)
            self.assertEqual(
                expected_number_external_annotations,
                len(neoantigen_annotation.annotations),
            )
            for a in neoantigen_annotation.annotations:
                self.assertIsInstance(a, Annotation)
                self.assertNotEmpty(a.name)
                if a.name == "VAF_RNA_limits":
                    self.assertIsNone(a.value)
                if a.name == "MHC_II_epitope_(WT)":
                    self.assertIsNotNone(a.value)

    def test_csv_neoantigens2model(self):
        neoantigens_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_model.txt"
        )
        data = pd.read_csv(neoantigens_file, sep="\t")
        neoantigens, external_annotations = ModelConverter.parse_neoantigens_dataframe(
            data
        )
        self.assertEqual(5, len(neoantigens))
        for n in neoantigens:
            self.assertTrue(isinstance(n, Neoantigen))
            self.assertNotEmpty(n.mutation)
            self.assertNotEmpty(n.patient_identifier)
            self.assertNotEmpty(n.gene)
            self.assertTrue(isinstance(n.mutation, Mutation))
            self.assertNotEmpty(n.mutation.mutated_xmer)
            self.assertNotEmpty(n.mutation.wild_type_xmer)
            self.assertIsNotNone(n.mutation.position)

        # test external annotations
        self._assert_external_annotations(
            expected_number_external_annotations=2,
            external_annotations=external_annotations,
        )

    def test_json_neoantigens2model(self):
        neoantigens_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_json.json"
        )
        neoantigens = ModelConverter.parse_neoantigens_json_file(neoantigens_file)
        self.assertEqual(5, len(neoantigens))
        for n in neoantigens:
            self.assertTrue(isinstance(n, Neoantigen))
            self.assertNotEmpty(n.mutation)
            self.assertNotEmpty(n.patient_identifier)
            self.assertNotEmpty(n.rna_expression)
            self.assertNotEmpty(n.rna_variant_allele_frequency)
            self.assertNotEmpty(n.dna_variant_allele_frequency)
            self.assertTrue(isinstance(n.mutation, Mutation))
            self.assertNotEmpty(n.mutation.position)

    def assertNotEmpty(self, value):
        self.assertIsNotNone(value)
        self.assertNotEqual(value, "")

    def test_overriding_patient_id(self):
        candidate_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data.txt"
        )
        with open(candidate_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens, external_annotations = ModelConverter().parse_candidate_file(
            candidate_file, patient_id="patientX"
        )
        for n in neoantigens:
            self.assertEqual(n.patient_identifier, "patientX")
        neoantigens, external_annotations = ModelConverter().parse_candidate_file(
            candidate_file
        )
        for n in neoantigens:
            self.assertEqual(n.patient_identifier, "Ptx")

    def test_patients_csv_file2model(self):
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/alleles.Pt29.csv"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertEqual(3, len(patients[0].mhc1))
        self.assertEqual(6, len([a for m in patients[0].mhc1 for a in m.alleles]))
        self.assertEqual(3, len(patients[0].mhc2))
        self.assertEqual(
            9, len([a for m in patients[0].mhc2 for g in m.genes for a in g.alleles])
        )
        self.assertEqual(patients[0].is_rna_available, False)

    def test_patients_csv_file2model2(self):
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/patient.Pt29.csv"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertEqual(3, len(patients[0].mhc1))
        self.assertEqual(6, len([a for m in patients[0].mhc1 for a in m.alleles]))
        self.assertEqual(3, len(patients[0].mhc2))
        self.assertEqual(
            9, len([a for m in patients[0].mhc2 for g in m.genes for a in g.alleles])
        )
        self.assertEqual(patients[0].is_rna_available, True)

    def test_patients_csv_file2model3(self):
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_patient_info.txt"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Ptx")
        self.assertEqual(3, len(patients[0].mhc1))
        self.assertEqual(6, len([a for m in patients[0].mhc1 for a in m.alleles]))
        self.assertEqual(3, len(patients[0].mhc2))
        self.assertEqual(
            10, len([a for m in patients[0].mhc2 for g in m.genes for a in g.alleles])
        )
        self.assertTrue(
            "HLA-A*03:01" in [a.name for m in patients[0].mhc1 for a in m.alleles]
        )
        self.assertTrue(
            "HLA-DQA1*04:01"
            in [a.name for m in patients[0].mhc2 for g in m.genes for a in g.alleles]
        )
        self.assertTrue(patients[0].is_rna_available)

    def test_patients_csv_file2model_without_mhc1(self):
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/patient.Pt29.without_mhc1.csv"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertIsNone(patients[0].mhc1)
        self.assertEqual(3, len(patients[0].mhc2))
        self.assertEqual(
            9, len([a for m in patients[0].mhc2 for g in m.genes for a in g.alleles])
        )
        self.assertEqual(patients[0].is_rna_available, True)

    def test_patients_csv_file2model_without_mhc2(self):
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/balachandran_supplementary_table1_patients.tsv"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 58)

    def test_annotations2short_wide_df(self):
        annotations = [
            NeoantigenAnnotations(
                neoantigen_identifier="12345",
                annotations=[
                    Annotation(name="this_name", value="this_value"),
                    Annotation(name="that_name", value="that_value"),
                    Annotation(name="diese_name", value="diese_value"),
                    Annotation(name="das_name", value="das_value"),
                ],
            ),
            NeoantigenAnnotations(
                neoantigen_identifier="6789",
                annotations=[
                    Annotation(name="this_name", value="0"),
                    Annotation(name="that_name", value="1"),
                    Annotation(name="diese_name", value="2"),
                    Annotation(name="das_name", value="3"),
                ],
            ),
        ]
        neoantigens = [
            Neoantigen(
                identifier="12345",
                mutation=Mutation(wild_type_xmer="AAAAAAA", mutated_xmer="AAACAAA", position=[]),
            ),
            Neoantigen(
                identifier="6789",
                mutation=Mutation(wild_type_xmer="AAAGAAA", mutated_xmer="AAAZAAA", position=[1, 2, 3]),
            ),
        ]
        df = ModelConverter.annotations2short_wide_table(
            neoantigen_annotations=annotations, neoantigens=neoantigens
        )
        self.assertEqual(df.shape[0], 2)
        self.assertEqual(df.shape[1], 13)
        self.assertEqual(0, df[df["mutation.position"].transform(lambda x: isinstance(x, list))].shape[0])

        df_annotations = ModelConverter.annotations2tall_skinny_table(annotations)
        self.assertEqual(df_annotations.shape[0], 8)
        self.assertEqual(df_annotations.shape[1], 3)

    def test_parse_mhc1_heterozygous_alleles(self):
        mhc1s = ModelConverter.parse_mhc1_alleles(
            [
                "HLA-A*01:01",
                "HLA-A*01:02",
                "HLA-B*07:02",
                "HLA-B*07:03",
                "HLA-C*01:02",
                "HLA-C*01:03",
            ], self.hla_database
        )
        self.assertEqual(3, len(mhc1s))
        for mhc1 in mhc1s:
            self.assertEqual(Zygosity.HETEROZYGOUS, mhc1.zygosity)
            self.assertEqual(2, len(mhc1.alleles))

    def test_parse_mhc1_homozygous_alleles(self):
        mhc1s = ModelConverter.parse_mhc1_alleles(
            [
                "HLA-A*01:01",
                "HLA-A*01:01",
                "HLA-B*07:02",
                "HLA-B*07:02",
                "HLA-C*01:02",
                "HLA-C*01:02",
            ], self.hla_database
        )
        self.assertEqual(3, len(mhc1s))
        for mhc1 in mhc1s:
            self.assertEqual(Zygosity.HOMOZYGOUS, mhc1.zygosity)
            self.assertEqual(1, len(mhc1.alleles))

    def test_parse_mhc1_hemizygous_alleles(self):
        mhc1s = ModelConverter.parse_mhc1_alleles(
            ["HLA-A*01:01", "HLA-B*07:02", "HLA-C*01:02"], self.hla_database
        )
        self.assertEqual(3, len(mhc1s))
        for mhc1 in mhc1s:
            self.assertEqual(Zygosity.HEMIZYGOUS, mhc1.zygosity)
            self.assertEqual(1, len(mhc1.alleles))

    def test_parse_mhc1_loss_alleles(self):
        mhc1s = ModelConverter.parse_mhc1_alleles([], self.hla_database)
        self.assertEqual(3, len(mhc1s))
        for mhc1 in mhc1s:
            self.assertEqual(Zygosity.LOSS, mhc1.zygosity)
            self.assertEqual(0, len(mhc1.alleles))

    def test_parse_mhc1_bad_format_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelConverter.parse_mhc1_alleles,
            [
                "c.12345C>G",
                "HLA-A*01:02",
                "HLA-A*01:01",
                "HLA-B*01:02",
                "HLA-C*01:01",
                "HLA-C*01:02",
            ],
            self.hla_database
        )

    def test_parse_mhc1_too_many_alleles_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelConverter.parse_mhc1_alleles,
            [
                "HLA-A*01:01",
                "HLA-A*01:02",
                "HLA-A*01:01",
                "HLA-B*01:02",
                "HLA-C*01:01",
                "HLA-C*01:02",
            ],
            self.hla_database
        )

    def test_parse_mhc1_bad_gene_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelConverter.parse_mhc1_alleles,
            [
                "HLA-A*01:01",
                "HLA-A*01:02",
                "HLA-G*01:01",
                "HLA-B*01:02",
                "HLA-C*01:01",
                "HLA-C*01:02",
            ],
            self.hla_database
        )

    def test_parse_mhc1_non_existing_allele_does_not_fail(self):
        mhc1s = ModelConverter.parse_mhc1_alleles(
            [
                "HLA-A*01:01",
                "HLA-A*01:01",
                "HLA-B*999:01",     # this one does not exist
                "HLA-B*07:02",
                "HLA-C*01:02",
                "HLA-C*01:02",
            ], self.hla_database
        )
        self.assertEqual(3, len(mhc1s))

    def test_parse_mhc2_heterozygous_alleles(self):
        mhc2s = ModelConverter.parse_mhc2_alleles(
            [
                "HLA-DRB1*01:01",
                "HLA-DRB1*01:02",
                "HLA-DPA1*01:03",
                "HLA-DPA1*01:04",
                "HLA-DPB1*02:01",
                "HLA-DPB1*02:02",
                "HLA-DQA1*01:01",
                "HLA-DQA1*01:02",
                "HLA-DQB1*02:01",
                "HLA-DQB1*02:02",
            ],
            self.hla_database
        )
        self.assertEqual(3, len(mhc2s))
        for mhc2 in mhc2s:
            self.assertEqual(1 if mhc2.name == Mhc2Name.DR else 2, len(mhc2.genes))
            for gene in mhc2.genes:
                self.assertEqual(2, len(gene.alleles))
                self.assertEqual(Zygosity.HETEROZYGOUS, gene.zygosity)
            self.assertEqual(2 if mhc2.name == Mhc2Name.DR else 4, len(mhc2.isoforms))
            self._assert_isoforms(mhc2)

    def test_parse_mhc2_homozygous_alleles(self):
        mhc2s = ModelConverter.parse_mhc2_alleles(
            [
                "HLA-DRB1*01:01",
                "HLA-DRB1*01:01",
                "HLA-DPA1*01:03",
                "HLA-DPA1*01:03",
                "HLA-DPB1*02:01",
                "HLA-DPB1*02:01",
                "HLA-DQA1*01:01",
                "HLA-DQA1*01:01",
                "HLA-DQB1*02:01",
                "HLA-DQB1*02:01",
            ],
            self.hla_database
        )
        self.assertEqual(3, len(mhc2s))
        for mhc2 in mhc2s:
            self.assertEqual(1 if mhc2.name == Mhc2Name.DR else 2, len(mhc2.genes))
            for gene in mhc2.genes:
                self.assertEqual(1, len(gene.alleles))
                self.assertEqual(Zygosity.HOMOZYGOUS, gene.zygosity)
            self.assertEqual(1 if mhc2.name == Mhc2Name.DR else 1, len(mhc2.isoforms))
            self._assert_isoforms(mhc2)

    def test_parse_mhc2_hemizygous_alleles(self):
        mhc2s = ModelConverter.parse_mhc2_alleles(
            [
                "HLA-DRB1*01:01",
                "HLA-DPA1*01:03",
                "HLA-DPB1*02:01",
                "HLA-DQA1*01:01",
                "HLA-DQB1*02:01",
            ],
            self.hla_database
        )
        self.assertEqual(3, len(mhc2s))
        for mhc2 in mhc2s:
            self.assertEqual(1 if mhc2.name == Mhc2Name.DR else 2, len(mhc2.genes))
            for gene in mhc2.genes:
                self.assertEqual(1, len(gene.alleles))
                self.assertEqual(Zygosity.HEMIZYGOUS, gene.zygosity)
            self.assertEqual(1 if mhc2.name == Mhc2Name.DR else 1, len(mhc2.isoforms))
            self._assert_isoforms(mhc2)

    def test_parse_mhc2_hetero_and_homozygous_alleles(self):
        mhc2s = ModelConverter.parse_mhc2_alleles(
            [
                "HLA-DRB1*01:01",
                "HLA-DRB1*01:01",
                "HLA-DPA1*01:03",
                "HLA-DPA1*01:03",
                "HLA-DPB1*02:01",
                "HLA-DPB1*02:02",
                "HLA-DQA1*01:01",
                "HLA-DQA1*01:01",
                "HLA-DQB1*02:01",
                "HLA-DQB1*02:02",
            ],
            self.hla_database
        )
        self.assertEqual(3, len(mhc2s))
        for mhc2 in mhc2s:
            self.assertEqual(1 if mhc2.name == Mhc2Name.DR else 2, len(mhc2.genes))
            self.assertEqual(1 if mhc2.name == Mhc2Name.DR else 2, len(mhc2.isoforms))

    def test_parse_mhc2_loss(self):
        mhc2s = ModelConverter.parse_mhc2_alleles([], self.hla_database)
        self.assertEqual(3, len(mhc2s))
        for mhc2 in mhc2s:
            self.assertEqual(1 if mhc2.name == Mhc2Name.DR else 2, len(mhc2.genes))
            for gene in mhc2.genes:
                self.assertEqual(Zygosity.LOSS, gene.zygosity)
                self.assertEqual(0, len(gene.alleles))
            self.assertEqual(0, len(mhc2.isoforms))

    def test_parse_mhc2_bad_format_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelConverter.parse_mhc2_alleles,
            [
                "c.12345C>G",
                "HLA-DRB1*01:01",
                "HLA-DPA1*01:01",
                "HLA-DPA1*01:01",
                "HLA-DPB1*01:01",
                "HLA-DPB1*01:02",
                "HLA-DQA1*01:01",
                "HLA-DQA1*01:01",
                "HLA-DQB1*01:01",
                "HLA-DQB1*01:02",
            ],
            self.hla_database
        )

    def test_parse_mhc2_too_many_alleles_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelConverter.parse_mhc2_alleles,
            [
                "HLA-DRB1*01:01",
                "HLA-DRB1*01:02",
                "HLA-DRB1*01:03",
                "HLA-DPA1*01:01",
                "HLA-DPA1*01:01",
                "HLA-DPB1*01:01",
                "HLA-DPB1*01:02",
                "HLA-DQA1*01:01",
                "HLA-DQA1*01:01",
                "HLA-DQB1*01:01",
            ],
            self.hla_database
        )

    def test_parse_mhc2_bad_gene_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelConverter.parse_mhc2_alleles,
            [
                "HLA-A*01:01",
                "HLA-A*01:02",
                "HLA-G*01:01",
                "HLA-B*01:02",
                "HLA-C*01:01",
                "HLA-C*01:02",
            ],
            self.hla_database
        )

    def test_parse_mhc2_non_existing_allele_does_not_fail(self):
        mhc2s = ModelConverter.parse_mhc2_alleles(
            [
                "HLA-DRB1*999:01",      # this one does not exist
                "HLA-DPA1*01:03",
                "HLA-DPB1*01:01",
                "HLA-DQA1*01:01",
                "HLA-DQB1*02:01",
            ],
            self.hla_database
        )
        self.assertEqual(3, len(mhc2s))

    def _assert_isoforms(self, mhc2):
        for isoform in mhc2.isoforms:
            if mhc2.name == Mhc2Name.DR:
                self.assertEqual("", isoform.alpha_chain.full_name)
            else:
                self.assertIsNot("", isoform.alpha_chain.full_name)
            self.assertIsNotNone(isoform.beta_chain)

    def _assert_lists_equal(self, neoantigens, neoantigens2):
        self.assertEqual(len(neoantigens), len(neoantigens2))
        for (
            n1,
            n2,
        ) in zip(neoantigens, neoantigens2):
            self.assertEqual(n1, n2)
