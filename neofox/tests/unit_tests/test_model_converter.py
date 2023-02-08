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
from neofox.model.conversion import ModelConverter
from neofox.model.neoantigen import (
    Neoantigen,
    Patient,
    Annotation,
    Annotations,
    Zygosity,
    Mhc2Name,
)
from neofox.model.factories import MhcFactory, NeoantigenFactory
from neofox.model.validation import ModelValidator
from neofox.references.references import ORGANISM_HOMO_SAPIENS
from neofox.tests.fake_classes import FakeHlaDatabase, FakeH2Database
from neofox.tests.tools import get_random_neoantigen


class ModelConverterTest(TestCase):

    def setUp(self) -> None:
        self.hla_database = FakeHlaDatabase()
        self.h2_database = FakeH2Database()

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
        csv_data = ModelConverter._objects2dataframe([neoantigen])
        self.assertIsNotNone(csv_data)
        self.assertIsInstance(csv_data, pd.DataFrame)
        self.assertEqual(
            neoantigen.dna_variant_allele_frequency,
            csv_data.iloc[0].dnaVariantAlleleFrequency,
        )

    def test_model2csv2model(self):
        neoantigen = get_random_neoantigen()
        csv_data = ModelConverter._objects2dataframe([neoantigen])
        neoantigen2 = ModelConverter._neoantigens_csv2objects(csv_data)[0]
        neoantigen.external_annotations = None
        neoantigen2.external_annotations = None
        neoantigen.neofox_annotations = None
        neoantigen2.neofox_annotations = None
        self.assertEqual(neoantigen, neoantigen2)

    def test_neoantigen_annotations(self):
        annotations = Annotations()
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

    def test_candidate_neoantigens2model_with_dot_in_column_name(self):
        candidate_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_with_dot_in_column_name.txt"
        )
        with open(candidate_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = ModelConverter().parse_candidate_file(candidate_file)
        for n in neoantigens:
            self.assertTrue(sum(a.name == "my.annotation.with.dots" for a in n.external_annotations) > 0)

    def test_candidate_neoantigens2model_with_mutation_in_column_name(self):
        candidate_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_with_mutation_in_column_name.txt"
        )
        with open(candidate_file) as f:
            self.count_lines = len(f.readlines())
        neoantigens = ModelConverter().parse_candidate_file(candidate_file)
        for n in neoantigens:
            self.assertTrue(sum(a.name == "mutation" for a in n.external_annotations) > 0)

    def _assert_external_annotations(
        self, expected_external_annotations, neoantigen_annotations, non_nullable_annotations=[]
    ):
        self.assertEqual(
            len(expected_external_annotations),
            len(neoantigen_annotations),
        )
        found_external_annotation = {}
        for a in neoantigen_annotations:
            self.assertIsInstance(a, Annotation)
            self.assertNotEmpty(a.name)
            if a.name in expected_external_annotations:
                found_external_annotation[a.name] = True
                if a.name in non_nullable_annotations:
                    self.assertIsNotNone(a.value, "{} should not be none".format(a.name))
        for ea in expected_external_annotations:
            self.assertTrue(found_external_annotation.get(ea, False), "Not found {} external annotation".format(ea))

    def test_csv_neoantigens2model(self):
        neoantigens_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_model.txt"
        )
        data = pd.read_csv(neoantigens_file, sep="\t")
        neoantigens = ModelConverter._neoantigens_csv2objects(data)
        self.assertEqual(5, len(neoantigens))
        for n in neoantigens:
            self.assertTrue(isinstance(n, Neoantigen))
            self.assertNotEmpty(n.patient_identifier)
            self.assertNotEmpty(n.gene)
            self.assertNotEmpty(n.mutated_xmer)
            self.assertNotEmpty(n.wild_type_xmer)
            self.assertIsNotNone(n.position)

            # test external annotations
            self._assert_external_annotations(
                expected_external_annotations=["external_annotation_1", "external_annotation_2"],
                neoantigen_annotations=n.external_annotations,
                non_nullable_annotations=["external_annotation_1", "external_annotation_2"]
            )
            for a in n.external_annotations:
                if a.name == "external_annotation_1":
                    self.assertEqual(a.value, "blah1")
                elif a.name == "external_annotation_2":
                    self.assertEqual(a.value, "blah2")

    def test_json_neoantigens2model(self):
        neoantigens_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_json.json"
        )
        neoantigens = ModelConverter.parse_neoantigens_json_file(neoantigens_file)
        self.assertEqual(5, len(neoantigens))
        for n in neoantigens:
            self.assertTrue(isinstance(n, Neoantigen))
            self.assertNotEmpty(n.patient_identifier)
            self.assertNotEmpty(n.rna_expression)
            self.assertNotEmpty(n.rna_variant_allele_frequency)
            self.assertNotEmpty(n.dna_variant_allele_frequency)
            self.assertNotEmpty(n.position)

    def assertNotEmpty(self, value):
        self.assertIsNotNone(value)
        self.assertNotEqual(value, "")

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

    def test_patients_without_mhc2(self):
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/alleles.Pt29_without_mhc2.csv"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 2)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertEqual(3, len(patients[0].mhc1))
        self.assertEqual(6, len([a for m in patients[0].mhc1 for a in m.alleles]))
        self.assertEqual(0, len(patients[0].mhc2))

    def test_patients_csv_file2model_mouse(self):
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/alleles.Pt29_mouse.csv"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.h2_database)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 1)
        self.assertIsInstance(patients[0], Patient)
        self.assertEqual(patients[0].identifier, "Pt29")
        self.assertEqual(3, len(patients[0].mhc1))
        self.assertEqual(6, len([a for m in patients[0].mhc1 for a in m.alleles]))
        self.assertEqual(2, len(patients[0].mhc2))
        self.assertEqual(
            3, len([a for m in patients[0].mhc2 for g in m.genes for a in g.alleles])
        )

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
        for m in patients[0].mhc1:
            self.assertEqual(m.zygosity, Zygosity.LOSS)
        self.assertEqual(3, len(patients[0].mhc2))
        self.assertEqual(
            9, len([a for m in patients[0].mhc2 for g in m.genes for a in g.alleles])
        )

    def test_patients_csv_file2model_without_mhc2(self):
        patients_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/balachandran_supplementary_table1_patients.tsv"
        )
        patients = ModelConverter.parse_patients_file(patients_file, self.hla_database)
        self.assertIsNotNone(patients)
        self.assertIsInstance(patients, list)
        self.assertTrue(len(patients) == 58)

    def test_annotations2short_wide_df(self):

        neoantigens = [
            Neoantigen(
                wild_type_xmer="AAAAAAA",
                mutated_xmer="AAACAAA",
                position=[],
                neofox_annotations=Annotations(
                    annotations=[
                        Annotation(name="this_name", value="this_value"),
                        Annotation(name="that_name", value="that_value"),
                        Annotation(name="diese_name", value="diese_value"),
                        Annotation(name="das_name", value="das_value"),
                    ]
                )
            ),
            Neoantigen(
                wild_type_xmer="AAAGAAA",
                mutated_xmer="AAAZAAA",
                position=[1, 2, 3],
                neofox_annotations=Annotations(
                    annotations=[
                        Annotation(name="this_name", value="0"),
                        Annotation(name="that_name", value="1"),
                        Annotation(name="diese_name", value="2"),
                        Annotation(name="das_name", value="3"),
                    ],
                )
            ),
        ]
        df = ModelConverter.annotations2neoantigens_table(neoantigens=neoantigens)
        self.assertEqual(df.shape[0], 2)
        self.assertEqual(df.shape[1], 13)
        self.assertEqual(0, df[df["position"].transform(lambda x: isinstance(x, list))].shape[0])

    def test_parse_mhc1_heterozygous_alleles(self):
        mhc1s = MhcFactory.build_mhc1_alleles(
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

    def test_parse_mhc1_heterozygous_alleles_mouse(self):
        mhc1s = MhcFactory.build_mhc1_alleles(
            [
                "H2Kd",
                "H2Kp",
                "H2Dd",
                "H2Dp",
                "H2Ld",
                "H2Lp"
            ], self.h2_database
        )
        self.assertEqual(3, len(mhc1s))
        for mhc1 in mhc1s:
            self.assertEqual(Zygosity.HETEROZYGOUS, mhc1.zygosity)
            self.assertEqual(2, len(mhc1.alleles))

    def test_parse_mhc1_homozygous_alleles(self):
        mhc1s = MhcFactory.build_mhc1_alleles(
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

    def test_parse_mhc1_homozygous_alleles_mouse(self):
        mhc1s = MhcFactory.build_mhc1_alleles(
            [
                "H2Kd",
                "H2Kd",
                "H2Dd",
                "H2Dd",
                "H2Ld",
                "H2Ld"
            ], self.h2_database
        )
        self.assertEqual(3, len(mhc1s))
        for mhc1 in mhc1s:
            self.assertEqual(Zygosity.HOMOZYGOUS, mhc1.zygosity)
            self.assertEqual(1, len(mhc1.alleles))

    def test_parse_mhc1_hemizygous_alleles(self):
        mhc1s = MhcFactory.build_mhc1_alleles(
            ["HLA-A*01:01", "HLA-B*07:02", "HLA-C*01:02"], self.hla_database
        )
        self.assertEqual(3, len(mhc1s))
        for mhc1 in mhc1s:
            self.assertEqual(Zygosity.HEMIZYGOUS, mhc1.zygosity)
            self.assertEqual(1, len(mhc1.alleles))

    def test_parse_mhc1_hemizygous_alleles_mouse(self):
        mhc1s = MhcFactory.build_mhc1_alleles(
            [
                "H2Kd",
                "H2Dd",
                "H2Ld"
            ], self.h2_database
        )
        self.assertEqual(3, len(mhc1s))
        for mhc1 in mhc1s:
            self.assertEqual(Zygosity.HEMIZYGOUS, mhc1.zygosity)
            self.assertEqual(1, len(mhc1.alleles))

    def test_parse_mhc1_loss_alleles(self):
        mhc1s = MhcFactory.build_mhc1_alleles([], self.hla_database)
        self.assertEqual(0, len(mhc1s))

    def test_parse_mhc1_loss_alleles_mouse(self):
        mhc1s = MhcFactory.build_mhc1_alleles([], self.h2_database)
        self.assertEqual(0, len(mhc1s))

    def test_parse_mhc1_bad_format_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc1_alleles,
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

    def test_parse_mhc1_bad_format_fails_mouse(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc1_alleles,
            [
                "c.12345C>G",
                "H2Kd",
                "H2Dd",
                "H2Ld"
            ],
            self.h2_database
        )

    def test_parse_mhc1_too_many_alleles_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc1_alleles,
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

    def test_parse_mhc1_too_many_alleles_fails_mouse(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc1_alleles,
            [
                "H2Kd",
                "H2Kd",
                "H2Kd",
                "H2Dd",
                "H2Ld"
            ],
            self.h2_database
        )

    def test_parse_mhc1_bad_gene_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc1_alleles,
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

    def test_parse_mhc1_bad_gene_fails_mouse(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc1_alleles,
            [
                "H2Kd",
                "H2Kd",
                "H2Pd",
                "H2Dd",
                "H2Ld"
            ],
            self.h2_database
        )

    def test_parse_mhc1_non_existing_allele_does_not_fail(self):
        mhc1s = MhcFactory.build_mhc1_alleles(
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

    def test_parse_mhc1_non_existing_allele_does_not_fail_mouse(self):
        mhc1s = MhcFactory.build_mhc1_alleles(
            [
                "H2Kd",
                "H2Kz",     # this one does not exist
                "H2Dd",
                "H2Ld"
            ], self.h2_database
        )
        self.assertEqual(3, len(mhc1s))

    def test_parse_mhc2_heterozygous_alleles(self):
        mhc2s = MhcFactory.build_mhc2_alleles(
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

    def test_parse_mhc2_heterozygous_alleles_mouse(self):
        mhc2s = MhcFactory.build_mhc2_alleles(
            [
                "H2Ad",
                "H2Ap",
                "H2Ed",
                "H2Ep"
            ],
            self.h2_database
        )
        self.assertEqual(2, len(mhc2s))
        for mhc2 in mhc2s:
            self.assertEqual(1, len(mhc2.genes))
            for gene in mhc2.genes:
                self.assertEqual(2, len(gene.alleles))
                self.assertEqual(Zygosity.HETEROZYGOUS, gene.zygosity)
            self.assertEqual(2, len(mhc2.isoforms))
            self._assert_isoforms(mhc2)

    def test_parse_mhc2_homozygous_alleles(self):
        mhc2s = MhcFactory.build_mhc2_alleles(
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

    def test_parse_mhc2_homozygous_alleles_mouse(self):
        mhc2s = MhcFactory.build_mhc2_alleles(
            [
                "H2Ad",
                "H2Ad",
                "H2Ed",
                "H2Ed"
            ],
            self.h2_database
        )
        self.assertEqual(2, len(mhc2s))
        for mhc2 in mhc2s:
            self.assertEqual(1, len(mhc2.genes))
            for gene in mhc2.genes:
                self.assertEqual(1, len(gene.alleles))
                self.assertEqual(Zygosity.HOMOZYGOUS, gene.zygosity)
            self.assertEqual(1, len(mhc2.isoforms))
            self._assert_isoforms(mhc2)

    def test_parse_mhc2_hemizygous_alleles(self):
        mhc2s = MhcFactory.build_mhc2_alleles(
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

    def test_parse_mhc2_hemizygous_alleles_mouse(self):
        mhc2s = MhcFactory.build_mhc2_alleles(
            [
                "H2Ad",
                "H2Ed",
            ],
            self.h2_database
        )
        self.assertEqual(2, len(mhc2s))
        for mhc2 in mhc2s:
            self.assertEqual(1, len(mhc2.genes))
            for gene in mhc2.genes:
                self.assertEqual(1, len(gene.alleles))
                self.assertEqual(Zygosity.HEMIZYGOUS, gene.zygosity)
            self.assertEqual(1 if mhc2.name == Mhc2Name.DR else 1, len(mhc2.isoforms))
            self._assert_isoforms(mhc2)

    def test_parse_mhc2_hetero_and_homozygous_alleles(self):
        mhc2s = MhcFactory.build_mhc2_alleles(
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
        mhc2s = MhcFactory.build_mhc2_alleles([], self.hla_database)
        self.assertEqual(0, len(mhc2s))

    def test_parse_mhc2_loss_mouse(self):
        mhc2s = MhcFactory.build_mhc2_alleles([], self.h2_database)
        self.assertEqual(0, len(mhc2s))

    def test_parse_mhc2_bad_format_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc2_alleles,
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

    def test_parse_mhc2_bad_format_fails_mouse(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc2_alleles,
            [
                "c.12345C>G",
                "H2Ad",
                "H2Ed",
            ],
            self.h2_database
        )

    def test_parse_mhc2_too_many_alleles_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc2_alleles,
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

    def test_parse_mhc2_too_many_alleles_fails_mouse(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc2_alleles,
            [
                "H2Ad",
                "H2Ae",
                "H2Af",
                "H2Ed",
            ],
            self.h2_database
        )

    def test_parse_mhc2_bad_gene_fails(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc2_alleles,
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

    def test_parse_mhc2_bad_gene_fails_mouse(self):
        self.assertRaises(
            NeofoxDataValidationException,
            MhcFactory.build_mhc2_alleles,
            [
                "H2Ad",
                "H2Ae",
                "H2Zf",
                "H2Ed",
            ],
            self.h2_database
        )

    def test_parse_mhc2_non_existing_allele_does_not_fail(self):
        mhc2s = MhcFactory.build_mhc2_alleles(
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

    def test_parse_mhc2_non_existing_allele_does_not_fail_mouse(self):
        mhc2s = MhcFactory.build_mhc2_alleles(
            [
                "H2Ad",
                "H2Az", # this one does not exist
                "H2Ed",
                "H2Ed",
            ],
            self.h2_database
        )
        self.assertEqual(2, len(mhc2s))

    def test_candidate_neoepitopes2model(self):
        candidate_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_neoepitopes.txt"
        )
        with open(candidate_file) as f:
            self.count_lines = len(f.readlines())
        neoepitopes = ModelConverter().parse_candidate_neoepitopes_file(candidate_file, self.hla_database)
        self.assertIsNotNone(neoepitopes)
        self.assertEqual(self.count_lines -1, len(neoepitopes))
        for n in neoepitopes:
            ModelValidator.validate_neoepitope(n, ORGANISM_HOMO_SAPIENS)

    def test_candidate_neoepitopes2model_with_patients(self):
        candidate_file = pkg_resources.resource_filename(
            neofox.tests.__name__, "resources/test_data_neoepitopes_with_patients.txt"
        )

        with open(candidate_file) as f:
            self.count_lines = len(f.readlines())

        neoepitopes = ModelConverter().parse_candidate_neoepitopes_file(candidate_file, self.hla_database)
        self.assertIsNotNone(neoepitopes)
        self.assertEqual(self.count_lines -1, len(neoepitopes))
        for n in neoepitopes:
            ModelValidator.validate_neoepitope(n, ORGANISM_HOMO_SAPIENS)

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
