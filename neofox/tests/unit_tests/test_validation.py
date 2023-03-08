from unittest import TestCase

from neofox.exceptions import NeofoxDataValidationException
from neofox.model.mhc_parser import MhcParser
from neofox.model.neoantigen import (
    Neoantigen,
    Patient,
    Mhc1,
    Mhc1Name,
    Zygosity,
    Mhc2,
    Mhc2Name,
    Mhc2GeneName,
    Mhc2Gene,
    Mhc2Isoform, PredictedEpitope, MhcAllele,
)
from neofox.model.validation import ModelValidator
from neofox.references.references import ORGANISM_HOMO_SAPIENS, ORGANISM_MUS_MUSCULUS
from neofox.tests.fake_classes import FakeHlaDatabase, FakeH2Database


class TestModelValidator(TestCase):

    def setUp(self) -> None:
        self.hla_parser = MhcParser.get_mhc_parser(FakeHlaDatabase())
        self.h2_parser = MhcParser.get_mhc_parser(FakeH2Database())

    def test_bad_type_raises_exception(self):

        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate,
            Neoantigen(
                patient_identifier=1234,  # this should be a string instead of an integer
                rna_expression=0.45,
            ),
        )

        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate,
            Neoantigen(patient_identifier="1234", rna_expression="0.45"),
        )  # this should be a float)

        # TODO: make validation capture this data types errors!
        ModelValidator.validate(
            Neoantigen(
                patient_identifier=[
                    "12345"
                ],  # this should be a string instead of a list of strings
                rna_expression=0.45,
            )
        )

    def test_good_data_does_not_raise_exceptions(self):

        neoantigen = Neoantigen(patient_identifier="1234", rna_expression=0.45)
        ModelValidator.validate(neoantigen)

        patient = Patient(identifier="1234")
        ModelValidator.validate(patient)

    def test_enum_with_wrong_value(self):
        # fails when a string is passed to the enum
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate,
            Mhc1(name="SOMETHING"),
        )
        # FIXME: does not fail if any integer is passed!
        ModelValidator.validate(Mhc1(name=9))

    def test_valid_hla_i_genotype(self):
        self._assert_valid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.A,
                        zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[
                            self.hla_parser.parse_mhc_allele("HLA-A01:01"),
                            self.hla_parser.parse_mhc_allele("HLA-A01:02"),
                        ],
                    ),
                    Mhc1(
                        name=Mhc1Name.B,
                        zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[self.hla_parser.parse_mhc_allele("HLA-B01:01")],
                    ),
                    Mhc1(
                        name=Mhc1Name.C,
                        zygosity=Zygosity.HEMIZYGOUS,
                        alleles=[self.hla_parser.parse_mhc_allele("HLA-C01:01")],
                    ),
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )

    def test_valid_h2_i_genotype(self):
        self._assert_valid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.H2D,
                        zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[
                            self.h2_parser.parse_mhc_allele("H2Dd"),
                            self.h2_parser.parse_mhc_allele("H2Da"),
                        ],
                    ),
                    Mhc1(
                        name=Mhc1Name.H2L,
                        zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[self.h2_parser.parse_mhc_allele("H2Ld")],
                    ),
                    Mhc1(
                        name=Mhc1Name.H2K,
                        zygosity=Zygosity.HEMIZYGOUS,
                        alleles=[self.h2_parser.parse_mhc_allele("H2Kk")],
                    ),
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )

    def test_invalid_hla_i_genotype(self):
        # 3 alleles for gene A
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.A,
                        zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[
                            self.hla_parser.parse_mhc_allele("HLA-A01:01"),
                            self.hla_parser.parse_mhc_allele("HLA-A01:02"),
                            self.hla_parser.parse_mhc_allele("HLA-A01:03"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # 2 alleles for homozygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.A,
                        zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[
                            self.hla_parser.parse_mhc_allele("HLA-A01:01"),
                            self.hla_parser.parse_mhc_allele("HLA-A01:02"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # 1 alleles for heterozygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.A,
                        zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[self.hla_parser.parse_mhc_allele("HLA-A01:01")],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # 1 alleles for hemizygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.A,
                        zygosity=Zygosity.HEMIZYGOUS,
                        alleles=[
                            self.hla_parser.parse_mhc_allele("HLA-A01:01"),
                            self.hla_parser.parse_mhc_allele("HLA-A01:02"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # alleles referring to a different gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.A,
                        zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[
                            self.hla_parser.parse_mhc_allele("HLA-B01:01"),
                            self.hla_parser.parse_mhc_allele("HLA-B01:02"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )

    def test_invalid_h2_i_genotype(self):
        # 3 alleles for gene A
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.H2D,
                        zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[
                            self.h2_parser.parse_mhc_allele("H2Dd"),
                            self.h2_parser.parse_mhc_allele("H2Dk"),
                            self.h2_parser.parse_mhc_allele("H2Dp"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )
        # 2 alleles for homozygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.H2D,
                        zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[
                            self.h2_parser.parse_mhc_allele("H2Dd"),
                            self.h2_parser.parse_mhc_allele("H2Dk"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )
        # 1 alleles for heterozygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.H2D,
                        zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[
                            self.h2_parser.parse_mhc_allele("H2Dd"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )
        # 2 alleles for hemizygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.H2D,
                        zygosity=Zygosity.HEMIZYGOUS,
                        alleles=[
                            self.h2_parser.parse_mhc_allele("H2Dd"),
                            self.h2_parser.parse_mhc_allele("H2Dk"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )
        # alleles referring to a different gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc1=[
                    Mhc1(
                        name=Mhc1Name.H2D,
                        zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[
                            self.h2_parser.parse_mhc_allele("H2Kd"),
                            self.h2_parser.parse_mhc_allele("H2Kk"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )

    def _assert_invalid_patient(self, patient, organism):
        self.assertRaises(
            NeofoxDataValidationException, ModelValidator.validate_patient, patient, organism
        )

    def _assert_valid_patient(self, patient, organism):
        ModelValidator.validate_patient(patient, organism=organism)

    def test_valid_hla_ii_genotype(self):
        patient = Patient(identifier="123", mhc2=[Mhc2(name=Mhc2Name.DQ, genes=[
            Mhc2Gene(name=Mhc2GeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                     alleles=[self.hla_parser.parse_mhc_allele("HLA-DQA1*01:01"),
                              self.hla_parser.parse_mhc_allele("HLA-DQA1*01:02"), ], ),
            Mhc2Gene(name=Mhc2GeneName.DQB1, zygosity=Zygosity.HOMOZYGOUS,
                     alleles=[self.hla_parser.parse_mhc_allele("HLA-DQB1*01:01")], ), ], isoforms=[
            self.hla_parser.parse_mhc2_isoform("HLA-DQA1*01:01-DQB1*01:01"),
            self.hla_parser.parse_mhc2_isoform("HLA-DQA1*01:02-DQB1*01:01"), ], )], )
        self._assert_valid_patient(patient=patient, organism=ORGANISM_HOMO_SAPIENS)

        patient2 = Patient(identifier="123", mhc2=[Mhc2(name=Mhc2Name.DR, genes=[
            Mhc2Gene(name=Mhc2GeneName.DRB1, zygosity=Zygosity.HOMOZYGOUS,
                     alleles=[self.hla_parser.parse_mhc_allele("HLA-DRB1*01:01")],
                     )], isoforms=[self.hla_parser.parse_mhc2_isoform("HLA-DRB1*01:01")], )], )
        self._assert_valid_patient(patient=patient2, organism=ORGANISM_HOMO_SAPIENS)

        patient3 = Patient(identifier="123",
                           mhc2=[
                               Mhc2(
                                   name=Mhc2Name.DQ,
                                   genes=[
                                       Mhc2Gene(name=Mhc2GeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                                                alleles=[
                                                    self.hla_parser.parse_mhc_allele("HLA-DQA1*01:01"),
                                                    self.hla_parser.parse_mhc_allele("HLA-DQA1*01:02"),
                                                ], ),
                                       Mhc2Gene(name=Mhc2GeneName.DQB1, zygosity=Zygosity.HOMOZYGOUS,
                                                alleles=[
                                                    self.hla_parser.parse_mhc_allele("HLA-DQB1*01:01")
                                                ], ), ],
                                   isoforms=[
                                       self.hla_parser.parse_mhc2_isoform("HLA-DQA1*01:01-DQB1*01:01"),
                                       self.hla_parser.parse_mhc2_isoform("HLA-DQA1*01:02-DQB1*01:01"),
                                   ], )], )
        self._assert_valid_patient(patient=patient3, organism=ORGANISM_HOMO_SAPIENS)

    def test_valid_h2_ii_genotype(self):
        patient = Patient(
            identifier="123",
            mhc2=[Mhc2(name=Mhc2Name.H2A_molecule,
                       genes=[
                           Mhc2Gene(name=Mhc2GeneName.H2A, zygosity=Zygosity.HETEROZYGOUS,
                                    alleles=[self.h2_parser.parse_mhc_allele("H2Ad"),
                                             self.h2_parser.parse_mhc_allele("H2Ap"), ], )],
                       isoforms=[
                           self.h2_parser.parse_mhc2_isoform("H2Ad"),
                           self.h2_parser.parse_mhc2_isoform("H2Ap"), ], )], )
        self._assert_valid_patient(patient=patient, organism=ORGANISM_MUS_MUSCULUS)

        patient2 = Patient(
            identifier="123",
            mhc2=[Mhc2(name=Mhc2Name.H2E_molecule,
                       genes=[
                           Mhc2Gene(name=Mhc2GeneName.H2E, zygosity=Zygosity.HOMOZYGOUS,
                                    alleles=[self.h2_parser.parse_mhc_allele("H2Ed")])],
                       isoforms=[self.h2_parser.parse_mhc2_isoform("H2Ed")])])
        self._assert_valid_patient(patient=patient2, organism=ORGANISM_MUS_MUSCULUS)

    def test_invalid_hla_ii_genotype(self):
        # 3 alleles for gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.DQ,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.DQA1,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[
                                    self.hla_parser.parse_mhc_allele("HLA-DQA1*01:01"),
                                    self.hla_parser.parse_mhc_allele("HLA-DQA1*01:02"),
                                    self.hla_parser.parse_mhc_allele("HLA-DQA1*01:03"),
                                ],
                            ),
                            Mhc2Gene(
                                name=Mhc2GeneName.DQB1,
                                zygosity=Zygosity.HOMOZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQB1*01:01")],
                            ),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # 2 alleles for homozygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.DQ,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.DQA1,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[
                                    self.hla_parser.parse_mhc_allele("HLA-DQA1*01:01"),
                                    self.hla_parser.parse_mhc_allele("HLA-DQA1*01:02"),
                                ],
                            ),
                            Mhc2Gene(
                                name=Mhc2GeneName.DQB1,
                                zygosity=Zygosity.HOMOZYGOUS,
                                alleles=[
                                    self.hla_parser.parse_mhc_allele("HLA-DQB1*01:01"),
                                    self.hla_parser.parse_mhc_allele("HLA-DQB1*01:02"),
                                ],
                            ),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # 1 alleles for heterozygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.DQ,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.DQA1,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQA1*01:01")],
                            ),
                            Mhc2Gene(
                                name=Mhc2GeneName.DQB1,
                                zygosity=Zygosity.HOMOZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQB1*01:01")],
                            ),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # 1 alleles for hemizygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.DQ,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.DQA1,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQA1*01:01")],
                            ),
                            Mhc2Gene(
                                name=Mhc2GeneName.DQB1,
                                zygosity=Zygosity.HEMIZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQB1*01:01")],
                            ),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # alleles referring to a different gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.DQ,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.DQA1,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQB1*01:01")],
                            ),
                            Mhc2Gene(
                                name=Mhc2GeneName.DQB1,
                                zygosity=Zygosity.HEMIZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQA1*01:01")],
                            ),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # MHC and gene referring to different entities
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.DP,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.DQA1,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQA1*01:01")],
                            ),
                            Mhc2Gene(
                                name=Mhc2GeneName.DQB1,
                                zygosity=Zygosity.HEMIZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQB1*01:01")],
                            ),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )
        # Molecules refer to alleles not in the MHC
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.DQ,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.DQA1,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[
                                    self.hla_parser.parse_mhc_allele("HLA-DQA1*01:01"),
                                    self.hla_parser.parse_mhc_allele("HLA-DQA1*01:02"),
                                ],
                            ),
                            Mhc2Gene(
                                name=Mhc2GeneName.DQB1,
                                zygosity=Zygosity.HOMOZYGOUS,
                                alleles=[self.hla_parser.parse_mhc_allele("HLA-DQB1*01:01")],
                            ),
                        ],
                        isoforms=[
                            Mhc2Isoform(name="HLA-DQA1*01:04-DQB1*01:01"),
                            Mhc2Isoform(name="HLA-DQA1*01:08-DQB1*01:01"),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_HOMO_SAPIENS
        )

    def test_invalid_h2_ii_genotype(self):
        # 3 alleles for gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.H2E_molecule,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.H2E,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[
                                    self.h2_parser.parse_mhc_allele("H2Ea"),
                                    self.h2_parser.parse_mhc_allele("H2Eb"),
                                    self.h2_parser.parse_mhc_allele("H2Ec"),
                                ],
                            )
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )
        # 2 alleles for homozygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.H2A_molecule,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.H2A,
                                zygosity=Zygosity.HOMOZYGOUS,
                                alleles=[
                                    self.h2_parser.parse_mhc_allele("H2Aa"),
                                    self.h2_parser.parse_mhc_allele("H2Ab"),
                                ],
                            ),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )
        # 1 alleles for heterozygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.H2A_molecule,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.H2A,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[self.h2_parser.parse_mhc_allele("H2Aa")],
                            )
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )
        # 1 alleles for hemizygous gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.H2A_molecule,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.H2A,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[self.h2_parser.parse_mhc_allele("H2Aa")],
                            )
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )
        # alleles referring to a different gene
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.H2A_molecule,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.H2A,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[self.h2_parser.parse_mhc_allele("H2Ea")],
                            ),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )
        # MHC and gene referring to different entities
        self._assert_invalid_patient(
            patient=Patient(
                identifier="123",
                mhc2=[
                    Mhc2(
                        name=Mhc2Name.H2A_molecule,
                        genes=[
                            Mhc2Gene(
                                name=Mhc2GeneName.H2E,
                                zygosity=Zygosity.HETEROZYGOUS,
                                alleles=[self.h2_parser.parse_mhc_allele("H2Ea")],
                            ),
                        ],
                    )
                ],
            ),
            organism=ORGANISM_MUS_MUSCULUS
        )

    def test_empty_patient_identifier(self):
        patient = Patient(identifier=None)
        self.assertRaises(
            NeofoxDataValidationException, ModelValidator.validate_patient, patient, ORGANISM_HOMO_SAPIENS
        )
        patient = Patient(identifier="")
        self.assertRaises(
            NeofoxDataValidationException, ModelValidator.validate_patient, patient, ORGANISM_HOMO_SAPIENS
        )
        patient = Patient(identifier="   ")
        self.assertRaises(
            NeofoxDataValidationException, ModelValidator.validate_patient, patient, ORGANISM_HOMO_SAPIENS
        )

    def test_validate_neoepitope_mhci(self):
        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            wild_type_peptide="DILVIDQTR",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
        )
        ModelValidator.validate_neoepitope(neoepitope, ORGANISM_HOMO_SAPIENS)

    def test_validate_neoepitope_mhci_without_wt(self):
        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
        )
        ModelValidator.validate_neoepitope(neoepitope, ORGANISM_HOMO_SAPIENS)

    def test_validate_neoepitope_mhci_bad_length(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                mutated_peptide="DILVT",  # 5 aa < min 8 aa
                wild_type_peptide="DILVIDQTR",
                allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            ),
            ORGANISM_HOMO_SAPIENS
        )
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                mutated_peptide="DILVTAAAAAAAAAAAAAAAAA",  # 22 aa > max 14 aa
                wild_type_peptide="DILVIDQTR",
                allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            ),
            ORGANISM_HOMO_SAPIENS
        )
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                wild_type_peptide="DILVT",  # 5 aa < min 8 aa
                mutated_peptide="DILVIDQTR",
                allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            ),
            ORGANISM_HOMO_SAPIENS
        )
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                wild_type_peptide="DILVTAAAAAAAAAAAAAAAAA",  # 22 aa > max 14 aa
                mutated_peptide="DILVIDQTR",
                allele_mhc_i=self._get_test_mhci_allele('HLA-A*01:01'),
            ),
            ORGANISM_HOMO_SAPIENS
        )

    def test_validate_neoepitope_mhcii(self):
        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            wild_type_peptide="DILVIDQTR",
            isoform_mhc_i_i=self._get_test_mhcii_isoform('DRB1*01:01'),
        )
        ModelValidator.validate_neoepitope(neoepitope, ORGANISM_HOMO_SAPIENS)

    def test_validate_neoepitope_mhcii_without_wt(self):
        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            isoform_mhc_i_i=self._get_test_mhcii_isoform('DRB1*01:01'),

        )
        ModelValidator.validate_neoepitope(neoepitope, ORGANISM_HOMO_SAPIENS)

    def test_validate_neoepitope_mhcii_bad_length(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                mutated_peptide="DILVT",  # 5 aa < min 8 aa
                wild_type_peptide="DILVIDQTR",
                isoform_mhc_i_i=self._get_test_mhcii_isoform('DRB1*01:01'),
            ),
            ORGANISM_HOMO_SAPIENS
        )
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                wild_type_peptide="DILVT",  # 5 aa < min 8 aa
                mutated_peptide="DILVIDQTR",
                isoform_mhc_i_i=self._get_test_mhcii_isoform('DRB1*01:01'),
            ),
            ORGANISM_HOMO_SAPIENS
        )

    def test_validate_neoepitope_with_patient_id(self):
        neoepitope = PredictedEpitope(
            mutated_peptide="DILVTDQTR",
            wild_type_peptide="DILVIDQTR",
            patient_identifier="123",
        )
        ModelValidator.validate_neoepitope(neoepitope, ORGANISM_HOMO_SAPIENS)

    def test_validate_neoepitope_with_patient_id_bad_length(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                mutated_peptide="DILVT",  # 5 aa < min 8 aa
                wild_type_peptide="DILVIDQTR",
                patient_identifier="123",
            ),
            ORGANISM_HOMO_SAPIENS
        )
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                wild_type_peptide="DILVT",  # 5 aa < min 8 aa
                mutated_peptide="DILVIDQTR",
                patient_identifier="123",
            ),
            ORGANISM_HOMO_SAPIENS
        )

    def test_validate_neoepitope_mhci_bad_allele(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                mutated_peptide="DILVTDQTR",
                allele_mhc_i=MhcAllele(name="something"),
            ),
            ORGANISM_HOMO_SAPIENS
        )

    def test_validate_neoepitope_mhcii_bad_isoform(self):
        self.assertRaises(
            NeofoxDataValidationException,
            ModelValidator.validate_neoepitope,
            PredictedEpitope(
                mutated_peptide="DILVTDQTR",
                isoform_mhc_i_i=Mhc2Isoform(name="something"),
            ),
            ORGANISM_HOMO_SAPIENS
        )

    def _get_test_mhci_allele(self, allele) -> MhcAllele:

        return self.hla_parser.parse_mhc_allele(allele)

    def _get_test_mhcii_isoform(self, isoform) -> Mhc2Isoform:
        return self.hla_parser.parse_mhc2_isoform(isoform)
