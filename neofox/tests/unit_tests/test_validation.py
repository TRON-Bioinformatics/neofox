from unittest import TestCase

from neofox.exceptions import NeofoxDataValidationException
from neofox.model.neoantigen import Gene, Neoantigen, Patient, MhcAllele, MhcOne, MhcOneGeneName, MhcOneGene, Zygosity, \
    MhcTwo, MhcTwoName, MhcTwoGeneName, MhcTwoGene, MhcTwoMolecule
from neofox.model.validation import ModelValidator


class TestModelValidator(TestCase):

    def test_bad_type_raises_exception(self):
        gene = Gene(
            gene="BRCA2",
            transcript_identifier=12345,        # this should be a string instead of an integer
            assembly="hg19")
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate, gene)

        neoantigen = Neoantigen(
            patient_identifier="1234",
            rna_expression="0.45")              # this should be a float
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate, neoantigen)

        patient = Patient(
            identifier="1234",
            is_rna_available="Richtig")            # this should be a boolean
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate, patient)

        # TODO: make validation capture this data types errors!
        gene = Gene(
            gene="BRCA2",
            transcript_identifier=["12345"],    # this should be a string instead of a list of strings
            assembly="hg19")
        ModelValidator.validate(gene)

    def test_good_data_does_not_raise_exceptions(self):
        gene = Gene(
            gene="BRCA2",
            transcript_identifier="12345",
            assembly="hg19")
        ModelValidator.validate(gene)

        neoantigen = Neoantigen(
            patient_identifier="1234",
            rna_expression=0.45)
        ModelValidator.validate(neoantigen)

        patient = Patient(
            identifier="1234",
            is_rna_available=True)
        ModelValidator.validate(patient)

    def test_enum_with_wrong_value(self):
        # fails when a string is passed to the enum
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate, MhcOne(name="SOMETHING"))
        # FIXME: does not fail if any integer is passed!
        ModelValidator.validate(MhcOne(name=9))

    def test_mhc_i_allele_validation(self):
        # adds the star
        self._assert_allele_validation(expected="HLA-A*01:01", allele=MhcAllele(name="HLA-A01:01"))
        # adds the HLA-
        self._assert_allele_validation(expected="HLA-A*01:01", allele=MhcAllele(name="A01:01"))
        # adds the colon to homogenise representation
        self._assert_allele_validation(expected="HLA-A*01:01", allele=MhcAllele(name="HLA-A0101"))
        # does not modify an originally good representation
        self._assert_allele_validation(expected="HLA-A*01:01", allele=MhcAllele(name="HLA-A*01:01"))
        # removes further information
        self._assert_allele_validation(expected="HLA-A*01:01", allele=MhcAllele(name="HLA-A01:01:02:03N"))
        self._assert_allele_validation(expected="HLA-A*01:01", allele=MhcAllele(name="HLA-A01:01:02N"))
        self._assert_allele_validation(expected="HLA-A*01:01", allele=MhcAllele(name="HLA-A01:01N"))
        self._assert_allele_validation(expected="HLA-A*01:01", allele=MhcAllele(gene="A", group="01", protein="01"))

    def _assert_allele_validation(self, allele, expected):
        validated_allele = ModelValidator.validate_mhc_allele_representation(allele)
        self.assertEqual(expected, validated_allele.name)

    def test_mhc_ii_allele_validation(self):
        # add the star
        self._assert_allele_validation(expected="HLA-DPB1*01:01", allele=MhcAllele(name="HLA-DPB101:01"))
        # adds the HLA-
        self._assert_allele_validation(expected="HLA-DPB1*01:01", allele=MhcAllele(name="DPB1*01:01"))
        # adds the colon to homogenise representation
        self._assert_allele_validation(expected="HLA-DPA1*01:01", allele=MhcAllele(name="HLA-DPA10101"))
        # does not reove the star
        self._assert_allele_validation(expected="HLA-DPA1*01:01", allele=MhcAllele(name="HLA-DPA1*0101"))
        # removes further information
        self._assert_allele_validation(expected="HLA-DPA1*01:01", allele=MhcAllele(name="HLA-DPA101:01:02:03N"))
        self._assert_allele_validation(expected="HLA-DPA1*01:01", allele=MhcAllele(name="HLA-DPA101:01:02N"))
        self._assert_allele_validation(expected="HLA-DPB1*01:01", allele=MhcAllele(name="HLA-DPB101:01"))
        self._assert_allele_validation(expected="HLA-DPA1*01:01",
                                       allele=MhcAllele(gene="DPA1", group="01", protein="01"))

    def test_invalid_mhc_i_alleles(self):
        # P gene is not valid
        self._assert_invalid_allele(MhcAllele(name="HLA-P01:01"))
        # serotype 1 is not valid
        self._assert_invalid_allele(MhcAllele(name="HLA-A1:01"))
        # no protein information
        self._assert_invalid_allele(MhcAllele(name="HLA-A01"))
        # bad protein format
        self._assert_invalid_allele(MhcAllele(name="HLA-A01:ABC"))
        # wrong organism, only human
        self._assert_invalid_allele(MhcAllele(name="GOGO-A01:01"))
        # no gene
        self._assert_invalid_allele(MhcAllele(name="HLA-0123456"))
        # nonsense
        self._assert_invalid_allele(MhcAllele(name="nonsense"))
        # missing protein
        self._assert_invalid_allele(MhcAllele(gene="A", group="01"))
        # bad protein
        self._assert_invalid_allele(MhcAllele(gene="A", group="01", protein="NaN"))
        # bad group
        self._assert_invalid_allele(MhcAllele(gene="A", group="NaN", protein="01"))
        # non existing gene
        self._assert_invalid_allele(MhcAllele(gene="Z", group="01", protein="01"))
        # MHC I non classical
        self._assert_invalid_allele(MhcAllele(gene="E", group="01", protein="01"))
        self._assert_invalid_allele(MhcAllele(gene="F", group="01", protein="01"))
        self._assert_invalid_allele(MhcAllele(gene="G", group="01", protein="01"))

    def _assert_invalid_allele(self, allele):
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_mhc_allele_representation, allele)

    def test_invalid_mhc_ii_alleles(self):
        # P gene is not valid
        self._assert_invalid_allele(MhcAllele(name="HLA-DPR01:01"))
        # serotype 1 is not valid
        self._assert_invalid_allele(MhcAllele(name="HLA-DPA11:01"))
        # no protein information
        self._assert_invalid_allele(MhcAllele(name="HLA-DPA101"))
        # bad protein format
        self._assert_invalid_allele(MhcAllele(name="HLA-DPA101:ABC"))
        # wrong organism, only human
        self._assert_invalid_allele(MhcAllele(name="GOGO-DPA101:01"))
        # no gene
        self._assert_invalid_allele(MhcAllele(name="HLA-0123456"))
        # nonsense
        self._assert_invalid_allele(MhcAllele(name="nonsense"))
        # missing protein
        self._assert_invalid_allele(MhcAllele(gene="DPA1", group="01"))
        # bad protein
        self._assert_invalid_allele(MhcAllele(gene="DPA1", group="01", protein="NaN"))
        # bad group
        self._assert_invalid_allele(MhcAllele(gene="DPA1", group="NaN", protein="01"))
        # non existing gene
        self._assert_invalid_allele(MhcAllele(gene="DPA1ZZZZ", group="01", protein="01"))

    def test_valid_mhc_i_genotype(self):
        self._assert_valid_patient(patient=Patient(
            identifier="123",
            mhc_one=[MhcOne(name=MhcOneGeneName.A, gene=MhcOneGene(
                        name=MhcOneGeneName.A, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-A01:01"), MhcAllele(name="HLA-A01:02")])),
                     MhcOne(name=MhcOneGeneName.B, gene=MhcOneGene(
                         name=MhcOneGeneName.B, zygosity=Zygosity.HOMOZYGOUS,
                         alleles=[MhcAllele(name="HLA-B01:01")])),
                     MhcOne(name=MhcOneGeneName.C, gene=MhcOneGene(
                         name=MhcOneGeneName.C, zygosity=Zygosity.HEMIZYGOUS,
                         alleles=[MhcAllele(name="HLA-C01:01")]))
                     ]
        ))

    def test_invalid_mhc_i_genotype(self):
        # 3 alleles for gene A
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_one=[MhcOne(
                name=MhcOneGeneName.A, gene=MhcOneGene(
                    name=MhcOneGeneName.A, zygosity=Zygosity.HOMOZYGOUS,
                    alleles=[MhcAllele(name="HLA-A01:01"), MhcAllele(name="HLA-A01:02"), MhcAllele(name="HLA-A01:03")]
                )
            )]
        ))
        # 2 alleles for homozygous gene
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_one=[MhcOne(
                name=MhcOneGeneName.A, gene=MhcOneGene(
                    name=MhcOneGeneName.A, zygosity=Zygosity.HOMOZYGOUS,
                    alleles=[MhcAllele(name="HLA-A01:01"), MhcAllele(name="HLA-A01:02")]
                )
            )]
        ))
        # 1 alleles for heterozygous gene
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_one=[MhcOne(
                name=MhcOneGeneName.A, gene=MhcOneGene(
                    name=MhcOneGeneName.A, zygosity=Zygosity.HETEROZYGOUS,
                    alleles=[MhcAllele(name="HLA-A01:01")]
                )
            )]
        ))
        # 1 alleles for hemizygous gene
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_one=[MhcOne(
                name=MhcOneGeneName.A, gene=MhcOneGene(
                    name=MhcOneGeneName.A, zygosity=Zygosity.HEMIZYGOUS,
                    alleles=[MhcAllele(name="HLA-A01:01"), MhcAllele(name="HLA-A01:02")]
                )
            )]
        ))
        # alleles referring to a different gene
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_one=[MhcOne(
                name=MhcOneGeneName.A, gene=MhcOneGene(
                    name=MhcOneGeneName.A, zygosity=Zygosity.HETEROZYGOUS,
                    alleles=[MhcAllele(name="HLA-B01:01"), MhcAllele(name="HLA-B01:02")]
                )
            )]
        ))
        # MHC and gene referring to different entities
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_one=[MhcOne(
                name=MhcOneGeneName.A, gene=MhcOneGene(
                    name=MhcOneGeneName.B, zygosity=Zygosity.HETEROZYGOUS,
                    alleles=[MhcAllele(name="HLA-A01:01"), MhcAllele(name="HLA-A01:02")]
                )
            )]
        ))

    def _assert_invalid_patient(self, patient):
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)

    def _assert_valid_patient(self, patient):
        ModelValidator.validate_patient(patient)

    def test_valid_mhc_ii_genotype(self):
        self._assert_valid_patient(patient=Patient(
            identifier="123",
            mhc_two=[MhcTwo(
                name=MhcTwoName.DQ, genes=[
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQA1*01:01"), MhcAllele(name="HLA-DQA1*01:02")],
                    ),
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQB1, zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQB1*01:01")]
                    )
                ], molecules=[
                    MhcTwoMolecule(name="HLA-DQA1*01:01-DQB1*01:01"),
                    MhcTwoMolecule(name="HLA-DQA1*01:02-DQB1*01:01")
                ])
            ])
        )
        self._assert_valid_patient(patient=Patient(
            identifier="123",
            mhc_two=[MhcTwo(
                name=MhcTwoName.DQ, genes=[
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQA1*01:01"), MhcAllele(name="HLA-DQA1*01:02")],
                    ),
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQB1, zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQB1*01:01")]
                    )
                ], molecules=[
                    MhcTwoMolecule(
                        alpha_chain=MhcAllele(name="HLA-DQA1*01:01"), beta_chain=MhcAllele(name="HLA-DQB1*01:01")),
                    MhcTwoMolecule(
                        alpha_chain=MhcAllele(name="HLA-DQA1*01:02"), beta_chain=MhcAllele(name="HLA-DQB1*01:01"))
                ])
            ])
        )

    def test_invalid_mhc_ii_genotype(self):
        # 3 alleles for gene
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_two=[MhcTwo(
                name=MhcTwoName.DQ, genes=[
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQA1*01:01"), MhcAllele(name="HLA-DQA1*01:02"),
                                 MhcAllele(name="HLA-DQA1*01:03")]
                    ),
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQB1, zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQB1*01:01")]
                    )
                ])
            ])
        )
        # 2 alleles for homozygous gene
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_two=[MhcTwo(
                name=MhcTwoName.DQ, genes=[
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQA1*01:01"), MhcAllele(name="HLA-DQA1*01:02")]
                    ),
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQB1, zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQB1*01:01"), MhcAllele(name="HLA-DQB1*01:02")]
                    )
                ])
            ])
        )
        # 1 alleles for heterozygous gene
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_two=[MhcTwo(
                name=MhcTwoName.DQ, genes=[
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQA1*01:01")]
                    ),
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQB1, zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQB1*01:01")]
                    )
                ])
            ])
        )
        # 1 alleles for hemizygous gene
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_two=[MhcTwo(
                name=MhcTwoName.DQ, genes=[
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQA1*01:01")]
                    ),
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQB1, zygosity=Zygosity.HEMIZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQB1*01:01")]
                    )
                ])
            ])
        )
        # alleles referring to a different gene
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_two=[MhcTwo(
                name=MhcTwoName.DQ, genes=[
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQB1*01:01")]
                    ),
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQB1, zygosity=Zygosity.HEMIZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQA1*01:01")]
                    )
                ])
            ])
        )
        # MHC and gene referring to different entities
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_two=[MhcTwo(
                name=MhcTwoName.DP, genes=[
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQA1*01:01")]
                    ),
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQB1, zygosity=Zygosity.HEMIZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQB1*01:01")]
                    )
                ])
            ])
        )
        # Molecules refer to alleles not in the MHC
        self._assert_invalid_patient(patient=Patient(
            identifier="123",
            mhc_two=[MhcTwo(
                name=MhcTwoName.DQ, genes=[
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQA1, zygosity=Zygosity.HETEROZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQA1*01:01"), MhcAllele(name="HLA-DQA1*01:02")],
                    ),
                    MhcTwoGene(
                        name=MhcTwoGeneName.DQB1, zygosity=Zygosity.HOMOZYGOUS,
                        alleles=[MhcAllele(name="HLA-DQB1*01:01")]
                    )
                ], molecules=[
                    MhcTwoMolecule(name="HLA-DQA1*01:04-DQB1*01:01"),
                    MhcTwoMolecule(name="HLA-DQA1*01:08-DQB1*01:01")
                ])
            ])
        )

    def test_empty_patient_identifier(self):
        patient = Patient(identifier=None)
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        patient = Patient(identifier="")
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)
        patient = Patient(identifier="   ")
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient, patient)

    def test_bad_is_rna_available(self):
        ModelValidator.validate_patient(Patient(identifier="123", is_rna_available=True))
        ModelValidator.validate_patient(Patient(identifier="123", is_rna_available=False))
        self.assertRaises(NeofoxDataValidationException, ModelValidator.validate_patient,
                          Patient(identifier="123", is_rna_available="False"))
