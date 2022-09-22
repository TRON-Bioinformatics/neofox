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
import betterproto
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein
from Bio.Data import IUPACData
from neofox.exceptions import NeofoxDataValidationException
from logzero import logger
from neofox.model.mhc_parser import HLA_MOLECULE_PATTERN, HLA_DR_MOLECULE_PATTERN, \
    ALLELE_PATTERN_BY_ORGANISM, H2_MOLECULE_PATTERN
from neofox.model.neoantigen import (
    Neoantigen,
    Patient,
    Mhc2Name,
    Mhc2GeneName,
    Zygosity,
    Mhc2,
    Mhc2Isoform,
    MhcAllele,
    Mhc1, Mhc1Name, PredictedEpitope
)
from neofox.references.references import ORGANISM_HOMO_SAPIENS, MHC_I_GENES_BY_ORGANISM, MHC_II_GENES_BY_ORGANISM, \
    ORGANISM_MUS_MUSCULUS

EXTERNAL_ANNOTATIONS_NAME = "External"
FIELD_VAF_DNA = "VAF_in_tumor"
FIELD_VAF_RNA = "VAF_in_RNA"
FIELD_TRANSCRIPT_EXPRESSION = "transcript_expression"
FIELD_GENE = "gene"
FIELD_WILD_TYPE_XMER = "[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)"
FIELD_MUTATED_XMER = "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)"


GENES_BY_MOLECULE = {
    Mhc2Name.DR: [Mhc2GeneName.DRB1],
    Mhc2Name.DP: [Mhc2GeneName.DPA1, Mhc2GeneName.DPB1],
    Mhc2Name.DQ: [Mhc2GeneName.DQA1, Mhc2GeneName.DQB1],
    Mhc2Name.H2A_molecule: [Mhc2GeneName.H2A],
    Mhc2Name.H2E_molecule: [Mhc2GeneName.H2E]
}


class ModelValidator(object):

    @staticmethod
    def validate(model: betterproto.Message):
        # TODO: make this method capture appropriately validation issues when dealing with int and float
        try:
            model.__bytes__()
        except Exception as e:
            raise NeofoxDataValidationException(e)

    @staticmethod
    def validate_neoantigen(neoantigen: Neoantigen):

        # checks format consistency first
        ModelValidator.validate(neoantigen)

        try:
            assert neoantigen.patient_identifier is not None and len(neoantigen.patient_identifier) > 0, \
                "A patient identifier is missing. Please provide patientIdentifier in the input file"

            # checks mutation
            assert neoantigen.mutated_xmer is not None and len(neoantigen.mutated_xmer) > 0, \
                "Missing mutated peptide sequence in input (mutatedXmer) "

            for aa in neoantigen.mutated_xmer:
                ModelValidator._validate_aminoacid(aa)

            # avoids this validation when there is no wild type
            if neoantigen.wild_type_xmer:
                for aa in neoantigen.wild_type_xmer:
                    ModelValidator._validate_aminoacid(aa)

            assert neoantigen.position is not None and neoantigen.position != "", \
                "The position of the mutation is empty, please use EpitopeHelper.mut_position_xmer_seq() to fill it"

            # check the expression values
            ModelValidator._validate_expression_values(neoantigen)
        except AssertionError as e:
            logger.error(neoantigen.to_json(indent=3))
            raise NeofoxDataValidationException(e)

    @staticmethod
    def validate_neoepitope(neoepitope: PredictedEpitope, organism: str):

        # checks format consistency first
        ModelValidator.validate(neoepitope)

        try:
            has_patient_id = neoepitope.patient_identifier is not None and len(neoepitope.patient_identifier) > 0
            has_mhc_i = ModelValidator.is_mhci_epitope(neoepitope)
            has_mhc_ii = ModelValidator.is_mhcii_epitope(neoepitope)

            assert has_patient_id or has_mhc_i or has_mhc_ii, \
                "A patient identifier is missing for a neoepitope without MHC-I allele or MHC-II isoform."

            assert not (has_mhc_i and has_mhc_ii), \
                "Neoepitopes can only be associated to either an MHC-I allele or MHC-II isoform"

            # checks peptides
            for aa in neoepitope.mutated_peptide:
                ModelValidator._validate_aminoacid(aa)
            has_wt_peptide = neoepitope.wild_type_peptide is not None and neoepitope.wild_type_peptide != ""
            if has_wt_peptide:
                for aa in neoepitope.wild_type_peptide:
                    ModelValidator._validate_aminoacid(aa)

            # check lengths according to MHC I or II
            length_mutated_peptide = len(neoepitope.mutated_peptide)
            if has_mhc_i:
                ModelValidator.validate_mhc_allele_representation(neoepitope.allele_mhc_i, organism=organism)
                assert ModelValidator.is_mhci_peptide_length_valid(length_mutated_peptide), \
                    "Mutated MHC-I peptide has a non supported length of {}".format(length_mutated_peptide)
            elif has_mhc_ii:
                ModelValidator.validate_mhc2_isoform_representation(neoepitope.isoform_mhc_i_i, organism=organism)
                assert ModelValidator.is_mhcii_peptide_length_valid(length_mutated_peptide), \
                    "Mutated MHC-II peptide has a non supported length of {}".format(length_mutated_peptide)
            else:
                assert ModelValidator.is_mhci_peptide_length_valid(length_mutated_peptide) or \
                       ModelValidator.is_mhcii_peptide_length_valid(length_mutated_peptide), \
                    "Mutated peptide has a non supported length of {}".format(length_mutated_peptide)

            if has_wt_peptide:
                length_wt_peptide = len(neoepitope.wild_type_peptide)
                if has_mhc_i:
                    assert ModelValidator.is_mhci_peptide_length_valid(length_wt_peptide), \
                        "Mutated MHC-I peptide has a non supported length of {}".format(length_wt_peptide)
                elif has_mhc_ii:
                    assert ModelValidator.is_mhcii_peptide_length_valid(length_wt_peptide), \
                        "Mutated MHC-II peptide has a non supported length of {}".format(length_wt_peptide)
                else:
                    assert ModelValidator.is_mhci_peptide_length_valid(length_wt_peptide) or \
                           ModelValidator.is_mhcii_peptide_length_valid(length_wt_peptide), \
                        "Mutated peptide has a non supported length of {}".format(length_wt_peptide)

            # check the expression values
            ModelValidator._validate_expression_values(neoepitope)
        except AssertionError as e:
            logger.error(neoepitope.to_json(indent=3))
            raise NeofoxDataValidationException(e)

    @staticmethod
    def is_mhcii_peptide_length_valid(length_mutated_peptide):
        return 9 <= length_mutated_peptide <= 20000

    @staticmethod
    def is_mhci_peptide_length_valid(length_mutated_peptide):
        return 8 <= length_mutated_peptide <= 14

    @staticmethod
    def is_mhcii_epitope(neoepitope):
        return neoepitope.isoform_mhc_i_i is not None and neoepitope.isoform_mhc_i_i.name != ''

    @staticmethod
    def is_mhci_epitope(neoepitope):
        return neoepitope.allele_mhc_i is not None and neoepitope.allele_mhc_i.name != ''

    @staticmethod
    def validate_patient(patient: Patient, organism=ORGANISM_HOMO_SAPIENS):

        # checks format consistency first
        ModelValidator.validate(patient)

        try:
            # checks that patient id is not empty considering white spaces

            patient_id = patient.identifier.strip() if patient.identifier else patient.identifier
            assert patient_id is not None and patient_id != "", "A patient identifier is missing"
            assert patient.identifier == patient.identifier.strip(), \
                "Patient identifier contains white spaces at start or end: {}".format(patient.identifier)

            # checks MHC I
            if patient.mhc1:
                for m in patient.mhc1:
                    ModelValidator._validate_mhc1(m, organism=organism)
            # checks MHC II
            if patient.mhc2:
                for m in patient.mhc2:
                    ModelValidator._validate_mhc2(m, organism=organism)

        except AssertionError as e:
            logger.error(patient.to_json(indent=3))
            raise NeofoxDataValidationException(e)

    @staticmethod
    def _validate_mhc1(mhc1: Mhc1, organism: str):
        assert mhc1.name in MHC_I_GENES_BY_ORGANISM.get(organism), "Invalid MHC I name"
        assert mhc1.zygosity in Zygosity, "Invalid zygosity"
        alleles = mhc1.alleles
        if mhc1.zygosity in [Zygosity.HOMOZYGOUS, Zygosity.HEMIZYGOUS]:
            assert (
                len(alleles) == 1
            ), "A homozygous or hemizygous gene must have 1 allele and not {}".format(len(alleles))
        elif mhc1.zygosity == Zygosity.HETEROZYGOUS:
            assert (
                len(alleles) == 2
            ), "A heterozygous gene must have 2 alleles and not {}".format(
                len(alleles)
            )
        elif mhc1.zygosity == Zygosity.LOSS:
            assert (
                len(alleles) == 0
            ), "A lost gene must have 0 alleles and not {}".format(len(alleles))
        for allele in alleles:
            ModelValidator.validate_mhc_allele_representation(allele, organism=organism)
            assert (
                allele.gene == mhc1.name.name
            ), "The allele referring to gene {} is inside gene {}".format(
                allele.gene, mhc1.name.name
            )

    @staticmethod
    def _validate_mhc2(mhc2: Mhc2, organism: str):
        assert mhc2.name in Mhc2Name, "Invalid MHC II name"
        genes = mhc2.genes
        for gene in genes:
            assert gene.name in Mhc2GeneName, "Invalid gene name from MHC II"
            assert gene.name in GENES_BY_MOLECULE.get(
                mhc2.name
            ), "Gene {} referring to isoform {}".format(gene.name, mhc2.name.name)
            assert gene.zygosity in Zygosity, "Invalid zygosity"
            alleles = gene.alleles
            if gene.zygosity in [Zygosity.HOMOZYGOUS, Zygosity.HEMIZYGOUS]:
                assert (
                        len(alleles) == 1
                ), "A homozygous or hemizygous gene must have 1 allele and not {}".format(
                    len(alleles)
                )
            elif gene.zygosity == Zygosity.HETEROZYGOUS:
                assert (
                        len(alleles) == 2
                ), "A heterozygous gene must have 2 alleles and not {}".format(
                    len(alleles)
                )
            elif gene.zygosity == Zygosity.LOSS:
                assert (
                        len(alleles) == 0
                ), "A lost gene must have 0 alleles and not {}".format(len(alleles))
            for allele in alleles:
                ModelValidator.validate_mhc_allele_representation(allele, organism)
                assert allele.gene == gene.name.name, \
                    "The allele {} is inside gene {}".format(allele.name, gene.name.name)
        for isoform in mhc2.isoforms:
            ModelValidator.validate_mhc2_isoform_representation(isoform, organism)
            if mhc2.name != Mhc2Name.DR:
                assert isoform.alpha_chain.name in [
                    a.name for g in genes for a in g.alleles
                ], "Alpha chain allele not present in the list of alleles"
            if mhc2.name not in [Mhc2Name.H2A_molecule, Mhc2Name.H2E_molecule]:
                assert isoform.beta_chain.name in [
                    a.name for g in genes for a in g.alleles
                ], "Beta chain allele not present in the list of alleles"
        return mhc2

    @staticmethod
    def _validate_expression_values(object: Neoantigen or PredictedEpitope):
        assert (
                object.rna_expression is None or object.rna_expression >= 0
        ), "RNA expression should be a positive integer or zero {}".format(
            object.rna_expression
        )
        ModelValidator._validate_vaf(object.dna_variant_allele_frequency)
        ModelValidator._validate_vaf(object.rna_variant_allele_frequency)

    @staticmethod
    def _validate_vaf(vaf):
        assert (
            vaf is None or vaf == -1.0 or 0.0 <= vaf <= 1.0
        ), "VAF should be a decimal number in the range [0.0, 1.0], or else -1.0 for missing values. " \
           "Provided value {}".format(vaf)

    @staticmethod
    def _validate_aminoacid(aminoacid):
        assert aminoacid is not None, "Amino acid field cannot be empty"
        assert isinstance(aminoacid, str), "Amino acid has to be a string"
        # this chunk is unused but let's leave in case it is handy in the future
        if len(aminoacid) == 3:
            assert (
                aminoacid in IUPACData.protein_letters_3to1_extended.keys()
            ), "Non existing 3 letter amino acid {}".format(aminoacid)
            assert aminoacid != "X", "Unknown amino acid X is not supported. Please, remove neoantigens containing an X."
            aminoacid = IUPACData.protein_letters_3to1_extended.get(aminoacid)
        if len(aminoacid) == 1:
            assert (
                aminoacid in ExtendedIUPACProtein.letters
            ), "Non existing aminoacid {}".format(aminoacid)
        else:
            assert False, "Invalid aminoacid {}".format(aminoacid)


    @staticmethod
    def has_peptide_rare_amino_acids(peptide: str):
        has_rare_amino_acid = False
        for aa in peptide:
            has_rare_amino_acid |= aa not in list(IUPACData.protein_letters_3to1.values())
        return has_rare_amino_acid

    @staticmethod
    def validate_mhc_allele_representation(allele: MhcAllele, organism: str):
        try:
            allele_pattern = ALLELE_PATTERN_BY_ORGANISM.get(organism)
            valid_genes = [g.name for g in MHC_I_GENES_BY_ORGANISM.get(organism) + MHC_II_GENES_BY_ORGANISM.get(organism)]

            assert allele_pattern.match(allele.name) is not None, \
                "Allele name does not match expected pattern: {}".format(allele.name)
            assert allele.gene in valid_genes, "MHC gene {} not from classic MHC for organism {}".format(
                allele.gene, organism)
            assert isinstance(allele.protein, str), \
                "The field protein in MHC allele model has the value {} and wrong type but must be a character " \
                "instead of {}".format(allele.protein, type(allele.protein))
            if organism == ORGANISM_HOMO_SAPIENS:
                assert isinstance(allele.group, str), \
                    "The field group in MHC allele model has the value {} and wrong type but must be a character " \
                    "instead of {}".format(allele.group, type(allele.group))
            elif organism == ORGANISM_MUS_MUSCULUS:
                assert allele.group is None or allele.group == "", \
                    "Provided group for H2 allele"
            else:
                raise NeofoxDataValidationException("Not supported organism {}".format(organism))
        except AssertionError as e:
            logger.error(allele.to_json(indent=3))
            raise NeofoxDataValidationException(e)

    @staticmethod
    def validate_mhc2_isoform_representation(isoform: Mhc2Isoform, organism: str):
        try:
            if organism == ORGANISM_HOMO_SAPIENS:
                match_molecule = HLA_MOLECULE_PATTERN.match(isoform.name)
                match_single_allele = HLA_DR_MOLECULE_PATTERN.match(isoform.name)
                assert match_molecule or match_single_allele, "MHC II isoform not following molecule pattern"
                ModelValidator.validate_mhc_allele_representation(isoform.beta_chain, organism)
                if match_molecule:
                    # the DR molecule does not have alpha chain
                    ModelValidator.validate_mhc_allele_representation(isoform.alpha_chain, organism)
            elif organism == ORGANISM_MUS_MUSCULUS:
                match = H2_MOLECULE_PATTERN.match(isoform.name)
                if match:
                    ModelValidator.validate_mhc_allele_representation(isoform.alpha_chain, organism)
                    #ModelValidator.validate_mhc_allele_representation(isoform.beta_chain, organism)
                else:
                    raise NeofoxDataValidationException(
                        "Transformed MHC II molecule name does not match H2 isoform pattern {}".format(isoform.name))
            else:
                raise NeofoxDataValidationException("Not supported organism {}".format(organism))

        except AssertionError as e:
            logger.error(isoform.to_json(indent=3))
            raise NeofoxDataValidationException(e)

    @staticmethod
    def validate_mhc1_gene(allele: MhcAllele):
        assert allele.gene in Mhc1Name.__members__, \
            "MHC I allele is not valid {} at {}".format(allele.gene, allele.full_name)

    @staticmethod
    def validate_mhc2_gene(allele: MhcAllele):
        assert allele.gene in Mhc2GeneName.__members__, \
            "MHC II allele is not valid {} at {}".format(
                allele.gene, allele.full_name) if allele.full_name != "" else "Gene from MHC II allele is empty"
