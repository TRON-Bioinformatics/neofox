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
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.exceptions import NeofoxDataValidationException
from logzero import logger
from neofox.model.mhc_parser import HLA_MOLECULE_PATTERN, HLA_ALLELE_PATTERN, HLA_DR_MOLECULE_PATTERN, \
    ALLELE_PATTERN_BY_ORGANISM
from neofox.model.neoantigen import (
    Neoantigen,
    Mutation,
    Patient,
    Mhc2Name,
    Mhc2GeneName,
    Zygosity,
    Mhc2,
    Mhc2Isoform,
    MhcAllele,
    Mhc1Name,
    Mhc1
)
from neofox.model.wrappers import get_mhc2_isoform_name
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
}


class ModelValidator(object):
    @staticmethod
    def validate(model: betterproto.Message):
        # TODO: make this method capture appropriately validation issues when dealing with int and float
        try:
            model.__bytes__()
        except Exception as e:
            raise NeofoxDataValidationException(e)

    # TODO: add patient validation: validate GTEx tissue and MHC alleles

    @staticmethod
    def validate_neoantigen(neoantigen: Neoantigen) -> Neoantigen:

        # checks format consistency first
        ModelValidator.validate(neoantigen)

        try:
            assert neoantigen.patient_identifier is not None and len(neoantigen.patient_identifier) > 0, \
                "A patient identifier is missing. Please provide patientIdentifier in the input file"

            # checks mutation
            neoantigen.mutation = ModelValidator._validate_mutation(neoantigen.mutation)

            # check the expression values
            ModelValidator._validate_expression_values(neoantigen)
        except AssertionError as e:
            logger.error(neoantigen.to_json(indent=3))
            raise NeofoxDataValidationException(e)

        return neoantigen

    @staticmethod
    def validate_patient(patient: Patient, organism=ORGANISM_HOMO_SAPIENS) -> Patient:

        # checks format consistency first
        ModelValidator.validate(patient)

        try:
            # checks that patient id is not empty considering white spaces
            patient_id = (
                patient.identifier.strip() if patient.identifier else patient.identifier
            )
            assert (
                patient_id is not None and patient_id != ""
            ), "A patient identifier is missing"
            patient.identifier = patient_id

            # TODO: validate new model with isoforms, genes and alleles
            # checks MHC I
            validated_mhc1s = []
            if patient.mhc1:
                for m in patient.mhc1:
                    validated_mhc1s.append(ModelValidator._validate_mhc1(m, organism=organism))
            patient.mhc1 = validated_mhc1s
            # checks MHC II
            validated_mhc2s = []
            if patient.mhc2:
                for m in patient.mhc2:
                    validated_mhc2s.append(ModelValidator._validate_mhc2(m, organism=organism))
            patient.mhc2 = validated_mhc2s

        except AssertionError as e:
            logger.error(patient.to_json(indent=3))
            raise NeofoxDataValidationException(e)

        return patient

    @staticmethod
    def _validate_mhc1(mhc1: Mhc1, organism: str) -> Mhc1:
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
        validated_alleles = []
        for allele in alleles:
            validated_allele = ModelValidator.validate_mhc_allele_representation(allele, organism=organism)
            validated_alleles.append(validated_allele)
            assert (
                validated_allele.gene == mhc1.name.name
            ), "The allele referring to gene {} is inside gene {}".format(
                validated_allele.gene, mhc1.name.name
            )
        mhc1.alleles = validated_alleles
        return mhc1

    @staticmethod
    def _validate_mhc2(mhc2: Mhc2, organism: str) -> Mhc2:
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
            validated_alleles = []
            for allele in alleles:
                validated_allele = ModelValidator.validate_mhc_allele_representation(
                    allele, organism
                )
                validated_alleles.append(validated_allele)
                assert (
                    validated_allele.gene == gene.name.name
                ), "The allele referring to gene {} is inside gene {}".format(
                    validated_allele.gene, gene.name.name
                )
            gene.alleles = validated_alleles
        isoforms = mhc2.isoforms
        validated_isoforms = []
        for isoform in isoforms:
            validated_isoform = ModelValidator.validate_mhc2_isoform_representation(
                isoform, organism
            )
            validated_isoforms.append(validated_isoform)
            if mhc2.name != Mhc2Name.DR:
                assert validated_isoform.alpha_chain.name in [
                    a.name for g in genes for a in g.alleles
                ], "Alpha chain allele not present in th list of alleles"
            assert validated_isoform.beta_chain.name in [
                a.name for g in genes for a in g.alleles
            ], "Beta chain allele not present in th list of alleles"
        mhc2.isoforms = validated_isoforms
        return mhc2

    @staticmethod
    def _validate_expression_values(neoantigen):
        assert (
            neoantigen.rna_expression is None or neoantigen.rna_expression >= 0
        ), "RNA expression should be a positive integer or zero {}".format(
            neoantigen.rna_expression
        )
        ModelValidator._validate_vaf(neoantigen.dna_variant_allele_frequency)
        ModelValidator._validate_vaf(neoantigen.rna_variant_allele_frequency)

    @staticmethod
    def _validate_mutation(mutation: Mutation) -> Mutation:
        assert mutation.mutated_xmer is not None and len(mutation.mutated_xmer) > 0, \
            "Missing mutated peptide sequence in input (mutation.mutatedXmer) "
        mutation.mutated_xmer = "".join(
            [ModelValidator._validate_aminoacid(aa) for aa in mutation.mutated_xmer]
        )
        # avoids this validation when there is no wild type
        if mutation.wild_type_xmer:
            mutation.wild_type_xmer = "".join(
                [
                    ModelValidator._validate_aminoacid(aa)
                    for aa in mutation.wild_type_xmer
                ]
            )
            mutation.position = EpitopeHelper.mut_position_xmer_seq(mutation=mutation)
        return mutation

    @staticmethod
    def _validate_vaf(vaf):
        assert (
            vaf is None or vaf == -1.0 or 0.0 <= vaf <= 1.0
        ), "VAF should be a decimal number in the range [0.0, 1.0], or else -1.0 for missing values. " \
           "Provided value {}".format(vaf)

    @staticmethod
    def _validate_aminoacid(aminoacid):
        assert aminoacid is not None, "Amino acid field cannot be empty"
        aminoacid = aminoacid.strip()
        assert isinstance(aminoacid, str), "Amino acid has to be a string"
        # this chunk is unused but let's leave in case it is handy in the future
        if len(aminoacid) == 3:
            assert (
                aminoacid in IUPACData.protein_letters_3to1_extended.keys()
            ), "Non existing 3 letter amino acid {}".format(aminoacid)
            assert aminoacid != "X", "Unknown amino acid X is not supported. Please, remove neoantigens containing an X."
            aminoacid = IUPACData.protein_letters_3to1_extended.get(aminoacid)
        if len(aminoacid) == 1:
            aminoacid = aminoacid.upper()
            assert (
                aminoacid in ExtendedIUPACProtein.letters
            ), "Non existing aminoacid {}".format(aminoacid)
        else:
            assert False, "Invalid aminoacid {}".format(aminoacid)
        return aminoacid


    @staticmethod
    def has_peptide_rare_amino_acids(peptide: str):
        has_rare_amino_acid = False
        for aa in peptide:
            has_rare_amino_acid |= aa not in list(IUPACData.protein_letters_3to1.values())
        return has_rare_amino_acid

    @staticmethod
    def validate_mhc_allele_representation(allele: MhcAllele, organism: str) -> MhcAllele:
        try:
            full_name = None
            allele_pattern = ALLELE_PATTERN_BY_ORGANISM.get(organism)
            valid_genes = MHC_I_GENES_BY_ORGANISM.get(organism) + MHC_II_GENES_BY_ORGANISM.get(organism)

            if organism == ORGANISM_HOMO_SAPIENS:
                if allele.full_name:
                    # infers gene, group and protein from the name
                    match = allele_pattern.match(allele.full_name)
                    assert (
                            match is not None
                    ), "Allele does not match HLA allele pattern {}".format(allele.name)
                    gene = match.group(1)
                    group = match.group(2)
                    protein = match.group(3)
                    full_name = allele.full_name
                elif allele.name:
                    # infers gene, group and protein from the name
                    match = allele_pattern.match(allele.name)
                    assert (
                            match is not None
                    ), "Allele does not match HLA allele pattern {}".format(allele.name)
                    gene = match.group(1)
                    group = match.group(2)
                    protein = match.group(3)
                elif allele.gene and allele.group and allele.protein:
                    # infers name from gene, group and protein
                    gene = allele.gene
                    group = allele.group
                    protein = allele.protein
                else:
                    logger.error(allele.to_json(indent=3))
                    raise NeofoxDataValidationException(
                        "HLA allele missing required fields, either name or gene, group and protein must be provided"
                    )
                # builds the final allele representation and validates it just in case
                name = "HLA-{gene}*{serotype}:{protein}".format(
                    gene=gene, serotype=group, protein=protein
                )
            elif organism == ORGANISM_MUS_MUSCULUS:
                group = ""
                if allele.full_name:
                    # infers gene, group and protein from the name
                    match = allele_pattern.match(allele.full_name)
                    assert (
                            match is not None
                    ), "Allele does not match H2 allele pattern {}".format(allele.name)
                    gene = match.group(1)
                    protein = match.group(2)
                    full_name = allele.full_name
                elif allele.name:
                    # infers gene, group and protein from the name
                    match = allele_pattern.match(allele.name)
                    assert (
                            match is not None
                    ), "Allele does not match H2 allele pattern {}".format(allele.name)
                    gene = match.group(1)
                    protein = match.group(2)
                elif allele.gene and allele.protein:
                    # infers name from gene, group and protein
                    gene = allele.gene
                    protein = allele.protein
                else:
                    logger.error(allele.to_json(indent=3))
                    raise NeofoxDataValidationException(
                        "H2 allele missing required fields, either name or gene and protein must be provided"
                    )
                # builds the final allele representation and validates it just in case
                name = "{gene}{protein}".format(gene=gene, protein=protein)
            else:
                raise NeofoxDataValidationException("Not supported organism {}".format(organism))

            assert gene in [g.name for g in valid_genes], "Gene not from classic MHC: {}".format(gene)
            match = allele_pattern.match(name)
            assert (
                    match is not None
            ), "Allele does not match allele pattern {}".format(name)
        except AssertionError as e:
            logger.error(allele.to_json(indent=3))
            raise NeofoxDataValidationException(e)

        return MhcAllele(
            full_name=full_name if full_name else name,
            name=name,
            gene=gene,
            group=group,
            protein=protein,
        )

    @staticmethod
    def validate_mhc2_isoform_representation(isoform: Mhc2Isoform, organism: str) -> Mhc2Isoform:
        try:
            if organism == ORGANISM_HOMO_SAPIENS:
                if isoform.name:
                    # infers alpha and beta chains
                    match = HLA_MOLECULE_PATTERN.match(isoform.name)
                    if match:
                        alpha_chain = ModelValidator.validate_mhc_allele_representation(
                            MhcAllele(name=match.group(1)), organism
                        )
                        beta_chain = ModelValidator.validate_mhc_allele_representation(
                            MhcAllele(name=match.group(2)), organism
                        )
                    else:
                        match = HLA_DR_MOLECULE_PATTERN.match(isoform.name)
                        assert (
                                match is not None
                        ), "Molecule does not match HLA isoform pattern {}".format(
                            isoform.name
                        )
                        alpha_chain = MhcAllele()
                        beta_chain = ModelValidator.validate_mhc_allele_representation(
                            MhcAllele(name=match.group(1)), organism
                        )
                elif isoform.alpha_chain and isoform.beta_chain:
                    # infers name from gene, group and protein
                    alpha_chain = ModelValidator.validate_mhc_allele_representation(
                        isoform.alpha_chain, organism
                    )
                    beta_chain = ModelValidator.validate_mhc_allele_representation(
                        isoform.beta_chain, organism
                    )
                else:
                    logger.error(isoform.to_json(indent=3))
                    raise NeofoxDataValidationException(
                        "HLA isoform missing required fields"
                    )

                # builds the final allele representation and validates it just in case
                name = get_mhc2_isoform_name(alpha_chain, beta_chain)
                match = HLA_MOLECULE_PATTERN.match(name)
                match2 = HLA_DR_MOLECULE_PATTERN.match(name)
                assert (
                        match is not None or match2 is not None
                ), "Molecule does not match HLA isoform pattern {}".format(name)
            elif organism == ORGANISM_MUS_MUSCULUS:
                # TODO
                name = None
                alpha_chain = None
                beta_chain = None
            else:
                raise NeofoxDataValidationException("Not supported organism {}".format(organism))

        except AssertionError as e:
            logger.error(isoform.to_json(indent=3))
            raise NeofoxDataValidationException(e)

        return Mhc2Isoform(name=name, alpha_chain=alpha_chain, beta_chain=beta_chain)
