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
import re
import base64
import hashlib
from typing import List, Tuple

import betterproto
from logzero import logger
from neofox.references.references import AvailableAlleles

from neofox.exceptions import NeofoxDataValidationException
from neofox.model.neoantigen import Neoantigen, Mutation, Gene, Patient, HlaAllele
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, IUPACData

DQB1 = "DQB1"
DQA1 = "DQA1"
DPB1 = "DPB1"
DPA1 = "DPA1"
DRB1 = "DRB1"


class ModelValidator(object):

    HLA_ALLELE_PATTERN = re.compile(
        r"(?:HLA-)(\w+)\*?([0-9]{2}):?([0-9]{2,}):?([0-9]{2,})?:?([0-9]{2,})?([N|L|S|Q]{0,1})")
    VALID_MHC_I_GENES = ["A", "B", "C"]     # only MHC I classical
    VALID_MHC_II_GENES = [DRB1, DPA1, DPB1, DQA1, DQB1]

    @staticmethod
    def validate(model: betterproto.Message):
        # TODO: make this method capture appropriately validation issues whend dealing with int and float
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
            # checks gene
            # TODO: do we want to verify existence of gene and transcript id?
            neoantigen.gene = ModelValidator._validate_gene(neoantigen.gene)

            # checks mutation
            neoantigen.mutation = ModelValidator._validate_mutation(neoantigen.mutation)

            # check the expression values
            ModelValidator._validate_expression_values(neoantigen)
        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        # infer other fields from the model
        return ModelValidator._enrich_neoantigen(neoantigen)

    @staticmethod
    def validate_patient(patient: Patient) -> Patient:

        # checks format consistency first
        ModelValidator.validate(patient)

        try:
            # checks that patient id is not empty considering white spaces
            patient_id = patient.identifier.strip() if patient.identifier else patient.identifier
            assert patient_id is not None and patient_id != "", "Patient identifier is empty"
            patient.identifier = patient_id

            # checks MHC I alleles
            patient.mhc_i_alleles = ModelValidator._validate_mhc_alleles(
                patient.mhc_i_alleles, ModelValidator.VALID_MHC_I_GENES)

            # checks MHC II alleles
            patient.mhc_i_i_alleles = ModelValidator._validate_mhc_alleles(
                patient.mhc_i_i_alleles, ModelValidator.VALID_MHC_II_GENES)

        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        return patient

    @staticmethod
    def _validate_mhc_alleles(alleles: List[HlaAllele], valid_genes: List[str]):

        # parses the individual alleles and perform some individual validation
        parsed_alleles = [ModelValidator._validate_mhc_allele_representation(
            a, valid_genes=valid_genes) for a in alleles]

        # checks the genotypes for MHC I genes
        for g in valid_genes:
            assert len(list(filter(lambda x: x.gene == g, parsed_alleles))) <= 2, \
                "MHC I gene {} has more than 2 alleles".format(g)
        return parsed_alleles

    @staticmethod
    def _validate_mhc_allele_representation(allele: HlaAllele, valid_genes: List[str]) -> HlaAllele:
        if allele.name:
            # infers gene, group and protein from the name
            match = ModelValidator.HLA_ALLELE_PATTERN.match(allele.name)
            assert match is not None, "Allele does not match HLA allele pattern {}".format(allele.name)
            gene = match.group(1)
            assert gene in valid_genes, "HLA gene not valid {}".format(gene)
            group = match.group(2)
            protein = match.group(3)
        elif allele.gene and allele.group and allele.protein:
            # infers name from gene, group and protein
            assert allele.gene in valid_genes, "HLA gene not valid {}".format(allele.gene)
            gene = allele.gene
            group = allele.group
            protein = allele.protein
        else:
            raise NeofoxDataValidationException("HLA allele missing required fields, either name or gene, group and "
                                                "protein must be provided")

        # builds the final allele representation and validates it just in case
        name = "HLA-{gene}*{serotype}:{protein}".format(gene=gene, serotype=group, protein=protein)
        match = ModelValidator.HLA_ALLELE_PATTERN.match(name)
        assert match is not None, "Allele does not match HLA allele pattern {}".format(name)

        return HlaAllele(name=name, gene=gene, group=group, protein=protein)

    @staticmethod
    def _validate_expression_values(neoantigen):
        assert neoantigen.rna_expression is None or neoantigen.rna_expression >= 0, \
            "RNA expression should be a positive integer or zero {}".format(neoantigen.rna_expression)
        ModelValidator._validate_vaf(neoantigen.dna_variant_allele_frequency)
        ModelValidator._validate_vaf(neoantigen.rna_variant_allele_frequency)

    @staticmethod
    def _validate_mutation(mutation: Mutation) -> Mutation:
        # checks aminoacids
        mutation.mutated_aminoacid = ModelValidator._validate_aminoacid(mutation.mutated_aminoacid)
        mutation.wild_type_aminoacid = ModelValidator._validate_aminoacid(mutation.wild_type_aminoacid)

        # checks left and right flanking regions
        assert mutation.left_flanking_region is not None, "Empty left flanking region"
        mutation.left_flanking_region = mutation.left_flanking_region.strip()
        assert len(mutation.left_flanking_region) > 0, "Empty left flanking region"
        for aa in mutation.left_flanking_region:
            ModelValidator._validate_aminoacid(aa)

        assert mutation.right_flanking_region is not None, "Empty right flanking region"
        mutation.right_flanking_region = mutation.right_flanking_region.strip()
        assert len(mutation.right_flanking_region) > 0, "Empty right flanking region"
        for aa in mutation.right_flanking_region:
            ModelValidator._validate_aminoacid(aa)

        # checks the position
        assert mutation.position is not None, "Empty position"
        assert isinstance(mutation.position, int), "Position must be an integer"
        assert mutation.position > 0, "Position must be a 1-based positive integer"
        return mutation

    @staticmethod
    def _validate_gene(gene: Gene) -> Gene:

        # TODO: validate that gene symbol exists
        gene_name = gene.gene.strip() if gene.gene else gene.gene
        assert gene_name is not None and len(gene_name) > 0, "Empty gene symbol"
        gene.gene = gene_name

        # TODO: validate that transcript identifier exists
        transcript_identifier = gene.transcript_identifier.strip() if gene.transcript_identifier \
            else gene.transcript_identifier
        assert transcript_identifier is not None and len(transcript_identifier) > 0, \
            "Empty transcript identifier"
        gene.transcript_identifier = transcript_identifier

        # TODO: support other assemblies
        assembly = gene.assembly if gene.assembly else "hg19"
        assert assembly == "hg19", "Other reference genome than hg19 is not supported"
        gene.assembly = assembly

        return gene

    @staticmethod
    def _validate_vaf(vaf):
        assert vaf is None or 0.0 <= vaf <= 1.0, "VAF should be a positive integer or zero {}".format(vaf)

    @staticmethod
    def _enrich_neoantigen(neoantigen: Neoantigen) -> Neoantigen:
        neoantigen.mutation.wild_type_xmer = "".join([
            neoantigen.mutation.left_flanking_region,
            neoantigen.mutation.wild_type_aminoacid,
            neoantigen.mutation.right_flanking_region])
        neoantigen.mutation.mutated_xmer = "".join([
            neoantigen.mutation.left_flanking_region,
            neoantigen.mutation.mutated_aminoacid,
            neoantigen.mutation.right_flanking_region])
        neoantigen.mutation.size_left_flanking_region = len(neoantigen.mutation.left_flanking_region)
        neoantigen.mutation.size_right_flanking_region = len(neoantigen.mutation.right_flanking_region)
        neoantigen.identifier = ModelValidator.generate_neoantigen_identifier(neoantigen)
        return neoantigen

    @staticmethod
    def _validate_aminoacid(aminoacid):
        assert aminoacid is not None, "Aminoacid field cannot be empty"
        aminoacid = aminoacid.strip()
        assert isinstance(aminoacid, str), "Aminoacid has to be a string"
        if len(aminoacid) == 3:
            assert aminoacid in IUPACData.protein_letters_3to1_extended.keys(), \
                "Non existing 3 letter aminoacid {}".format(aminoacid)
            aminoacid = IUPACData.protein_letters_3to1_extended.get(aminoacid)
        if len(aminoacid) == 1:
            assert aminoacid in ExtendedIUPACProtein.letters, "Non existing aminoacid {}".format(aminoacid)
        else:
            assert False, "Invalid aminoacid {}".format(aminoacid)
        return aminoacid

    @staticmethod
    def generate_neoantigen_identifier(neoantigen: Neoantigen) -> str:
        neoantigen.identifier = None    # this needs to be done otherwise we cannot rreproduce the id after it is set
        return base64.b64encode(hashlib.md5(neoantigen.to_json().encode('utf8')).digest()).decode('utf8')
