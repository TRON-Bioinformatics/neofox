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
import base64
import hashlib

import betterproto

from neofox.model.neoantigen import Neoantigen, Mutation, Gene
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, IUPACData


class ModelValidator(object):

    @staticmethod
    def validate(model: betterproto.Message):
        # TODO: make this method capture appropriately validation issues whend dealing with int and float
        return model.__bytes__()

    # TODO: add patient validation: validate GTEx tissue and MHC alleles

    @staticmethod
    def validate_neoantigen(neoantigen: Neoantigen) -> Neoantigen:

        # checks format consistency first
        ModelValidator.validate(neoantigen)

        # checks gene
        # TODO: do we want to verify existence of gene and transcript id?
        ModelValidator._validate_gene(neoantigen.gene)

        # checks mutation
        neoantigen.mutation = ModelValidator._validate_mutation(neoantigen.mutation)

        # check the expression values
        ModelValidator._validate_expression_values(neoantigen)

        # infer other fields from the model
        return ModelValidator._enrich_neoantigen(neoantigen)

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
        assert len(mutation.left_flanking_region) > 0, "Empty left flanking region"
        for aa in mutation.left_flanking_region:
            ModelValidator._validate_aminoacid(aa)
        assert mutation.right_flanking_region is not None, "Empty right flanking region"
        assert len(mutation.right_flanking_region) > 0, "Empty right flanking region"
        for aa in mutation.right_flanking_region:
            ModelValidator._validate_aminoacid(aa)
        # checks the position
        assert mutation.position is not None, "Empty position"
        assert isinstance(mutation.position, int), "Position must be an integer"
        assert mutation.position > 0, "Position must be a 1-based positive integer"
        return mutation

    @staticmethod
    def _validate_gene(gene: Gene):
        assert gene.gene is not None and len(gene.gene) > 0, "Empty gene symbol"
        assert gene.transcript_identifier is not None and len(gene.transcript_identifier) > 0, \
            "Empty transcript identifier"
        assert gene.assembly is None or gene.assembly == "hg19", "Other reference genome than hg19 is not supported"

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
