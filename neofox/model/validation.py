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
from neofox.model.conversion import ModelConverter
from neofox.exceptions import NeofoxDataValidationException
from neofox.model.neoantigen import Neoantigen, Mutation, Transcript, Patient, MhcAllele, Mhc2Isoform, Mhc1, \
    Mhc1Name, Zygosity, Mhc2, Mhc2Name, Mhc2GeneName
from Bio.Alphabet.IUPAC import ExtendedIUPACProtein, IUPACData


class ModelValidator(object):

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
            neoantigen.transcript = ModelValidator._validate_transcript(neoantigen.transcript)

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

            # TODO: validate new model with isoforms, genes and alleles
            # checks MHC I
            validated_mhc1s = []
            for m in patient.mhc1:
                validated_mhc1s.append(ModelValidator._validate_mhc1(m))
            patient.mhc1 = validated_mhc1s
            # checks MHC II
            validated_mhc2s = []
            for m in patient.mhc2:
                validated_mhc2s.append(ModelValidator._validate_mhc2(m))
            patient.mhc2 = validated_mhc2s

        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        return patient

    @staticmethod
    def _validate_mhc1(mhc1: Mhc1) -> Mhc1:
        assert mhc1.name in Mhc1Name, "Invalid MHC I name"
        assert mhc1.zygosity in Zygosity, "Invalid zygosity"
        alleles = mhc1.alleles
        if mhc1.zygosity in [Zygosity.HOMOZYGOUS, Zygosity.HEMIZYGOUS]:
            assert len(alleles) == 1, "An homozygous gene must have 1 allele and not {}".format(len(alleles))
        elif mhc1.zygosity == Zygosity.HETEROZYGOUS:
            assert len(alleles) == 2, "An heterozygous or hemizygous gene must have 2 alleles and not {}".format(
                len(alleles))
        elif mhc1.zygosity == Zygosity.LOSS:
            assert len(alleles) == 0, "An lost gene must have 0 alleles and not {}".format(
                len(alleles))
        validated_alleles = []
        for allele in alleles:
            validated_allele = ModelValidator.validate_mhc_allele_representation(allele)
            validated_alleles.append(validated_allele)
            assert validated_allele.gene == mhc1.name.name, \
                "The allele referring to gene {} is inside gene {}".format(validated_allele.gene, mhc1.name.name)
        mhc1.alleles = validated_alleles
        return mhc1

    @staticmethod
    def _validate_mhc2(mhc2: Mhc2) -> Mhc2:
        assert mhc2.name in Mhc2Name, "Invalid MHC II name"
        genes = mhc2.genes
        for gene in genes:
            assert gene.name in Mhc2GeneName, "Invalid gene name from MHC II"
            assert gene.name in ModelConverter.GENES_BY_MOLECULE.get(mhc2.name), \
                "Gene {} referring to isoform {}".format(gene.name, mhc2.name.name)
            assert gene.zygosity in Zygosity, "Invalid zygosity"
            alleles = gene.alleles
            if gene.zygosity == Zygosity.HOMOZYGOUS:
                assert len(alleles) == 1, "An homozygous gene must have 1 allele and not {}".format(len(alleles))
            elif gene.zygosity in [Zygosity.HETEROZYGOUS, Zygosity.HEMIZYGOUS]:
                assert len(alleles) == 2, "An heterozygous or hemizygous gene must have 2 alleles and not {}".format(
                    len(alleles))
            elif gene.zygosity == Zygosity.LOSS:
                assert len(alleles) == 0, "An lost gene must have 0 alleles and not {}".format(
                    len(alleles))
            validated_alleles = []
            for allele in alleles:
                validated_allele = ModelValidator.validate_mhc_allele_representation(allele)
                validated_alleles.append(validated_allele)
                assert validated_allele.gene == gene.name.name, \
                    "The allele referring to gene {} is inside gene {}".format(validated_allele.gene, gene.name.name)
            gene.alleles = validated_alleles
        isoforms = mhc2.isoforms
        validated_isoforms = []
        for isoform in isoforms:
            validated_isoform = ModelValidator.validate_mhc2_isoform_representation(isoform)
            validated_isoforms.append(validated_isoform)
            if mhc2.name != Mhc2Name.DR:
                assert validated_isoform.alpha_chain.name in [a.name for g in genes for a in g.alleles], \
                    "Alpha chain allele not present in th list of alleles"
            assert validated_isoform.beta_chain.name in [a.name for g in genes for a in g.alleles], \
                "Beta chain allele not present in th list of alleles"
        mhc2.isoforms = validated_isoforms
        return mhc2

    @staticmethod
    def validate_mhc_allele_representation(allele: MhcAllele) -> MhcAllele:
        try:
            full_name = None
            if allele.full_name:
                # infers gene, group and protein from the name
                match = ModelConverter.HLA_ALLELE_PATTERN.match(allele.full_name)
                assert match is not None, "Allele does not match HLA allele pattern {}".format(allele.name)
                gene = match.group(1)
                group = match.group(2)
                protein = match.group(3)
                full_name = allele.full_name
            elif allele.name:
                # infers gene, group and protein from the name
                match = ModelConverter.HLA_ALLELE_PATTERN.match(allele.name)
                assert match is not None, "Allele does not match HLA allele pattern {}".format(allele.name)
                gene = match.group(1)
                group = match.group(2)
                protein = match.group(3)
            elif allele.gene and allele.group and allele.protein:
                # infers name from gene, group and protein
                gene = allele.gene
                group = allele.group
                protein = allele.protein
            else:
                raise NeofoxDataValidationException(
                    "HLA allele missing required fields, either name or gene, group and protein must be provided")

            assert gene in list(Mhc1Name.__members__.keys()) + list(Mhc2GeneName.__members__.keys()), \
                "Gene not from classic MHC: {}".format(gene)
            # builds the final allele representation and validates it just in case
            name = "HLA-{gene}*{serotype}:{protein}".format(gene=gene, serotype=group, protein=protein)
            match = ModelConverter.HLA_ALLELE_PATTERN.match(name)
            assert match is not None, "Allele does not match HLA allele pattern {}".format(name)
        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        return MhcAllele(full_name=full_name if full_name else name, name=name, gene=gene, group=group, protein=protein)

    @staticmethod
    def validate_mhc2_isoform_representation(isoform: Mhc2Isoform) -> Mhc2Isoform:
        try:
            if isoform.name:
                # infers alpha and beta chains
                match = ModelConverter.HLA_MOLECULE_PATTERN.match(isoform.name)
                if match:
                    alpha_chain = ModelValidator.validate_mhc_allele_representation(MhcAllele(name=match.group(1)))
                    beta_chain = ModelValidator.validate_mhc_allele_representation(MhcAllele(name=match.group(2)))
                else:
                    match = ModelConverter.HLA_DR_MOLECULE_PATTERN.match(isoform.name)
                    assert match is not None, "Molecule does not match HLA isoform pattern {}".format(isoform.name)
                    alpha_chain = MhcAllele()
                    beta_chain = ModelValidator.validate_mhc_allele_representation(MhcAllele(name=match.group(1)))
            elif isoform.alpha_chain and isoform.beta_chain:
                # infers name from gene, group and protein
                alpha_chain = ModelValidator.validate_mhc_allele_representation(isoform.alpha_chain)
                beta_chain = ModelValidator.validate_mhc_allele_representation(isoform.beta_chain)
            else:
                raise NeofoxDataValidationException("HLA isoform missing required fields")

            # builds the final allele representation and validates it just in case
            name = ModelConverter.get_mhc2_isoform_name(alpha_chain, beta_chain)
            match = ModelConverter.HLA_MOLECULE_PATTERN.match(name)
            match2 = ModelConverter.HLA_DR_MOLECULE_PATTERN.match(name)
            assert match is not None or match2 is not None, \
                "Molecule does not match HLA isoform pattern {}".format(name)
        except AssertionError as e:
            raise NeofoxDataValidationException(e)

        return Mhc2Isoform(name=name, alpha_chain=alpha_chain, beta_chain=beta_chain)

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
    def _validate_transcript(transcript: Transcript) -> Transcript:

        # TODO: validate that gene symbol exists
        gene_name = transcript.gene.strip() if transcript.gene else transcript.gene
        assert gene_name is not None and len(gene_name) > 0, "Empty gene symbol"
        transcript.gene = gene_name

        # TODO: validate that transcript identifier exists
        transcript_identifier = transcript.identifier.strip() if transcript.identifier else transcript.identifier
        assert transcript_identifier is not None and len(transcript_identifier) > 0, "Empty transcript identifier"
        transcript.identifier = transcript_identifier

        # TODO: support other assemblies
        assembly = transcript.assembly if transcript.assembly else "hg19"
        assert assembly == "hg19", "Other reference genome than hg19 is not supported"
        transcript.assembly = assembly

        return transcript

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
