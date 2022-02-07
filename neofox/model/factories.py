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
from typing import List
import stringcase
from neofox import NOT_AVAILABLE_VALUE
from neofox.exceptions import NeofoxDataValidationException
from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.mhc_parser import MhcParser, get_mhc2_isoform_name
from neofox.model.neoantigen import Annotation, Patient, Mhc1, Zygosity, Mhc2, Mhc2Gene, Mhc2Name, Mhc2Isoform, \
    MhcAllele, Mhc2GeneName, Mhc1Name, Mutation, Neoantigen
from neofox.model.validation import ModelValidator, GENES_BY_MOLECULE
from neofox.references.references import MhcDatabase


class AnnotationFactory(object):
    @staticmethod
    def build_annotation(name, value):
        if isinstance(value, bool):
            # prints booleans as 0/1 strings
            value = "1" if value else "0"
        if not isinstance(value, str) and not isinstance(value, type(None)):
            if isinstance(value, float):
                value = "{0:.5g}".format(round(value, 5))
            else:
                value = str(value)
        if value is None:
            value = NOT_AVAILABLE_VALUE
        return Annotation(name=name, value=value)


class NeoantigenFactory(object):
    @staticmethod
    def build_neoantigen(wild_type_xmer=None, mutated_xmer=None, patient_identifier=None, gene=None,
                         rna_expression=None, rna_variant_allele_frequency=None, dna_variant_allele_frequency=None,
                         imputed_gene_expression=None, **kw):

        neoantigen = Neoantigen()
        neoantigen.patient_identifier = patient_identifier
        neoantigen.gene = gene
        neoantigen.rna_expression = rna_expression
        neoantigen.rna_variant_allele_frequency = rna_variant_allele_frequency
        neoantigen.dna_variant_allele_frequency = dna_variant_allele_frequency
        neoantigen.imputed_gene_expression = imputed_gene_expression

        mutation = Mutation()
        mutation.wild_type_xmer = wild_type_xmer
        mutation.mutated_xmer = mutated_xmer
        mutation.position = EpitopeHelper.mut_position_xmer_seq(mutation)
        neoantigen.mutation = mutation

        external_annotation_names = dict.fromkeys(
            nam for nam in kw.keys() if stringcase.snakecase(nam) not in set(Neoantigen.__annotations__.keys()))
        neoantigen.external_annotations = [
            Annotation(name=name, value=str(kw.get(name))) for name in external_annotation_names]

        ModelValidator.validate_neoantigen(neoantigen)

        return neoantigen


class PatientFactory(object):
    @staticmethod
    def build_patient(identifier, is_rna_available=False, tumor_type=None, mhc_alleles: List = [],
                      mhc2_alleles: List = [], mhc_database: MhcDatabase =None):
        patient = Patient(
            identifier=identifier,
            is_rna_available=is_rna_available,
            tumor_type=tumor_type,
            mhc1=MhcFactory.build_mhc1_alleles(mhc_alleles, mhc_database),
            mhc2=MhcFactory.build_mhc2_alleles(mhc2_alleles, mhc_database)
        )
        ModelValidator.validate_patient(patient=patient, organism=mhc_database.organism)
        return patient


class MhcFactory(object):

    @staticmethod
    def build_mhc1_alleles(alleles: List[str], mhc_database: MhcDatabase) -> List[Mhc1]:
        isoforms = []
        try:
            mhc_parser = MhcParser.get_mhc_parser(mhc_database)
            # NOTE: during the pandas parsing of empty columns empty lists become a list with one empty string
            parsed_alleles = list(map(mhc_parser.parse_mhc_allele, filter(lambda a: a != "", alleles)))
            for a in parsed_alleles:
                ModelValidator.validate_mhc1_gene(a)

            # do we need to validate genes anymore? add test creating MhcAllele with bad gene and see what happens
            for mhc1_gene in mhc_database.mhc1_genes:
                gene_alleles = list(
                    filter(lambda a: a.gene == mhc1_gene.name, parsed_alleles)
                )
                zygosity = MhcFactory._get_zygosity_from_alleles(gene_alleles)
                if zygosity == Zygosity.HOMOZYGOUS:
                    gene_alleles = [
                        gene_alleles[0]
                    ]  # we don't want repeated instances of the same allele
                isoforms.append(
                    Mhc1(name=mhc1_gene, zygosity=zygosity, alleles=gene_alleles)
                )
        except AssertionError as e:
            raise NeofoxDataValidationException(e)
        return list(filter(lambda i: i.zygosity != Zygosity.LOSS, isoforms))

    @staticmethod
    def build_mhc2_alleles(alleles: List[str], mhc_database: MhcDatabase) -> List[Mhc2]:
        mhc2s = []
        try:
            mhc_parser = MhcParser.get_mhc_parser(mhc_database)
            # NOTE: during the pandas parsing of empty columns empty lists become a list with one empty string
            parsed_alleles = list(map(mhc_parser.parse_mhc_allele, filter(lambda a: a != "", alleles)))
            for a in parsed_alleles:
                ModelValidator.validate_mhc2_gene(a)

            # do we need to validate genes anymore? add test creating MhcAllele with bad gene and see what happens
            for mhc2_isoform_name in mhc_database.mhc2_molecules:
                mhc2_isoform_genes = GENES_BY_MOLECULE.get(mhc2_isoform_name)
                isoform_alleles = list(
                    filter(lambda a: a.gene in [g.name for g in mhc2_isoform_genes], parsed_alleles)
                )
                genes = []
                for gene_name in mhc2_isoform_genes:
                    gene_alleles = list(
                        filter(lambda a: a.gene == gene_name.name, isoform_alleles)
                    )
                    zygosity = MhcFactory._get_zygosity_from_alleles(gene_alleles)
                    if zygosity == Zygosity.HOMOZYGOUS:
                        gene_alleles = [
                            gene_alleles[0]
                        ]  # we don't want repeated instances of the same allele
                    genes.append(
                        Mhc2Gene(
                            name=gene_name, zygosity=zygosity, alleles=gene_alleles
                        )
                    )
                isoforms = MhcFactory._get_mhc2_isoforms(mhc2_isoform_name, genes)
                mhc2s.append(Mhc2(name=mhc2_isoform_name, genes=genes, isoforms=isoforms))
        except AssertionError as e:
            raise NeofoxDataValidationException(e)
        return list(filter(lambda m: all(map(lambda g: g.zygosity != Zygosity.LOSS, m.genes)), mhc2s))

    @staticmethod
    def _get_mhc2_isoforms(
            isoform_name: Mhc2Name, genes: List[Mhc2Gene]
    ) -> List[Mhc2Isoform]:
        isoforms = []
        if isoform_name == Mhc2Name.DR:
            assert len(genes) <= 1, "More than one gene provided for MHC II DR"
            # alpha chain of the MHC II DR is not modelled as it is constant
            isoforms = [
                Mhc2Isoform(name=a.name, alpha_chain=MhcAllele(), beta_chain=a)
                for g in genes
                for a in g.alleles
            ]
        elif isoform_name == Mhc2Name.DP:
            assert len(genes) <= 2, "More than two genes provided for MHC II DP"
            alpha_alleles = [
                a for g in genes if g.name == Mhc2GeneName.DPA1 for a in g.alleles
            ]
            beta_alleles = [
                a for g in genes if g.name == Mhc2GeneName.DPB1 for a in g.alleles
            ]
            isoforms = [
                Mhc2Isoform(
                    name=get_mhc2_isoform_name(a, b), alpha_chain=a, beta_chain=b
                )
                for a in alpha_alleles
                for b in beta_alleles
            ]
        elif isoform_name == Mhc2Name.DQ:
            assert len(genes) <= 2, "More than two genes provided for MHC II DQ"
            alpha_alleles = [
                a for g in genes if g.name == Mhc2GeneName.DQA1 for a in g.alleles
            ]
            beta_alleles = [
                a for g in genes if g.name == Mhc2GeneName.DQB1 for a in g.alleles
            ]
            isoforms = [
                Mhc2Isoform(
                    name=get_mhc2_isoform_name(a, b), alpha_chain=a, beta_chain=b
                )
                for a in alpha_alleles
                for b in beta_alleles
            ]
        # mouse MHC II molecules do not act as pairs
        elif isoform_name == Mhc2Name.H2A_molecule:
            assert len(genes) <= 2, "More than two genes provided for H2A"
            isoforms = [
                Mhc2Isoform(name=a.name, alpha_chain=a, beta_chain=MhcAllele())
                for g in genes if g.name == Mhc2GeneName.H2A for a in g.alleles
            ]
        elif isoform_name == Mhc2Name.H2E_molecule:
            assert len(genes) <= 2, "More than two genes provided for H2E"
            isoforms = [
                Mhc2Isoform(name=a.name, alpha_chain=a, beta_chain=MhcAllele())
                for g in genes if g.name == Mhc2GeneName.H2E for a in g.alleles
            ]
        return isoforms

    @staticmethod
    def _get_zygosity_from_alleles(alleles: List[MhcAllele]) -> Zygosity:
        assert (
                len(set([a.gene for a in alleles])) <= 1
        ), "Trying to get zygosity from alleles of different genes"
        assert len(alleles) <= 2, "More than 2 alleles for gene {}".format(
            alleles[0].gene
        )
        if len(alleles) == 2 and alleles[0].name == alleles[1].name:
            zygosity = Zygosity.HOMOZYGOUS
        elif len(alleles) == 2:
            zygosity = Zygosity.HETEROZYGOUS
        elif len(alleles) == 1:
            zygosity = Zygosity.HEMIZYGOUS
        else:
            zygosity = Zygosity.LOSS
        return zygosity
