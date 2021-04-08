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
from neofox.exceptions import NeofoxDataValidationException
from neofox.model.neoantigen import Mhc2GeneName, MhcAllele, Mhc2Isoform, Mhc1Name
import re
from logzero import logger

from neofox.model.wrappers import get_mhc2_isoform_name
from neofox.references.references import HlaDatabase

HLA_ALLELE_PATTERN_WITHOUT_SEPARATOR = re.compile(
    r"(?:HLA-)?((?:A|B|C|DPA1|DPB1|DQA1|DQB1|DRB1))[\*|_]?([0-9]{2,})([0-9]{2,3})[:|_]?([0-9]{2,})?[:|_]?([0-9]{2,})?([N|L|S|Q]{0,1})"
)
HLA_ALLELE_PATTERN = re.compile(
    r"(?:HLA-)?((?:A|B|C|DPA1|DPB1|DQA1|DQB1|DRB1))[\*|_]?([0-9]{2,})[:|_]?([0-9]{2,3})[:|_]?([0-9]{2,})?[:|_]?([0-9]{2,})?([N|L|S|Q]{0,1})"
)
HLA_MOLECULE_PATTERN = re.compile(
    r"(?:HLA-)?((?:DPA1|DPB1|DQA1|DQB1|DRB1)[\*|_]?[0-9]{2,}[:|_]?[0-9]{2,})[-|_]{1,2}"
    r"((?:DPA1|DPB1|DQA1|DQB1|DRB1)[\*|_]?[0-9]{2,}[:|_]?[0-9]{2,})"
)
HLA_DR_MOLECULE_PATTERN = re.compile(r"(?:HLA-)?(DRB1[\*|_]?[0-9]{2,}[:|_]?[0-9]{2,})")


class MhcParser:

    def __init__(self, hla_database: HlaDatabase):
        self.hla_database = hla_database

    def parse_mhc_allele(self, allele: str) -> MhcAllele:
        match = HLA_ALLELE_PATTERN_WITHOUT_SEPARATOR.match(allele)
        if match is not None:
            # allele without separator, controls for ambiguities
            gene = match.group(1)
            group = match.group(2)
            protein = match.group(3)
            default_allele_exists = self.hla_database.exists(gene, group, protein)
            if not default_allele_exists:
                # if default allele does not exist, tries alternative
                protein = group[-1:] + protein
                group = group[0: -1]
        else:
            # infers gene, group and protein from the name
            match = HLA_ALLELE_PATTERN.match(allele)
            assert match is not None, "Allele does not match HLA allele pattern {}".format(
                allele) if allele != "" else "Please check the format of provided alleles. An empty allele is provided"
            gene = match.group(1)
            group = match.group(2)
            protein = match.group(3)

        # controls for existence in the HLA database and warns the user
        if not self.hla_database.exists(gene, group, protein):
            logger.warning("Allele {} does not exist in the HLA database".format(allele))

        # builds a normalized representation of the allele
        name = "HLA-{gene}*{serotype}:{protein}".format(
            gene=gene, serotype=group, protein=protein
        )
        # ensures that full name stores the complete allele as provided but normalizes
        # its representation
        full_name = name
        six_digits_id = match.group(4)
        if six_digits_id is not None and six_digits_id != "":
            full_name = full_name + ":{}".format(six_digits_id)
            eight_digits_id = match.group(5)
            if eight_digits_id is not None and eight_digits_id != "":
                full_name = full_name + ":{}".format(eight_digits_id)
                expression_change = match.group(6)
                if expression_change is not None and expression_change != "":
                    full_name = full_name + expression_change
        return MhcAllele(
            full_name=full_name, name=name, gene=gene, group=group, protein=protein
        )

    def parse_mhc2_isoform(self, isoform: str) -> Mhc2Isoform:
        # TODO: this method currently fails for netmhc2pan alleles which are like 'HLA-DQA10509-DQB10630'
        # infers gene, group and protein from the name
        match = HLA_MOLECULE_PATTERN.match(isoform)
        if match:
            alpha_chain = self.parse_mhc_allele(match.group(1))
            beta_chain = self.parse_mhc_allele(match.group(2))
        else:
            match = HLA_DR_MOLECULE_PATTERN.match(isoform)
            assert (
                    match is not None
            ), "Molecule does not match HLA isoform pattern {}".format(isoform)
            alpha_chain = MhcAllele()
            beta_chain = self.parse_mhc_allele(match.group(1))
        # builds the final allele representation and validates it just in case
        name = get_mhc2_isoform_name(alpha_chain, beta_chain)
        return Mhc2Isoform(name=name, alpha_chain=alpha_chain, beta_chain=beta_chain)

    @staticmethod
    def validate_mhc_allele_representation(allele: MhcAllele) -> MhcAllele:
        try:
            full_name = None
            if allele.full_name:
                # infers gene, group and protein from the name
                match = HLA_ALLELE_PATTERN.match(allele.full_name)
                assert (
                        match is not None
                ), "Allele does not match HLA allele pattern {}".format(allele.name)
                gene = match.group(1)
                group = match.group(2)
                protein = match.group(3)
                full_name = allele.full_name
            elif allele.name:
                # infers gene, group and protein from the name
                match = HLA_ALLELE_PATTERN.match(allele.name)
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

            assert gene in list(Mhc1Name.__members__.keys()) + list(
                Mhc2GeneName.__members__.keys()
            ), "Gene not from classic MHC: {}".format(gene)
            # builds the final allele representation and validates it just in case
            name = "HLA-{gene}*{serotype}:{protein}".format(
                gene=gene, serotype=group, protein=protein
            )
            match = HLA_ALLELE_PATTERN.match(name)
            assert (
                    match is not None
            ), "Allele does not match HLA allele pattern {}".format(name)
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
    def validate_mhc2_isoform_representation(isoform: Mhc2Isoform) -> Mhc2Isoform:
        try:
            if isoform.name:
                # infers alpha and beta chains
                match = HLA_MOLECULE_PATTERN.match(isoform.name)
                if match:
                    alpha_chain = MhcParser.validate_mhc_allele_representation(
                        MhcAllele(name=match.group(1))
                    )
                    beta_chain = MhcParser.validate_mhc_allele_representation(
                        MhcAllele(name=match.group(2))
                    )
                else:
                    match = HLA_DR_MOLECULE_PATTERN.match(isoform.name)
                    assert (
                            match is not None
                    ), "Molecule does not match HLA isoform pattern {}".format(
                        isoform.name
                    )
                    alpha_chain = MhcAllele()
                    beta_chain = MhcParser.validate_mhc_allele_representation(
                        MhcAllele(name=match.group(1))
                    )
            elif isoform.alpha_chain and isoform.beta_chain:
                # infers name from gene, group and protein
                alpha_chain = MhcParser.validate_mhc_allele_representation(
                    isoform.alpha_chain
                )
                beta_chain = MhcParser.validate_mhc_allele_representation(
                    isoform.beta_chain
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
        except AssertionError as e:
            logger.error(isoform.to_json(indent=3))
            raise NeofoxDataValidationException(e)

        return Mhc2Isoform(name=name, alpha_chain=alpha_chain, beta_chain=beta_chain)
