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
from abc import ABC, abstractmethod
from typing import List

from neofox.exceptions import NeofoxInputParametersException, NeofoxDataValidationException
from neofox.model.neoantigen import MhcAllele, Mhc2Isoform, Mhc2GeneName, Mhc2
import re
from logzero import logger

from neofox.references.references import MhcDatabase, ORGANISM_HOMO_SAPIENS, ORGANISM_MUS_MUSCULUS

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

H2_ALLELE_PATTERN = re.compile(r"(H2K|H2D|H2L|H2A|H2E)([a-z][0-9]?)")
H2_NETMHCPAN_ALLELE_PATTERN = re.compile(r"H-2-I?(K|D|L|A|E)([a-z][0-9]?)")
H2_MOLECULE_PATTERN = re.compile(r"(H2A|H2E)([a-z][0-9]?)")

ALLELE_PATTERN_BY_ORGANISM = {
    ORGANISM_HOMO_SAPIENS: HLA_ALLELE_PATTERN,
    ORGANISM_MUS_MUSCULUS: H2_ALLELE_PATTERN,
}


class MhcParser(ABC):

    def __init__(self, mhc_database: MhcDatabase):
        super().__init__()
        self.mhc_database = mhc_database

    @abstractmethod
    def parse_mhc_allele(self, allele: str):
        pass

    @abstractmethod
    def parse_mhc2_isoform(self, allele: str):
        pass

    @abstractmethod
    def get_netmhcpan_representation(self, allele: MhcAllele):
        pass

    @abstractmethod
    def get_netmhc2pan_representation(self, isoform: Mhc2Isoform):
        pass

    @staticmethod
    def get_mhc_parser(mhc_database: MhcDatabase):
        if mhc_database.is_homo_sapiens():
            mhc_parser = HlaParser(mhc_database=mhc_database)
        elif mhc_database.is_mus_musculus():
            mhc_parser = H2Parser(mhc_database=mhc_database)
        else:
            raise NeofoxInputParametersException("Organism not supported {}".format(mhc_database.organism))
        return mhc_parser


class H2Parser(MhcParser):

    def parse_mhc_allele(self, allele: str, pattern=H2_ALLELE_PATTERN) -> MhcAllele:
        match = H2_NETMHCPAN_ALLELE_PATTERN.match(allele)
        if match:
            # this ensures that netmhcpan output is normalized
            allele = "H2{gene}{protein}".format(gene=match.group(1), protein=match.group(2))
        match = H2_ALLELE_PATTERN.match(allele)
        if match is None:
            raise NeofoxDataValidationException(
                "Allele does not match H2 allele pattern {}".format(allele) if allele != "" else
                "Please check the format of provided alleles. An empty allele is provided")

        gene = match.group(1)
        protein = match.group(2)

        # controls for existence in the HLA database and warns the user
        mhc_allele = MhcAllele(gene=gene, protein=protein)
        if not self.mhc_database.exists(mhc_allele):
            logger.warning("Allele {} does not exist in the H2 database".format(allele))

        # builds a normalized representation of the allele
        name = "{gene}{protein}".format(gene=gene, protein=protein)

        # full name is the same as name in this case as the pattern does not allow variability
        mhc_allele.name = name
        mhc_allele.full_name = name
        return mhc_allele

    def parse_mhc2_isoform(self, allele: str) -> Mhc2Isoform:
        # MHC II molecules in H2 lab mouse are represented as single chain proteins
        # NOTE: by convention we represent this allele in both the alpha and beta chains
        match = H2_NETMHCPAN_ALLELE_PATTERN.match(allele)
        if match:
            # this ensures that netmhcpan output is normalized
            allele = "H2{gene}{protein}".format(gene=match.group(1), protein=match.group(2))
        allele = self.parse_mhc_allele(allele=allele, pattern=H2_MOLECULE_PATTERN)
        return Mhc2Isoform(name=allele.name, alpha_chain=allele, beta_chain=allele)

    def get_netmhcpan_representation(self, allele: MhcAllele):
        return "H-2-{gene}{protein}".format(gene=allele.gene.strip("H2"), protein=allele.protein)

    def get_netmhc2pan_representation(self, isoform: Mhc2Isoform):
        return "H-2-I{gene}{protein}".format(
            gene=isoform.alpha_chain.gene.strip("H2"), protein=isoform.alpha_chain.protein)


class HlaParser(MhcParser):

    def parse_mhc_allele(self, allele: str) -> MhcAllele:
        match = HLA_ALLELE_PATTERN_WITHOUT_SEPARATOR.match(allele)
        if match is not None:
            # allele without separator, controls for ambiguities
            gene = match.group(1)
            group = match.group(2)
            protein = match.group(3)
            default_allele_exists = self.mhc_database.exists(MhcAllele(gene=gene, group=group, protein=protein))
            if not default_allele_exists:
                # if default allele does not exist, tries alternative
                protein = group[-1:] + protein
                group = group[0: -1]
        else:
            # infers gene, group and protein from the name
            match = HLA_ALLELE_PATTERN.match(allele)
            if match is None:
                raise NeofoxDataValidationException(
                    "Allele does not match HLA allele pattern {}".format(allele) if allele != "" else
                    "Please check the format of provided alleles. An empty allele is provided")
            gene = match.group(1)
            group = match.group(2)
            protein = match.group(3)

        # controls for existence in the HLA database and warns the user
        mhc_allele = MhcAllele(gene=gene, group=group, protein=protein)
        if not self.mhc_database.exists(mhc_allele):
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
        mhc_allele.name = name
        mhc_allele.full_name = full_name
        return mhc_allele

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

    def get_netmhcpan_representation(self, allele: MhcAllele):
        return "HLA-{gene}{group}:{protein}".format(gene=allele.gene, group=allele.group, protein=allele.protein)

    def get_netmhc2pan_representation(self, isoform: Mhc2Isoform):
        if isoform.beta_chain.gene == Mhc2GeneName.DRB1.name:
            return "{gene}_{group}{protein}".format(
                gene=isoform.beta_chain.gene, group=isoform.beta_chain.group, protein=isoform.beta_chain.protein
            )
        else:
            return "HLA-{gene_a}{group_a}{protein_a}-{gene_b}{group_b}{protein_b}".format(
                gene_a=isoform.alpha_chain.gene,
                group_a=isoform.alpha_chain.group,
                protein_a=isoform.alpha_chain.protein,
                gene_b=isoform.beta_chain.gene,
                group_b=isoform.beta_chain.group,
                protein_b=isoform.beta_chain.protein,
            )


def get_alleles_by_gene(
    mhc_isoforms: List[Mhc2], gene: Mhc2GeneName
) -> List[MhcAllele]:
    return [
        a for m in mhc_isoforms for g in m.genes if g.name == gene for a in g.alleles
    ]


def get_mhc2_isoform_name(a: MhcAllele, b: MhcAllele):
    # NOTE: this is needed as jus setting alpha chain to None wouldn't work with protobuf
    if a is not None and a.name:
        return "{}-{}".format(a.name, b.name.replace("HLA-", ""))
    else:
        return b.name
