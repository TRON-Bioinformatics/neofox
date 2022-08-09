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
import json
import os
from abc import ABCMeta, abstractmethod, ABC
from typing import List

from logzero import logger
import neofox
from neofox.exceptions import NeofoxConfigurationException
import pandas as pd

from neofox.model.neoantigen import Mhc1Name, Mhc2GeneName, MhcAllele, Mhc2Name, Resource

DEFAULT_MAKEBLASTDB = 'makeblastdb'
DEFAULT_PRIME = 'PRIME'
DEFAULT_NETMHCPAN = 'netMHCpan'
DEFAULT_NETMHC2PAN = 'netMHCIIpan'
DEFAULT_RSCRIPT = 'Rscript'
DEFAULT_MIXMHCPRED = "MixMHCpred"
DEFAULT_MIXMHC2PRED = "MixMHC2pred_unix"
DEFAULT_BLASTP = "blastp"

ORGANISM_HOMO_SAPIENS = 'human'
HOMO_SAPIENS_MHC_I_GENES = [Mhc1Name.A, Mhc1Name.B, Mhc1Name.C]
HOMO_SAPIENS_MHC_II_GENES = [Mhc2GeneName.DPA1, Mhc2GeneName.DPB1, Mhc2GeneName.DQA1, Mhc2GeneName.DQB1,
                             Mhc2GeneName.DRB1]
HOMO_SAPIENS_MHC_II_MOLECULES = [Mhc2Name.DP, Mhc2Name.DQ, Mhc2Name.DR]
ORGANISM_MUS_MUSCULUS = 'mouse'
MUS_MUSCULUS_MHC_I_GENES = [Mhc1Name.H2K, Mhc1Name.H2L, Mhc1Name.H2D]
MUS_MUSCULUS_MHC_II_GENES = [Mhc2GeneName.H2E, Mhc2GeneName.H2A]
MUS_MUSCULUS_MHC_II_MOLECULES = [Mhc2Name.H2E_molecule, Mhc2Name.H2A_molecule]

MHC_I_GENES_BY_ORGANISM = {
    ORGANISM_HOMO_SAPIENS: HOMO_SAPIENS_MHC_I_GENES,
    ORGANISM_MUS_MUSCULUS: MUS_MUSCULUS_MHC_I_GENES
}
MHC_II_GENES_BY_ORGANISM = {
    ORGANISM_HOMO_SAPIENS: HOMO_SAPIENS_MHC_II_GENES,
    ORGANISM_MUS_MUSCULUS: MUS_MUSCULUS_MHC_II_GENES
}

PREFIX_HOMO_SAPIENS = "homo_sapiens"
HOMO_SAPIENS_FASTA = "Homo_sapiens.fa"
HOMO_SAPIENS_PICKLE = "Homo_sapiens.pickle"
PREFIX_MUS_MUSCULUS = "mus_musculus"
MUS_MUSCULUS_FASTA = "Mus_musculus.fa"
MUS_MUSCULUS_PICKLE = "Mus_musculus.pickle"
PROTEOME_DB_FOLDER = "proteome_db"

IEDB_FASTA_HOMO_SAPIENS = "IEDB_homo_sapiens.fasta"
IEDB_BLAST_PREFIX_HOMO_SAPIENS = "iedb_blast_db_homo_sapiens"
IEDB_FASTA_MUS_MUSCULUS = "IEDB_mus_musculus.fasta"
IEDB_BLAST_PREFIX_MUS_MUSCULUS = "iedb_blast_db_mus_musculus"
IEDB_FOLDER = "iedb"

NETMHCPAN_AVAILABLE_ALLELES_FILE = "netmhcpan_available_alleles_human.txt"
NETMHCPAN_AVAILABLE_ALLELES_MICE_FILE = "netmhcpan_available_alleles_mice.txt"
NETMHC2PAN_AVAILABLE_ALLELES_FILE = "netmhc2pan_available_alleles_human.txt"
NETMHC2PAN_AVAILABLE_ALLELES_MICE_FILE = "netmhc2pan_available_alleles_mice.txt"
HLA_DATABASE_AVAILABLE_ALLELES_FILE = "hla_database_allele_list.csv"
H2_DATABASE_AVAILABLE_ALLELES_FILE = "h2_database_allele_list.csv"
MIXMHCPRED_AVAILABLE_ALLELES_FILE = "allele_list.txt"
MIXMHC2PRED_AVAILABLE_ALLELES_FILE = "Alleles_list.txt"
PRIME_AVAILABLE_ALLELES_FILE = "alleles.txt"

RESOURCES_VERSIONS = "resources_versions.json"


class AbstractDependenciesConfiguration:
    def _check_and_load_binary(self, variable_name, default_value=None, optional=False, path_search=True):
        """
        Fetches the binary from the provided environment variable, if not available uses the default.
        If the binary is not a path to an executable, it searches for it in the PATH
        """

        program = None
        variable_value = os.environ.get(variable_name, default=default_value)

        if variable_value != '':
            fpath, _ = os.path.split(variable_value)
            if fpath:
                # makes sure that the provided path is absolute
                if not os.path.isabs(variable_value):
                    raise NeofoxConfigurationException(
                        "Please, use an absolute path in the environment variable ${}!".format(variable_name))
                # checks that it is executable
                if os.path.isfile(variable_value) and os.access(variable_value, os.X_OK):
                    program = variable_value
            elif path_search:
                # if no path searches for this command in the path
                for path in os.environ.get("PATH", "").split(os.pathsep):
                    exe_file = os.path.join(path, variable_value)
                    if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
                        program = variable_value
                        break

        if not optional and program is None:
            raise NeofoxConfigurationException(
                "Please, set the environment variable ${} pointing to the right binary and make sure to have "
                "execution permissions!".format(
                    variable_name
                )
            )

        return program


class DependenciesConfiguration(AbstractDependenciesConfiguration):
    def __init__(self):
        self.blastp = self._check_and_load_binary(neofox.NEOFOX_BLASTP_ENV, default_value=DEFAULT_BLASTP)
        self.mix_mhc2_pred = self._check_and_load_binary(
            neofox.NEOFOX_MIXMHC2PRED_ENV, default_value=DEFAULT_MIXMHC2PRED, optional=True, path_search=False)
        if self.mix_mhc2_pred is not None:
            self.mix_mhc2_pred_alleles_list = os.path.join(
                os.path.dirname(self.mix_mhc2_pred), MIXMHC2PRED_AVAILABLE_ALLELES_FILE)
        else:
            self.mix_mhc2_pred_alleles_list = None
        self.mix_mhc_pred = self._check_and_load_binary(
            neofox.NEOFOX_MIXMHCPRED_ENV, default_value=DEFAULT_MIXMHCPRED, optional=True, path_search=False)
        if self.mix_mhc_pred is not None:
            self.mix_mhc_pred_alleles_list = os.path.join(
                os.path.dirname(self.mix_mhc_pred), "lib", MIXMHCPRED_AVAILABLE_ALLELES_FILE)
        else:
            self.mix_mhc_pred_alleles_list = None
        self.rscript = self._check_and_load_binary(neofox.NEOFOX_RSCRIPT_ENV, default_value=DEFAULT_RSCRIPT)
        self.net_mhc2_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHC2PAN_ENV, default_value=DEFAULT_NETMHC2PAN)
        self.net_mhc_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHCPAN_ENV, default_value=DEFAULT_NETMHCPAN)
        self.prime = self._check_and_load_binary(neofox.NEOFOX_PRIME_ENV, default_value=DEFAULT_PRIME, optional=True,
                                                 path_search=False)
        if self.prime:
            self.prime_alleles_list = os.path.join(
                os.path.dirname(self.prime), "lib", PRIME_AVAILABLE_ALLELES_FILE
            )


class DependenciesConfigurationForInstaller(AbstractDependenciesConfiguration):
    def __init__(self):
        self.net_mhc2_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHC2PAN_ENV, default_value=DEFAULT_NETMHC2PAN)
        self.net_mhc_pan = self._check_and_load_binary(neofox.NEOFOX_NETMHCPAN_ENV, default_value=DEFAULT_NETMHCPAN)
        self.make_blastdb = self._check_and_load_binary(neofox.NEOFOX_MAKEBLASTDB_ENV, default_value=DEFAULT_MAKEBLASTDB)
        self.rscript = self._check_and_load_binary(neofox.NEOFOX_RSCRIPT_ENV, default_value=DEFAULT_RSCRIPT)


class MhcDatabase(ABC):

    organism = None
    mhc1_genes = None
    mhc2_genes = None
    mhc2_molecules = None

    def __init__(self, database_filename: str):
        super().__init__()
        self.alleles = self._load_alleles(database_filename)

    @abstractmethod
    def _load_alleles(self, hla_database_filename: str):
        pass

    @abstractmethod
    def exists(self, allele: MhcAllele):
        pass

    def is_homo_sapiens(self):
        return self.organism == ORGANISM_HOMO_SAPIENS

    def is_mus_musculus(self):
        return self.organism == ORGANISM_MUS_MUSCULUS


class HlaDatabase(MhcDatabase):

    organism = ORGANISM_HOMO_SAPIENS
    mhc1_genes = HOMO_SAPIENS_MHC_I_GENES
    mhc2_genes = HOMO_SAPIENS_MHC_II_GENES
    mhc2_molecules = HOMO_SAPIENS_MHC_II_MOLECULES

    def _load_alleles(self, hla_database_filename: str):
        # Assumes there is a column named Allele which contains the HLA alleles following the HLA nomenclature
        hla_database = pd.read_csv(hla_database_filename, comment="#")
        hla_database["parsed_allele"] = hla_database.Allele.transform(self._parse_allele)
        alleles = hla_database["parsed_allele"].dropna()
        alleles = alleles[
            alleles.str.startswith(Mhc1Name.A.name) | alleles.str.startswith(Mhc1Name.B.name) |
            alleles.str.startswith(Mhc1Name.C.name) |
            alleles.str.startswith(Mhc2GeneName.DPA1.name) | alleles.str.startswith(Mhc2GeneName.DPB1.name) |
            alleles.str.startswith(Mhc2GeneName.DQA1.name) | alleles.str.startswith(Mhc2GeneName.DQB1.name) |
            alleles.str.startswith(Mhc2GeneName.DRB1.name)]
        return alleles.unique()

    def _parse_allele(self, allele):
        values = allele.split(":")
        if len(values) > 1:
            return values[0] + ":" + values[1]
        else:
            return None

    def exists(self, allele: MhcAllele):
        return "{gene}*{group}:{protein}".format(
            gene=allele.gene, group=allele.group, protein=allele.protein) in self.alleles


class H2Database(MhcDatabase):

    organism = ORGANISM_MUS_MUSCULUS
    mhc1_genes = MUS_MUSCULUS_MHC_I_GENES
    mhc2_genes = MUS_MUSCULUS_MHC_II_GENES
    mhc2_molecules = MUS_MUSCULUS_MHC_II_MOLECULES

    def _load_alleles(self, h2_database_filename: str):
        # Assumes there is a column named Allele which contains the HLA alleles following the HLA nomenclature
        # the format of H2 alleles is gene and then allele without separator
        # (eg: H2Kb where H2K is the gene and b is the allele)
        h2_database = pd.read_csv(h2_database_filename, comment="#")
        alleles = h2_database["Allele"].dropna()
        alleles = alleles[
            alleles.str.startswith(Mhc1Name.H2D.name) |
            alleles.str.startswith(Mhc1Name.H2K.name) |
            alleles.str.startswith(Mhc1Name.H2L.name) |
            alleles.str.startswith(Mhc2GeneName.H2A.name) |
            alleles.str.startswith(Mhc2GeneName.H2E.name)]
        return alleles.unique()

    def exists(self, allele: MhcAllele):
        # because there is no nomenclature distinguishing alleles by their serological classification and then at
        # the sequence level we will here use the protein field only and skip the group
        return "{gene}{protein}".format(gene=allele.gene, protein=allele.protein) in self.alleles


class ReferenceFolder(object):

    def __init__(self, organism=ORGANISM_HOMO_SAPIENS, verbose=False):
        self.organism = organism
        if not organism == ORGANISM_HOMO_SAPIENS and not organism == ORGANISM_MUS_MUSCULUS:
            raise NeofoxConfigurationException(
                "Non supported organism: {}. Use {} or {}".format(
                    organism, ORGANISM_HOMO_SAPIENS, ORGANISM_MUS_MUSCULUS))
        self.reference_genome_folder = self._check_reference_genome_folder()

        self.iedb = self._get_reference_file_name(IEDB_FOLDER)
        self.proteome_db = self._get_reference_file_name(PROTEOME_DB_FOLDER)
        self.uniprot = self.get_proteome_fasta()
        self.uniprot_pickle = self.get_proteome_pickle()
        if self.organism == ORGANISM_HOMO_SAPIENS:
            self.mhc_database_filename = self._get_reference_file_name(HLA_DATABASE_AVAILABLE_ALLELES_FILE)
            self.available_mhc_ii = self._get_reference_file_name(NETMHC2PAN_AVAILABLE_ALLELES_FILE)
            self.available_mhc_i = self._get_reference_file_name(NETMHCPAN_AVAILABLE_ALLELES_FILE)
        elif self.organism == ORGANISM_MUS_MUSCULUS:
            self.mhc_database_filename = self._get_reference_file_name(H2_DATABASE_AVAILABLE_ALLELES_FILE)
            self.available_mhc_ii = self._get_reference_file_name(NETMHC2PAN_AVAILABLE_ALLELES_MICE_FILE)
            self.available_mhc_i = self._get_reference_file_name(NETMHCPAN_AVAILABLE_ALLELES_MICE_FILE)
        else:
            raise NeofoxConfigurationException("No support for organism {}".format(self.organism))

        self.resources_versions_file = self._get_reference_file_name(RESOURCES_VERSIONS)

        self.resources = [
            self.available_mhc_ii,
            self.available_mhc_i,
            self.iedb,
            self.proteome_db,
            self.uniprot,
            self.get_iedb_fasta(),
            self.mhc_database_filename,
            self.resources_versions_file
        ]
        self._check_resources()
        self.resources_versions = self.get_resources_versions()
        if verbose:
            self._log_configuration()
        self.__available_alleles = None
        self.__mhc_database = None

    def get_resources_versions(self):
        try:
            with open(self.resources_versions_file) as fd:
                resources_version = json.loads(fd.read())
        except FileNotFoundError:
            # NOTE: capturing this error is here to make unit tests easier, otherwise we need to create this resources
            # file always.
            resources_version = []
        return [Resource().from_dict(r) for r in resources_version]

    def get_available_alleles(self):
        # this enforces lazy initialisation (useful for testing)
        if not self.__available_alleles:
            self.__available_alleles = AvailableAlleles(self)
        return self.__available_alleles

    def get_mhc_database(self) -> MhcDatabase:
        # this enforces lazy initialisation (useful for testing)
        if not self.__mhc_database:
            if self.organism == ORGANISM_HOMO_SAPIENS:
                self.__mhc_database = HlaDatabase(self.mhc_database_filename)
            elif self.organism == ORGANISM_MUS_MUSCULUS:
                self.__mhc_database = H2Database(self.mhc_database_filename)
            else:
                raise NeofoxConfigurationException("No support for organism {}".format(self.organism))
        return self.__mhc_database

    def get_proteome_database(self):
        return os.path.join(
            self.proteome_db,
            PREFIX_HOMO_SAPIENS if self.organism == ORGANISM_HOMO_SAPIENS else PREFIX_MUS_MUSCULUS)

    def get_proteome_fasta(self):
        return os.path.join(
            self.proteome_db,
            HOMO_SAPIENS_FASTA if self.organism == ORGANISM_HOMO_SAPIENS
            else MUS_MUSCULUS_FASTA)

    def get_proteome_pickle(self):
        return os.path.join(
            self.proteome_db,
            HOMO_SAPIENS_PICKLE if self.organism == ORGANISM_HOMO_SAPIENS
            else MUS_MUSCULUS_PICKLE)

    def get_iedb_database(self):
        return os.path.join(
            self.iedb,
            IEDB_BLAST_PREFIX_HOMO_SAPIENS if self.organism == ORGANISM_HOMO_SAPIENS
            else IEDB_BLAST_PREFIX_MUS_MUSCULUS)

    def get_iedb_fasta(self):
        return os.path.join(
            self.iedb,
            IEDB_FASTA_HOMO_SAPIENS if self.organism == ORGANISM_HOMO_SAPIENS
            else IEDB_FASTA_MUS_MUSCULUS)

    def _check_reference_genome_folder(self):
        reference_genome_folder = os.environ.get(neofox.REFERENCE_FOLDER_ENV)
        if reference_genome_folder is None:
            raise NeofoxConfigurationException(
                "Please, set the environment variable ${} pointing to the reference genome folder!".format(
                    neofox.REFERENCE_FOLDER_ENV
                )
            )
        if not os.path.exists(reference_genome_folder):
            raise NeofoxConfigurationException(
                "The provided reference genome '{}' in ${} does not exist!".format(
                    reference_genome_folder, neofox.REFERENCE_FOLDER_ENV
                )
            )
        return reference_genome_folder

    def _check_resources(self):
        # check existence of all resources explicitly defined
        missing_resources = []
        for r in self.resources:
            if not os.path.exists(r):
                missing_resources.append(r)
        if len(missing_resources) > 0:
            raise NeofoxConfigurationException(
                "Missing resources in the reference folder: {}".format(
                    str(missing_resources)
                )
            )

    def _log_configuration(self):
        logger.info("Reference genome folder: {}".format(self.reference_genome_folder))
        logger.info("Resources:")
        for r in self.resources:
            logger.info(r)
        logger.info("Annotation resources:")
        for r in self.resources_versions:
            logger.info(r.to_json())

    def _get_reference_file_name(self, file_name_suffix):
        return os.path.join(self.reference_genome_folder, file_name_suffix)


class AvailableAlleles(object):
    def __init__(self, references):
        self.available_mhc_i = self._load_available_mhc_alleles(
            mhc=neofox.MHC_I, references=references
        )
        self.available_mhc_ii = self._load_available_mhc_alleles(
            mhc=neofox.MHC_II, references=references
        )

    def _load_available_mhc_alleles(self, references: ReferenceFolder, mhc=neofox.MHC_I) -> List:
        """
        loads file with available hla alllels for netmhcpan4/netmhcIIpan prediction, returns set
        """
        if mhc == neofox.MHC_II:
            fileMHC = references.available_mhc_ii
        else:
            fileMHC = references.available_mhc_i
        set_available_mhc = set()
        with open(fileMHC) as f:
            for line in f:
                set_available_mhc.add(line.strip())
        return set_available_mhc

    def get_available_mhc_i(self):
        return self.available_mhc_i

    def get_available_mhc_ii(self):
        return self.available_mhc_ii
