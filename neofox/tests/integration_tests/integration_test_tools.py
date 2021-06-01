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
import random
import tempfile

import dotenv
from Bio.Alphabet.IUPAC import IUPACData

from neofox.model.conversion import ModelConverter
from neofox.references.references import ReferenceFolder, DependenciesConfiguration


def load_references():
    dotenv.load_dotenv(override=True)
    return ReferenceFolder(), DependenciesConfiguration()


def create_temp_aminoacid_fasta_file():
    fastafile = tempfile.NamedTemporaryFile(mode="w", delete=False)
    with fastafile as f:
        f.write(get_random_kmer())
    return fastafile


def get_random_kmer(k=25):
    return "".join(random.choices(list(IUPACData.protein_letters), k=k))


def get_mhc_one_test(hla_database):
    return ModelConverter.parse_mhc1_alleles(
        [
            "HLA-A*24:02",
            "HLA-A*02:01",
            "HLA-B*15:01",
            "HLA-B*44:02",
            "HLA-C*07:02",
            "HLA-C*05:01",
        ], hla_database
    )


def get_mhc_two_test(hla_database):
    return ModelConverter.parse_mhc2_alleles(
            [
                "HLA-DRB1*04:02",
                "HLA-DRB1*08:01",
                "HLA-DQA1*03:01",
                "HLA-DQA1*04:01",
                "HLA-DQB1*03:02",
                "HLA-DQB1*04:02",
                "HLA-DPA1*01:03",
                "HLA-DPA1*02:01",
                "HLA-DPB1*13:01",
                "HLA-DPB1*04:01",
            ], hla_database
        )


mutations_with_rare_aminoacids = [
            ("UTTDSDGKF", "UTTDSWGKF"),  # this is an epitope from IEDB of length 9
            ("XTTDSDGKF", "XTTDSWGKF"),
            ("BTTDSDGKF", "BTTDSWGKF"),
            ("JTTDSDGKF", "JTTDSWGKF"),
            ("OTTDSDGKF", "OTTDSWGKF"),
            ("ZTTDSDGKF", "ZTTDSWGKF"),
            # only present in the wild type
            ("UTTDSDGKF", "TTTDSWGKF"),  # this is an epitope from IEDB of length 9
            ("XTTDSDGKF", "TTTDSWGKF"),
            ("BTTDSDGKF", "TTTDSWGKF"),
            ("JTTDSDGKF", "TTTDSWGKF"),
            ("OTTDSDGKF", "TTTDSWGKF"),
            ("ZTTDSDGKF", "TTTDSWGKF"),
            # only present in the mutation
            ("TTTDSDGKF", "UTTDSWGKF"),  # this is an epitope from IEDB of length 9
            ("TTTDSDGKF", "XTTDSWGKF"),
            ("TTTDSDGKF", "BTTDSWGKF"),
            ("TTTDSDGKF", "JTTDSWGKF"),
            ("TTTDSDGKF", "OTTDSWGKF"),
            ("TTTDSDGKF", "ZTTDSWGKF")
        ]
