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
import os
import random
import numpy as np
from Bio.Data import IUPACData
from mock import Mock
from neofox.model.neoantigen import Neoantigen, Mutation, Patient, MhcAllele


def mock_file_existence(existing_files=[], non_existing_files=[]):
    original_os_path_exists = os.path.exists

    def side_effect(filename):
        if filename in existing_files:
            return True
        elif filename in non_existing_files:
            return False
        else:
            return original_os_path_exists(filename)

    os.path.exists = Mock(side_effect=side_effect)


def mock_file_is_executable(executable_files=[], non_executable_files=[]):
    original_os_access = os.access

    def side_effect(filename, mode):
        if filename in executable_files and mode == os.X_OK:
            return True
        elif filename in non_executable_files and mode == os.X_OK:
            return False
        else:
            return original_os_access(filename, mode)

    os.access = Mock(side_effect=side_effect)


def get_random_neoantigen():
    neoantigen = Neoantigen()
    neoantigen.variant_allele_frequency = np.random.uniform(0, 1)
    neoantigen.expression_value = np.random.uniform(0, 50)
    mutation = Mutation()
    mutated_aminoacid = random.choices(list(IUPACData.protein_letters), k=1)[0]
    wild_type_aminoacid = random.choices(list(IUPACData.protein_letters), k=1)[0]
    left_flanking_region = "".join(random.choices(list(IUPACData.protein_letters), k=5))
    right_flanking_region = "".join(
        random.choices(list(IUPACData.protein_letters), k=5)
    )
    mutation.mutated_xmer = (
        left_flanking_region + mutated_aminoacid + right_flanking_region
    )
    mutation.wild_type_xmer = (
        left_flanking_region + wild_type_aminoacid + right_flanking_region
    )
    neoantigen.mutation = mutation
    neoantigen.gene = "BRCA2"
    return neoantigen
