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
from neofox.model.neoantigen import Neoantigen, Mutation, Transcript, Patient, MhcAllele


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


def head(file_name, n=10):
    with open(file_name) as myfile:
        try:
            for x in range(n):
                print((next(myfile)))
        except StopIteration:
            pass


def print_and_delete(filename):
    head(filename)
    os.remove(filename)


def get_random_neoantigen():
    neoantigen = Neoantigen()
    neoantigen.variant_allele_frequency = np.random.uniform(0, 1)
    neoantigen.expression_value = np.random.uniform(0, 50)
    mutation = Mutation()
    mutation.mutated_aminoacid = random.choices(list(IUPACData.protein_letters), k=1)[0]
    mutation.wild_type_aminoacid = random.choices(list(IUPACData.protein_letters), k=1)[0]
    mutation.left_flanking_region = "".join(random.choices(list(IUPACData.protein_letters), k=5))
    mutation.right_flanking_region = "".join(random.choices(list(IUPACData.protein_letters), k=5))
    mutation.position = np.random.randint(0, 1000)
    neoantigen.mutation = mutation
    transcript = Transcript()
    transcript.gene = "BRCA2"
    transcript.identifier = "ENST1234567"
    transcript.assembly = "hg19"
    neoantigen.gene = transcript
    return neoantigen


def get_random_patient():
    patient = Patient()
    patient.estimated_tumor_content = np.random.uniform(0, 1)
    patient.is_rna_available = np.random.choice([True, False], 1)[0]
    patient.identifier = 'Pt12345'
    patient.mhc_i_alleles = [
        MhcAllele(name='HLA-A:01:01'), MhcAllele(name='HLA-B:01:01')]
    patient.mhc_i_i_alleles = [
        MhcAllele(name='HLA-DPA1:01:01'), MhcAllele(name='HLA-DPB1:01:01')]
    patient.tissue = 'skin'
    return patient