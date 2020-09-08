import os

# TODO: change this import when we move to python3
# from unittest.mock import Mock
import random

import numpy as np
from Bio.Data import IUPACData
from mock import Mock

from neofox.model.neoantigen import Neoantigen, Mutation, Gene, Patient


def _mock_file_existence(existing_files=[], unexisting_files=[]):
    original_os_path_exists = os.path.exists

    def side_effect(filename):
        if filename in existing_files:
            return True
        elif filename in unexisting_files:
            return False
        else:
            return original_os_path_exists(filename)

    os.path.exists = Mock(side_effect=side_effect)


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
    mutation.position = np.random.randint(0, 1000)
    neoantigen.mutation = mutation
    gene = Gene()
    gene.gene = "BRCA2"
    gene.transcript_identifier = "ENST1234567"
    gene.assembly = "hg19"
    neoantigen.gene = gene
    return neoantigen


def get_random_patient():
    patient = Patient()
    patient.estimated_tumor_content = np.random.uniform(0, 1)
    patient.is_rna_available = np.random.choice([True, False], 1)[0]
    patient.identifier = 'Pt12345'
    patient.mhc_i_alleles = ['A', 'B', 'C']
    patient.mhc_i_i_alleles = ['X', 'Y']
    patient.tissue = 'skin'
    return patient