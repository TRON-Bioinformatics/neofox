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
import pickle

import numpy as np
import scipy.io as sio
import os

from Bio.Data import IUPACData

from neofox.helpers.epitope_helper import EpitopeHelper
from neofox.model.neoantigen import PredictedEpitope

from functools import lru_cache

MAT = "SIRdata.mat"
GENES_EXPRESSION_PICKLE = "genes-expression.pickle"
ACIDS_FEATURES_PICKLE = "amino-acids-features.pickle"


class Preprocessor(object):
    def __init__(self):
        self.load_data = sio.loadmat(
            os.path.join(os.path.abspath(os.path.dirname(__file__)), MAT)
        )
        with open(
            os.path.join(
                os.path.abspath(os.path.dirname(__file__)), GENES_EXPRESSION_PICKLE
            ),
            "rb",
        ) as handle:
            self.dict_expression = pickle.load(handle)
        with open(
            os.path.join(
                os.path.abspath(os.path.dirname(__file__)), ACIDS_FEATURES_PICKLE
            ),
            "rb",
        ) as handle:
            self.dict_data = pickle.load(handle)

    @staticmethod
    def seq2bin(seq):
        aa = list(IUPACData.protein_letters_3to1.values())
        dict_aa = dict((i, j) for j, i in enumerate(aa))
        arr = np.zeros((1, 9 * 20))
        for ii, letter in enumerate(seq):
            arr[0, ii * 20 + dict_aa.get(letter)] = 1
        return arr

    @staticmethod
    def get_hydrophbicity(x, dict_):
        pair_letters = [c for c in x if c.isupper()]
        res = np.abs(dict_[pair_letters[0]] - dict_[pair_letters[1]])
        return res

    @staticmethod
    def get_size(x, dict_):
        pair_letters = [c for c in x if c.isupper()]
        res = np.abs(dict_[pair_letters[0]] - dict_[pair_letters[1]])
        return res

    @staticmethod
    def get_charge_change(x, dict_):
        pair_letters = [c for c in x if c.isupper()]
        if dict_[pair_letters[0]] == dict_[pair_letters[1]]:
            return 0
        else:
            return 1

    @staticmethod
    def get_charge_abs(x, dict_):
        pair_letters = [c for c in x if c.isupper()]
        if dict_[pair_letters[0]] == dict_[pair_letters[1]]:
            return 0
        else:
            return 1

    @staticmethod
    def get_polar(x, dict_):
        pair_letters = [c for c in x if c.isupper()]
        res = np.abs(dict_[pair_letters[0]] - dict_[pair_letters[1]])
        return res

    @staticmethod
    def get_absolute(x, dict_):
        pair_letters = [c for c in x if c.isupper()]
        res = np.abs(dict_[pair_letters[0]] - dict_[pair_letters[1]])
        return res

    @staticmethod
    def get_diffetenet(x, dict_):
        pair_letters = [c for c in x if c.isupper()]
        if dict_[pair_letters[0]] == dict_[pair_letters[1]]:
            return 0
        else:
            return 1

    def get_gene_expression(self, gene):
        res = self.dict_expression.get(gene, 0.0)
        return res

    def get_properties(self, amino_substitution):
        return np.asarray(
            [
                self.get_diffetenet(amino_substitution, self.dict_data["Charge"]),
                self.get_absolute(amino_substitution, self.dict_data["Size"]),
                self.get_absolute(amino_substitution, self.dict_data["Hydro"]),
                self.get_absolute(amino_substitution, self.dict_data["Charge"]),
                self.get_diffetenet(amino_substitution, self.dict_data["Polar"]),
            ]
        )

    def main(self, gene: str, epitope: PredictedEpitope):

        seq_arr = self.seq2bin(epitope.mutated_peptide)
        # tap score
        tap_mat = self.load_data.get("tap")
        tap_score = tap_mat.dot(seq_arr.T).ravel()
        # cleavge score
        clv_mat = self.load_data.get("clv")
        clv_mat = clv_mat[0, 20:200]
        clv_score = clv_mat.dot(seq_arr.T).ravel()

        position_of_mutation = EpitopeHelper.position_of_mutation_epitope(epitope=epitope)
        wild_type_aminoacid = epitope.wild_type_peptide[position_of_mutation - 1]  # it is 1-based
        mutated_aminoacid = epitope.mutated_peptide[position_of_mutation - 1]
        features_aa = self.get_properties(wild_type_aminoacid + mutated_aminoacid)
        # expresion
        expression_value = self.get_gene_expression(gene)

        result = np.hstack((expression_value, features_aa, clv_score, tap_score))

        mat_features = np.asarray([result])

        return mat_features
