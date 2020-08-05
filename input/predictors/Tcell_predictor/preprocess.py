import pickle

import numpy as np
import scipy.io as sio
import os


MAT = "SIRdata.mat"
GENES_EXPRESSION_PICKLE = "genes-expression.pickle"
ACIDS_FEATURES_PICKLE = "amino-acids-features.pickle"


class Preprocessor(object):

    def __init__(self):
        self.load_data = sio.loadmat(os.path.join(os.path.abspath(os.path.dirname(__file__)), MAT))
        with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), GENES_EXPRESSION_PICKLE), 'rb') as handle:
            self.dict_expression = pickle.load(handle)
        with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), ACIDS_FEATURES_PICKLE), 'rb') as handle:
            self.dict_data = pickle.load(handle)

    @staticmethod
    def seq2bin(seq):
        aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
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
        return np.asarray([self.get_diffetenet(amino_substitution, self.dict_data['Charge']),
                           self.get_absolute(amino_substitution, self.dict_data['Size']),
                           self.get_absolute(amino_substitution, self.dict_data['Hydro']),
                           self.get_absolute(amino_substitution, self.dict_data['Charge']),
                           self.get_diffetenet(amino_substitution, self.dict_data['Polar'])])

    def main(self, f_name):
        lst_data = []
        with open(f_name, 'r') as f:
            for row in f:
                gene_name, sequence, aa_subs = row.split()
                seq_arr = self.seq2bin(sequence)
                # tap score
                tap_mat = self.load_data.get('tap')
                tap_score = tap_mat.dot(seq_arr.T).ravel()
                # cleavge score
                clv_mat = self.load_data.get('clv')
                clv_mat = clv_mat[0, 20:200]
                clv_score = clv_mat.dot(seq_arr.T).ravel()

                features_aa = self.get_properties(aa_subs)
                # expresion
                expression_value = self.get_gene_expression(gene_name)

                lst_data.append(np.hstack((expression_value, features_aa, clv_score, tap_score)))
                mat_features = np.asarray(lst_data)
        return mat_features
