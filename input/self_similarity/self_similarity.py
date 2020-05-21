#!/usr/bin/env python

from input import MHC_I, MHC_II

import math
import os

BETA = 0.11387
BLOSUM62_FILE_NAME = 'BLOSUM62-2.matrix.txt'


class SelfSimilarityCalculator():

    def __init__(self):
        blosum_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), BLOSUM62_FILE_NAME)
        blosum_dict = self._load_blosum(blosum_file)
        self.k1 = self._compute_k1(blosum_dict)

    def _compute_k1(self, blosum_dict):
        K1 = {}
        for i in list(blosum_dict.keys()):
            x = K1.get(i, {})
            for j in list(blosum_dict[i].keys()):
                x[j] = math.pow(blosum_dict[i][j], BETA)
            K1[i] = x
        return K1

    def _load_blosum(self, blosum):
        blosum_dict = {}
        colid = []
        rowid = []
        c = 0
        with open(blosum) as f:
            for line in f:
                c += 1
                if c == 1:
                    colid = line.strip("\n").split(" ")
                    continue
                w = line.strip("\n").split(" ")
                id = w[0]
                v = [float(x) for x in w[1:]]
                rowid.append(id)
                x = blosum_dict.get(id, {})
                for i, vi in enumerate(v):
                    x[colid[i]] = vi
                blosum_dict[id] = x
        return blosum_dict

    def compute_k_hat_3(self, x, y):  # K^3
        return self._compute_k3(x, y) / math.sqrt(self._compute_k3(x, x) * self._compute_k3(y, y))

    def _compute_k3(self, f, g):
        max_k = min(len(f), len(g))
        s = 0
        for k in range(1, max_k + 1):
            for i in range(len(f) - (k - 1)):
                u = f[i:i + k]
                for j in range(len(g) - (k - 1)):
                    v = g[j:j + k]
                    s += self._compute_k2k(u, v, self.k1)
        return s

    def _compute_k2k(self, u, v, K1):
        if len(u) != len(v):
            return None
        k = len(u)
        p = K1[u[0]][v[0]]
        for i in range(1, k):
            p = p * K1[u[i]][v[i]]
        return p


def get_self_similarity(mutation, wild_type):
    """Returns self-similiarity between mutated and wt epitope according to Bjerregard et al.,
    Argument mhc indicates if determination for MHC I or MHC II epitopes
    """
    self_similarity = 'NA'
    try:
        self_similarity = str(SelfSimilarityCalculator().compute_k_hat_3(mutation, wild_type))
    except ZeroDivisionError:
        pass
    return self_similarity


def is_improved_binder(props, mhc):
    '''
    This function checks if mutated epitope is improved binder according to Bjerregard et al.
    '''
    if mhc == MHC_I:
        sc_mut = props["best%Rank_netmhcpan4"]
        sc_wt = props["best%Rank_netmhcpan4_WT"]
    elif mhc == MHC_II:
        sc_mut = props["MHC_II_score_.best_prediction."].replace(",", ".")
        sc_wt = props["MHC_II_score_.WT."].replace(",", ".")

    try:
        improved_binder = float(sc_wt) / float(sc_mut) >= 1.2
    except (ZeroDivisionError, ValueError) as e:
        return "NA"
    return "1" if improved_binder else "0"


def selfsimilarity_of_conserved_binder_only(props):
    '''this function returns selfsimilarity for conserved binder but not for improved binder
    '''
    conserved_binder = props["ImprovedBinding_mhcI"]
    similiarity = props["Selfsimilarity_mhcI"]
    try:
        if conserved_binder == str(0):
            return similiarity
        else:
            return "NA"
    except (ZeroDivisionError, ValueError) as e:
        return "NA"


def position_of_mutation_epitope(wild_type, mutation):
    '''
    This function determines the position of the mutation within the epitope sequence.
    '''
    p1 = -1
    try:
        for i, aa in enumerate(mutation):
            if aa != wild_type[i]:
                p1 = i + 1
        return str(p1)
    except:
        return "NA"


def position_in_anchor_position(props, netMHCpan=False, nine_mer=False):
    '''
    This function determines if the mutation is located within an anchor position in mhc I.
    As an approximation, we assume that the second and the last position are anchor positions for all alleles.
    '''
    if netMHCpan:
        pos_mhcI = props["pos_MUT_MHCI_affinity_epi"]
        pep_len = len(props["best_epitope_netmhcpan4"])
    elif nine_mer:
        pos_mhcI = props["pos_MUT_MHCI_affinity_epi_9mer"]
        pep_len = 9
    else:
        pos_mhcI = props["pos_MUT_MHCI"]
        pep_len = props["MHC_I_peptide_length_.best_prediction."]

    anchor = 0
    try:
        anchor = int(pos_mhcI) == int(pep_len) or int(pos_mhcI) == 2
        return str(1) if anchor else str(0)
    except:
        return "NA"
