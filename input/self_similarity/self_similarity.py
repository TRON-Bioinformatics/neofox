#!/usr/bin/env python

import input.self_similarity.compute_self_similarity as compute_self_similarity
from input import MHC_I, MHC_II


def get_self_similarity(props, mhc, references):
    """Returns self-similiarity between mutated and wt epitope according to Bjerregard et al.,
    Argument mhc indicates if determination for MHC I or MHC II epitopes
    """
    if mhc == MHC_I:
        mhcI_mut = props["MHC_I_epitope_.best_prediction."]
        mhcI_wt = props["MHC_I_epitope_.WT."]
    elif mhc == MHC_II:
        mhcI_mut = props["MHC_II_epitope_.best_prediction."]
        mhcI_wt = props["MHC_II_epitope_.WT."]
    self_similarity = 'NA'
    try:
        self_similarity = str(compute_self_similarity.selfsim(references.blosum62).compute_K_hat_3(mhcI_mut, mhcI_wt))
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
        sc_mut = props["MHC_II_score_.best_prediction."].replace(",",".")
        sc_wt = props["MHC_II_score_.WT."].replace(",",".")

    try:
        improved_binder = float(sc_wt)/float(sc_mut) >= 1.2
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



def position_of_mutation_epitope(props, mhc):
    '''
    This function determines the position of the mutation within the epitope sequence.
    '''
    if mhc == MHC_I:
        mhc_mut = props["MHC_I_epitope_.best_prediction."]
        mhc_wt = props["MHC_I_epitope_.WT."]
    elif mhc == MHC_II:
        mhc_mut = props["MHC_II_epitope_.best_prediction."]
        mhc_wt = props["MHC_II_epitope_.WT."]
    p1 = -1
    try:
        for i,aa in enumerate(mhc_mut):
            if aa != mhc_wt[i]:
                p1 = i + 1
        return str(p1)
    except:
        return "NA"


def position_of_mutation_epitope_affinity(props, nine_mer = False):
    '''
    This function determines the position of the mutation within the epitope sequence.
    '''
    if nine_mer:
        mhc_mut = props["best_affinity_epitope_netmhcpan4_9mer"]
        mhc_wt = props["best_epitope_netmhcpan4_9mer_WT"]
    else:
        mhc_mut = props["best_affinity_epitope_netmhcpan4"]
        mhc_wt = props["best_affinity_epitope_netmhcpan4_WT"]
    p1 = -1
    try:
        for i,aa in enumerate(mhc_mut):
            if aa != mhc_wt[i]:
                p1 = i + 1
        return str(p1)
    except:
        return "NA"




def position_in_anchor_position(props, netMHCpan = False, nine_mer = False):
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
        anchor = int(pos_mhcI) == int(pep_len) or int(pos_mhcI)==2
        return str(1) if anchor else str(0)
    except:
        return "NA"
