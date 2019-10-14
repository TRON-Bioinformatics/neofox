#!/usr/bin/env python

import sys
import os.path
import compute_self_similarity

my_path = os.path.abspath(os.path.dirname(__file__))

def selfsimilarity(props, mhc):
    """Returns self-similiarity between mutated and wt epitope according to Bjerregard et al.,
    Argument mhc indicates if determination for MHC I or MHC II epitopes
    """
    if mhc == "mhcI":
        mhcI_mut = props["MHC_I_epitope_.best_prediction."]
        mhcI_wt = props["MHC_I_epitope_.WT."]
    elif mhc == "mhcII":
        mhcI_mut = props["MHC_II_epitope_.best_prediction."]
        mhcI_wt = props["MHC_II_epitope_.WT."]
    selfsim = compute_self_similarity.selfsim(os.path.join(my_path, "./BLOSUM62-2.matrix.txt"))
    #print mhcI_mut, mhcI_wt
    try:
        return str(selfsim.compute_K_hat_3(mhcI_mut, mhcI_wt))
    except ZeroDivisionError:
        return "NA"


def improved_binder(props, mhc):
    '''
    This function checks if mutated epitope is improved binder according to Bjerregard et al.
    '''
    if mhc == "mhcI":
        #sc_mut = props["MHC_I_score_.best_prediction."].replace(",",".")
        #sc_wt = props["MHC_I_score_.WT."].replace(",",".")
        score_mut = props["best%Rank_netmhcpan4"]
        score_wt = props["best%Rank_netmhcpan4_WT"]
    elif mhc == "mhcII":
        sc_mut = props["MHC_II_score_.best_prediction."].replace(",",".")
        sc_wt = props["MHC_II_score_.WT."].replace(",",".")
    imp_binder = 0
    try:
        imp_binder = float(sc_wt)/float(sc_mut) >= 1.2
        return str(1) if imp_binder else str(0)
    except (ZeroDivisionError, ValueError) as e:
        return "NA"


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
    if mhc == "mhcI":
        mhc_mut = props["MHC_I_epitope_.best_prediction."]
        mhc_wt = props["MHC_I_epitope_.WT."]
    elif mhc == "mhcII":
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




def position_in_anchor_position(props):
    '''
    This function determines if the mutation is located within an anchor position in mhc I.
    As an approximation, we assume that the second and the last position are anchor positions for all alleles.
    '''
    pos_mhcI = props["pos_MUT_MHCI"]
    pep_len = props["MHC_I_peptide_length_.best_prediction."]

    anchor = 0
    try:
        anchor = int(pos_mhcI) == int(pep_len) or int(pos_mhcI)==2
        return str(1) if anchor else str(0)
    except:
        return "NA"


if __name__ == '__main__':
    import sys
    basedir = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2"
    sys.path.append(basedir)
    import data_import

    file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/INPuT2/nonprogramm_files/20170713_IS_IM_data.complete.update_Dv10.csv"
    indel = False
    dat = data_import.import_dat_icam(file, indel)

    properties = {}
    for nam,char in zip(dat[0], dat[1][1]):
        properties[nam] = char


    print selfsimilarity(properties, "mhcI")
    print position_of_mutation_epitope(properties, "mhcI")
