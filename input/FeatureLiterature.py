#!/usr/bin/env python

# import modules
import math

from logzero import logger




def number_of_mismatches(epitope_wild_type, epitope_mutation):
    """
    This function calculates the number of mismatches between the wt and the mutated epitope
    """
    p1 = 0
    try:
        for i, aa in enumerate(epitope_mutation):
            if aa != epitope_wild_type[i]:
                p1 += 1
        return str(p1)
    except IndexError:
        return "NA"


def match_in_proteome(sequence, db):
    """
    This function checks if the mutated epitope has an exact match in a protein database (uniprot)
    Returns 0 if mutation is present in proteome and 1 if it not present
    """
    try:
        seq_in_db = [sequence in entry for entry in db]
        return "0" if any(seq_in_db) else "1"
    except:
        return "NA"


def calc_logistic_function(mhc_score):
    '''Calculates negative logistic function given mhc score
    '''
    try:
        return float(1 / (1 + math.exp(5 * (float(mhc_score) - 2))))
    except (OverflowError, ValueError) as e:
        return "NA"


def calc_priority_score(vaf_tumor, vaf_rna, transcript_expr, no_mismatch, score_mut, score_wt, mut_not_in_prot):
    """
    This function calculates the Priority Score using parameters for mhc I.
    """
    L_mut = calc_logistic_function(score_mut)
    L_wt = calc_logistic_function(score_wt)
    priority_score = "NA"
    try:
        if vaf_tumor not in ["-1", "NA"]:
            priority_score = (L_mut * float(vaf_tumor) * math.tanh(float(transcript_expr))) * (
                    float(mut_not_in_prot) * (1 - 2 ** (-float(no_mismatch)) * L_wt))
        else:
            priority_score = (L_mut * float(vaf_rna) * math.tanh(float(transcript_expr))) * (
                    float(mut_not_in_prot) * (1 - 2 ** (-float(no_mismatch)) * L_wt))
    except (TypeError, ValueError) as e:
        pass
    return str(priority_score)


