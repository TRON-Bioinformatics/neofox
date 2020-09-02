#!/usr/bin/python

import math


# TODO: review if this can be replaced by out of the box implementations

def mean_of_list(list_numbs):
    '''
    This function takes a list of numbers as neofox and calculates the mean
    '''
    mean_list = float(sum(list_numbs) / len(list_numbs))
    return mean_list


def sd_of_list(list_numbs):
    '''
    This function takes a list of numbers as neofox and calculates the standard deviation
    '''
    mean = mean_of_list(list_numbs)
    sd_list = float((sum([((x - mean) ** 2) for x in list_numbs]) / len(list_numbs)) ** 0.5)
    return (sd_list)


def fold_change(expression_tumor, expression_reference):
    """
    This function determines the classical fold change between tumour and reference transcript expression.
    Log2(expr in tumor / expr in reference)
    """
    try:
        return str(math.log(float(expression_tumor) / float(expression_reference), 2))
    except (TypeError, ValueError, ZeroDivisionError) as e:
        return "NA"


def percentile_calc(expression_tumor, expression_reference_sum):
    """
    This function calculates the expression difference between tumour and reference data in form of a percentile value.
    expr in tumor * 100 / (sum of expr in ref tissue + 1)
    """
    try:
        return str((float(expression_tumor) * 100) / (float(expression_reference_sum) + 1))
    except (TypeError, ValueError, ZeroDivisionError) as e:
        return "NA"


def pepper_calc(expression_tumor, expression_reference, expression_reference_sd):
    """
    This function calculates the expression difference between tumour and reference data based on Pepper publication, in a z-score similar manner.
    expr in tumour - mean epxr in reference tissue / standard deviation of expression in refernce
    """
    try:
        return str((float(expression_tumor) - float(expression_reference)) / float(expression_reference_sd))
    except (TypeError, ValueError, ZeroDivisionError) as e:
        return "NA"
