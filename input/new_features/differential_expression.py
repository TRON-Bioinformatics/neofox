#!/usr/bin/python

import math


def mean_of_list(list_numbs):
    '''
    This function takes a list of numbers as input and calculates the mean
    '''
    mean_list = float(sum(list_numbs) / len(list_numbs))
    return mean_list


def sd_of_list(list_numbs):
    '''
    This function takes a list of numbers as input and calculates the standard deviation
    '''
    mean = mean_of_list(list_numbs)
    sd_list = float((sum([((x - mean) ** 2) for x in list_numbs]) / len(list_numbs)) ** 0.5)
    return (sd_list)


def add_rna_reference(props, reference_dat, i):
    '''
    This function takes the output of load_rna_expression_reference function and appends the values to the epitope data
    i = (0,1,2) --> ("mean_ref_expression","sd_ref_expression", "sum_ref_expression")
    '''
    if "gene" in props:
        gene = props["gene"]
    else:
        gene = props["gene.x"]

    if gene in reference_dat:
        return str(reference_dat[gene][i])
    else:
        return "NA"


def fold_change(props):
    '''
    This function determines the classical fold change between tumour and reference transcript expression.
    Log2(expr in tumor / expr in reference)
    '''
    expr_tumour = props["transcript_expression"]
    expr_reference = props["mean_ref_expression"]
    try:
        return str(math.log(float(expr_tumour) / float(expr_reference), 2))
    except (ValueError, ZeroDivisionError) as e:
        return "NA"


def percentile_calc(props):
    '''
    This function calculates the expression difference between tumour and reference data in form of a percentile value.
    expr in tumor * 100 / (sum of expr in ref tissue + 1)
    '''
    expr_tumour = props["transcript_expression"]
    expr_reference = props["sum_ref_expression"]
    try:
        return str((float(expr_tumour) * 100) / (float(expr_reference) + 1))
    except (ValueError, ZeroDivisionError) as e:
        return "NA"


def pepper_calc(props):
    '''
    This function calculates the expression difference between tumour and reference data based on Pepper publication, in a z-score similar manner.
    expr in tumour - mean epxr in reference tissue / standard deviation of expression in refernce
    '''
    expr_tumour = props["transcript_expression"]
    expr_reference = props["mean_ref_expression"]
    expr_reference_sd = props["sd_ref_expression"]
    try:
        return str((float(expr_tumour) - float(expr_reference)) / float(expr_reference_sd))
    except (ValueError, ZeroDivisionError) as e:
        return "NA"


if __name__ == '__main__':
    import sys
    import data_import

    ref_file = sys.argv[1]
    # "/projects/CM27_IND_patients/GTEX_normal_tissue_data/Skin .csv"
    ref_list = load_rna_expression_reference(ref_file)
    f = sys.argv[2]
    data = data_import.import_dat_icam(f)
    dat_merged = merge_data_reference(data, ref_list)
    print(wrapper_diff_expr(dat_merged)[0])

    # write_ouptut_to_file(dat_epi,header)
