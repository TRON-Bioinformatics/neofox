#!/usr/bin/python

from logzero import logger
import pandas as pd
from input.model.schema_conversion import SchemaConverter
from input.model.neoantigen import Patient


def import_dat_icam(in_file, indel=False):
    """
    This function imports a result table from the iCAM pipeline
    Parameter indel (True/False) indicates if output contains indels or SNVs
    """
    file_format = 0
    data = []
    header = []
    c = 0
    # check format of file
    if "txt" in in_file or "transcript" in in_file:
        file_format = "txt"
    elif "csv" in in_file:
        file_format = "csv"
    else:
        logger.warn("Table should be in csv or txt format!!")
    # read input file
    with open(in_file) as f:
        if file_format == "csv":
            for line in f:
                c += 1
                w = line.replace('"', "").replace(",", ".").strip("\n").split(";")
                if c == 1:
                    header = w
                    subst_col = header.index("substitution")
                    continue
                # skip indels --> on position 11 information about
                indel_inf = w[subst_col]
                is_indel_inf_null = has_indel_inf_null_value(indel_inf)
                if indel and is_indel_inf_null:
                    data.append(w)
                elif not is_indel_inf_null:
                    data.append(w)
        elif file_format == "txt":
            for line in f:
                if line.startswith("#"):
                    continue
                c += 1
                w = [x.replace('"', '').replace(",", ".") for x in line.strip("\n").split("\t")]
                if c == 1:
                    header = w
                    subst_col = header.index("substitution")
                    continue
                indel_inf = w[subst_col]
                is_indel_inf_null = has_indel_inf_null_value(indel_inf)
                if indel and is_indel_inf_null:
                    data.append(w)
                elif not is_indel_inf_null:
                    data.append(w)

    logger.info("reading input done {} items; {} columns".format(len(data), len(data[0])))
    header = [x.strip('"').strip("\r").strip('"') for x in header]
    data = [[y.strip('"').strip("\r").strip('"') for y in x] for x in data]
    return header, data


def has_indel_inf_null_value(indel_inf):
    return indel_inf == "" or "-" in indel_inf or "NA" in indel_inf


def import_dat_general(in_file):
    '''
    This function imports a csv or txt formated table
    '''
    file_format = 0
    data = []
    header = []
    c = 0
    # ceck format of file
    if in_file.endswith(".txt"):
        file_format = "txt"
    elif in_file.endswith(".csv"):
        file_format = "csv"
    else:
        logger.warn("Table should be in csv or txt format!!")
    # read input file
    with open(in_file) as f:
        if file_format == "csv":
            for line in f:
                c += 1
                w = line.replace('"', "").replace(",", ".").strip("\n").split(";")
                if c == 1:
                    header = w
                    continue
                else:
                    data.append(w)
        elif file_format == "txt":
            for line in f:
                if line.startswith("#"):
                    continue
                c += 1
                w = [x.replace('"', '').replace(",", ".") for x in line.strip("\n").split("\t")]
                if c == 1:
                    header = w
                    continue
                else:
                    data.append(w)
    logger.info("reading input done {} items; {}  columns".format(len(data), len(data[0])))
    header = [x.strip('"').strip("\r").strip('"') for x in header]
    data = [[y.strip('"').strip("\r").strip('"') for y in x] for x in data]
    return header, data


def change_col_names(header, data):
    """This function changes the names of columns if table hat not been importet to R before."""
    if "MHC_I_peptide_length_(best_prediction)" in header:
        mut_long = header.index("+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)")
        wt_long = header.index("[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)")
        lenI_ind = header.index("MHC_I_peptide_length_(best_prediction)")
        allI_ind = header.index("MHC_I_allele_(best_prediction)")
        epiI_ind = header.index("MHC_I_epitope_(best_prediction)")
        epiIwt_ind = header.index("MHC_I_epitope_(WT)")
        scI_ind = header.index("MHC_I_score_(best_prediction)")
        scIwt_ind = header.index("MHC_I_score_(WT)")
        lenII_ind = header.index("MHC_II_peptide_length_(best_prediction)")
        allII_ind = header.index("MHC_II_allele_(best_prediction)")
        epiII_ind = header.index("MHC_II_epitope_(best_prediction)")
        epiIIwt_ind = header.index("MHC_II_epitope_(WT)")
        scII_ind = header.index("MHC_II_score_(best_prediction)")
        scIIwt_ind = header.index("MHC_II_score_(WT)")

        header[mut_long] = "X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."
        header[wt_long] = "X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."
        header[lenI_ind] = "MHC_I_peptide_length_.best_prediction."
        header[allI_ind] = "MHC_I_allele_.best_prediction."
        header[epiI_ind] = "MHC_I_epitope_.best_prediction."
        header[epiIwt_ind] = "MHC_I_epitope_.WT."
        header[scI_ind] = "MHC_I_score_.best_prediction."
        header[scIwt_ind] = "MHC_I_score_.WT."
        header[lenII_ind] = "MHC_II_peptide_length_.best_prediction."
        header[allII_ind] = "MHC_II_allele_.best_prediction."
        header[epiII_ind] = "MHC_II_epitope_.best_prediction."
        header[epiIIwt_ind] = "MHC_II_epitope_.WT."
        header[scII_ind] = "MHC_II_score_.best_prediction."
        header[scIIwt_ind] = "MHC_II_score_.WT."

    return header, data


def import_patients_data(patients_file):
    """
    :param patients_file: the file to patients data CSV file
    :type patients_file: str
    :return: the parsed CSV into model objects
    :rtype: list[Patient]
    """
    split_comma_separated_list = lambda x: x.split(',')
    df = pd.read_csv(
        patients_file,
        sep='\t',
        converters={'mhcIAlleles': split_comma_separated_list,
                    'mhcIIAlleles': split_comma_separated_list,
                    # TODO: remove this conversion if this is fixed
                    #  https://github.com/danielgtaylor/python-betterproto/issues/96
                    'estimatedTumorContent': lambda x: float(x) if x != "NA" else x})
    return SchemaConverter.patient_metadata_csv2model(df)
