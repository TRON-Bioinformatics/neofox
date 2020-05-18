#!/usr/bin/python

from logzero import logger


def import_dat_icam(in_file, indel):
    """
    This function imports a result table from the iCAM pipeline
    Parameter indel (True/False) indicates if output contains indels or SNVs
    """
    file_format = 0
    c = 0
    data = []
    header = []
    c = 0
    # ceck format of file
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


def import_as_dict(in_file, key):
    '''Reads csv file and returns dictionary. specifiy what is key, remaining is value
    '''
    d = {}
    c = 0
    with open(in_file) as f:
        for line in f:
            w = line.replace('"', "").rstrip("\r\n").split(";")
            if c == 0:
                header = w
                keyIndex = int(header.index(key))
                c += 1
                continue
            if c > 0:
                d[w[keyIndex]] = w[:keyIndex] + w[keyIndex + 1:]
            c += 1
    return d


def import_allele_file(allele_file):
    '''imports allele.csv file in form of dictionary'''
    d = {}
    with open(allele_file) as f:
        for line in f:
            w = line.replace('"', "").rstrip("\r\n").rstrip(';').split(";")
            d[w[0] + "_" + w[1]] = w[2:]
    return d


def get_header_from_tuple(tuple_dat_head):
    '''
    get columnnames/header of data frame stored in tuple
    '''
    return tuple_dat_head[0]


def get_data_from_tuple(tuple_dat_head):
    '''
    get columnnames/header of data frame stored in tuple
    '''
    return tuple_dat_head[1]


def change_col_names(tuple_dat_head):
    """This function changes the names of columns if table hat not been importet to R before."""
    dat_new = tuple_dat_head[1]
    head_new = tuple_dat_head[0]
    if "MHC_I_peptide_length_(best_prediction)" in head_new:
        mut_long = head_new.index("+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)")
        wt_long = head_new.index("[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)")
        lenI_ind = head_new.index("MHC_I_peptide_length_(best_prediction)")
        allI_ind = head_new.index("MHC_I_allele_(best_prediction)")
        epiI_ind = head_new.index("MHC_I_epitope_(best_prediction)")
        epiIwt_ind = head_new.index("MHC_I_epitope_(WT)")
        scI_ind = head_new.index("MHC_I_score_(best_prediction)")
        scIwt_ind = head_new.index("MHC_I_score_(WT)")
        lenII_ind = head_new.index("MHC_II_peptide_length_(best_prediction)")
        allII_ind = head_new.index("MHC_II_allele_(best_prediction)")
        epiII_ind = head_new.index("MHC_II_epitope_(best_prediction)")
        epiIIwt_ind = head_new.index("MHC_II_epitope_(WT)")
        scII_ind = head_new.index("MHC_II_score_(best_prediction)")
        scIIwt_ind = head_new.index("MHC_II_score_(WT)")

        head_new[mut_long] = "X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."
        head_new[wt_long] = "X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."
        head_new[lenI_ind] = "MHC_I_peptide_length_.best_prediction."
        head_new[allI_ind] = "MHC_I_allele_.best_prediction."
        head_new[epiI_ind] = "MHC_I_epitope_.best_prediction."
        head_new[epiIwt_ind] = "MHC_I_epitope_.WT."
        head_new[scI_ind] = "MHC_I_score_.best_prediction."
        head_new[scIwt_ind] = "MHC_I_score_.WT."
        head_new[lenII_ind] = "MHC_II_peptide_length_.best_prediction."
        head_new[allII_ind] = "MHC_II_allele_.best_prediction."
        head_new[epiII_ind] = "MHC_II_epitope_.best_prediction."
        head_new[epiIIwt_ind] = "MHC_II_epitope_.WT."
        head_new[scII_ind] = "MHC_II_score_.best_prediction."
        head_new[scIIwt_ind] = "MHC_II_score_.WT."

        try:
            patid = head_new.index("patient")
            head_new[patid] = "patient.id"
        except ValueError:
            pass

        return head_new, dat_new

    else:
        return head_new, dat_new


def subst_semicolon(tuple_dat_head):
    '''
    This function substitutes any semilicon by "_", since output is in csv format --> problems when importing into R
    '''
    dat_new = tuple_dat_head[1]
    head_new = tuple_dat_head[0]
    data = []
    for ii, i in enumerate(dat_new):
        new = []
        for element in i:
            if ";" in element:
                new.append(element.replace(";", "_"))
            else:
                new.append(element)
        data.append(new)
    logger.info(" ';' substituted by '_' ")
    return head_new, data


def append_patient(tuple_dat_head, in_file):
    '''
    There is no column indicating the column name in Tesla icam outputs. This function adds a columns with with the patient id.
    '''
    dat_new = tuple_dat_head[1]
    head_new = tuple_dat_head[0]
    if "patient" in head_new or "patient.x" in head_new:
        return head_new, dat_new
    else:
        pat = f.split("/")
        pat = pat[len(pat) - 1]
        pat = pat.split("_")[0]
        for ii, i in enumerate(dat_new):
            dat_new[ii].append(pat)
        head_new.append("patient")
        return head_new, dat_new


if __name__ == '__main__':
    import sys

    f = sys.argv[1]
    dat = import_dat(f)
    dat = change_col_names(dat)
    dat = append_patient(dat, f)
    print(subst_semicolon(dat))
    # write_ouptut_to_file(dat)
else:
    import sys
