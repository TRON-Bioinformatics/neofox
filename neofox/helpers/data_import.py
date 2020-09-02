#!/usr/bin/python
from logzero import logger


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
    # read neofox file
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
    logger.info("reading neofox done {} items; {}  columns".format(len(data), len(data[0])))
    header = [x.strip('"').strip("\r").strip('"') for x in header]
    data = [[y.strip('"').strip("\r").strip('"') for y in x] for x in data]
    return header, data
