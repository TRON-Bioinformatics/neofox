#!/usr/bin/python
#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
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
    logger.info("reading input done {} items; {}  columns".format(len(data), len(data[0])))
    header = [x.strip('"').strip("\r").strip('"') for x in header]
    data = [[y.strip('"').strip("\r").strip('"') for y in x] for x in data]
    return header, data
