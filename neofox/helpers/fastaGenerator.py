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

from Bio.SeqIO.FastaIO import SimpleFastaParser


def write_single_seq_fasta(seq, id, file_name):
    '''Returns fasta file with a single sequences.
    '''
    with open(file_name, "w") as f:
        id = "".join([">", id, "\n"])
        f.write(id)
        seq = "".join([seq, "\n"])
        f.write(seq)


def read_multiple_seqs_simple(fasta_file):
    '''This function reads a fasta file using simplefasta parser and returns dictionary with gene names as keys and protein sequences as values.
    '''
    database = {}
    with open(fasta_file) as handle:
        for record in SimpleFastaParser(handle):
            # record[0] = fasta header; record[1] = protein sequence
            database[record[0]] = record[1]
    return database


if __name__ == '__main__':
    seq = "ABDSF"
    id = "M1"
    # write_single_seq_fasta(seq, id, "temp.fasta")
    read_multiple_seqs_simple(
        "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/indels/best_WT/BLAST_DB/ligandome/DB_ligandome_safe.fasta")
