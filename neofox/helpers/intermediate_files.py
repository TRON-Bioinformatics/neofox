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
import tempfile


def create_temp_file(prefix=None, suffix=None, dir=None):
    temp_file = tempfile.NamedTemporaryFile(
        prefix=prefix, suffix=suffix, dir=dir, delete=False
    )
    return temp_file.name


def create_temp_fasta(sequences, prefix=None, comment_prefix="seq"):
    """
    Writes seqs given in seqs list into fasta file
    """
    fasta_temp_file = create_temp_file(prefix=prefix, suffix=".fasta")
    counter = 1
    with open(fasta_temp_file, "w") as f:
        for seq in sequences:
            _id = ">{comment_prefix}{index}".format(
                comment_prefix=comment_prefix, index=counter
            )
            f.write(_id + "\n")
            f.write(seq + "\n")
            counter += 1
    return fasta_temp_file


def create_temp_peptide(sequences, prefix=None):
    """
    Writes seqs given in seqs list into PEPTIDE format
    """
    pep_temp_file = create_temp_file(prefix=prefix, suffix=".pep")
    with open(pep_temp_file, "w") as f:
        for seq in sequences:
            f.write(seq + "\n")
    return pep_temp_file
