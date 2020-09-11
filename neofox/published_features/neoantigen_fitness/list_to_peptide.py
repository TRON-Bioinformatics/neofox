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
import sys


def list_to_fa(data, col, outfile):
    n = []
    with open(outfile, "w") as f:
        for i, d in enumerate(data):
            id = ">M_" + str(i + 1)
            f.write(id + "\n")
            n.append(id)
            f.write(d[col] + "\n")
    return n


if __name__ == "__main__":
    c = -1
    with open(sys.argv[1]) as f:
        for line in f:
            c += 1
            if c == 0:
                continue
            w = line.strip("\t").split(",")
            print(">M_" + str(c))
            print(w[1])
