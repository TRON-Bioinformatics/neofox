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
import neofox


class AvailableAlleles(object):

    def __init__(self, references):
        self.available_mhc_i = self._load_available_hla_alleles(mhc=neofox.MHC_I, references=references)
        self.available_mhc_ii = self._load_available_hla_alleles(mhc=neofox.MHC_II, references=references)

    def _load_available_hla_alleles(self, references, mhc=neofox.MHC_I):
        """
        loads file with available hla alllels for netmhcpan4/netmhcIIpan prediction, returns set
        :type references: neofox.references.ReferenceFolder
        :type mhc: str
        :rtype list:
        """
        if mhc == neofox.MHC_II:
            fileMHC = references.available_mhc_ii
        else:
            fileMHC = references.available_mhc_i
        set_available_mhc = set()
        with open(fileMHC) as f:
            for line in f:
                set_available_mhc.add(line.strip())
        return set_available_mhc

    def get_available_mhc_i(self):
        return self.available_mhc_i

    def get_available_mhc_ii(self):
        return self.available_mhc_ii
