import input


class AvailableAlleles(object):

    def __init__(self, references):
        self.available_mhc_i = self._load_available_hla_alleles(mhc=input.MHC_I, references=references)
        self.available_mhc_ii = self._load_available_hla_alleles(mhc=input.MHC_II, references=references)

    def _load_available_hla_alleles(self, references, mhc=input.MHC_I):
        """
        loads file with available hla alllels for netmhcpan4/netmhcIIpan prediction, returns set
        :type references: input.references.ReferenceFolder
        :type mhc: str
        :rtype list:
        """
        if mhc == input.MHC_II:
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
