from typing import List
from neofox.model.neoantigen import Mhc1, Zygosity


class MhcHelper:

    @staticmethod
    def get_homozygous_mhc1_alleles(mhc_isoforms: List[Mhc1]) -> List[str]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise retunrs empty list
        """
        return [
            a.name
            for m in mhc_isoforms
            for a in m.alleles
            if m.zygosity == Zygosity.HOMOZYGOUS
        ]

    @staticmethod
    def get_heterozygous_or_hemizygous_mhc1_alleles(mhc_isoforms: List[Mhc1]) -> List[str]:
        """
        Returns alleles that occur more than one time in list of patient alleles and hence are homozygous alleles.
        Otherwise retunrs empty list
        """
        return [
            a.name
            for m in mhc_isoforms
            for a in m.alleles
            if m.zygosity in [Zygosity.HETEROZYGOUS, Zygosity.HEMIZYGOUS]
        ]
