from typing import List
from neofox.model.neoantigen import Mhc1, Zygosity, MhcAllele, PredictedEpitope, Neoantigen
from neofox.model.mhc_parser import MhcParser

import pandas as pd
from pandas.errors import EmptyDataError
from logzero import logger

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

class MixMhcHelper:

    def __init__(self, mhc_parser: MhcParser):
        self.mhc_parser = mhc_parser

    @staticmethod
    def get_mixmhc_allele_representation(file_alleles, mhc_alleles: List[MhcAllele]):
        """
        loads file with available HLA II alleles for MixMHCpred and PRIME prediction, returns set
        :return:
        """
        alleles = pd.read_csv(
            file_alleles, sep="\t"
        )
        available_alleles =  set(alleles["Allele"])

        converted_mhc_alleles = list(map(
                    lambda x: "{gene}{group}{protein}".format(gene=x.gene, group=x.group, protein=x.protein),
                    mhc_alleles))

        not_available_alleles = list(set(converted_mhc_alleles).difference(available_alleles))
        if len(not_available_alleles) > 0:
            logger.warning(
                "MHC I alleles {} are not supported by MixMHCpred.".format(",".join(not_available_alleles))
            )

        return list(
            set(converted_mhc_alleles).intersection(available_alleles)
        )

    def parse_mixmhcpred_prime_output(self, mixmhc_prime_result) -> List[PredictedEpitope]:
        RANK = "%Rank_"
        PEPTIDE = "Peptide"
        SCORE = "Score_"

        parsed_results = []

        try:
            results = pd.read_csv(mixmhc_prime_result, sep="\t", comment="#")
        except EmptyDataError:
            logger.error("Results from MixMHCpred are empty, something went wrong")
            results = pd.DataFrame()

        mhc_alleles = set()
        for col in results.columns:
            # take out alleles and eliminate the column Score_bestAllele out of the set
            if col.startswith(SCORE) and not col.endswith('e'):
                allele = col.split('_')[-1]
                mhc_alleles.add(allele)

        for _, row in results.iterrows():
            # when MixMHCpred returns no results it provides a row with the peptide and NAs for other fields
            # pandas reads NAs as float nan. Skip these
            for allele in mhc_alleles:
                if isinstance(row[PEPTIDE], str):
                    score = str(SCORE + allele)
                    rank = str(RANK + allele)

                    parsed_results.append(
                        PredictedEpitope(
                            allele_mhc_i=self.mhc_parser.parse_mhc_allele(allele),
                            mutated_peptide=row[PEPTIDE],
                            affinity_mutated=float(row[score]),
                            rank_mutated=float(row[rank]),
                        ))
        return parsed_results