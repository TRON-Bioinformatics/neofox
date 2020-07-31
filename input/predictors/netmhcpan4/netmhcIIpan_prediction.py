#!/usr/bin/env python

import tempfile

from logzero import logger

from input.helpers import data_import
from input.helpers.epitope_helper import EpitopeHelper
from input.predictors.netmhcpan4.abstract_netmhcpan_predictor import AbstractNetMhcPanPredictor


class NetMhcIIPanPredictor(EpitopeHelper, AbstractNetMhcPanPredictor):

    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.mhcII_score = "NA"
        self.epitopeII = "-"
        self.alleleII = "NA"
        self.affinityII = "NA"
        self.affinity_epitopeII = "-"
        self.affinity_alleleII = "NA"


    def check_format_allele(self, allele):
        """
        sometimes genotyping may be too detailed. (e.g. HLA-DRB1*04:01:01 should be HLA-DRB1*04:01)
        :param allele: HLA-allele
        :return: HLA-allele in correct format
        """
        if allele.count(":") > 1:
            allele_correct = ":".join(allele.split(":")[0:2])
        else:
            allele_correct = allele
        return allele_correct



    def generate_mhcII_alelles_combination_list(self, hla_alleles, set_available_mhc):
        ''' given list of HLA II alleles, returns list of HLA-DRB1 (2x), all possible HLA-DPA1/HLA-DPB1 (4x) and HLA-DQA1/HLA-DPQ1 (4x)
        '''
        allels_for_prediction = []
        dqa_alleles = []
        dpa_alleles = []
        dqb_alleles = []
        dpb_alleles = []
        for allele in hla_alleles:
            allele = self.check_format_allele(allele)
            if allele.startswith("HLA-DRB1"):
                allele = allele.replace("HLA-", "").replace("*", "_").replace(":", "")
                logger.info(allele)
                if allele in set_available_mhc:
                    allels_for_prediction.append(allele)
                else:
                    logger.info(allele + "not available")
            else:
                allele = allele.replace("*", "").replace(":", "")
                if allele.startswith("HLA-DPA"):
                    dpa_alleles.append(allele)
                elif allele.startswith("HLA-DPB"):
                    dpb_alleles.append(allele)
                elif allele.startswith("HLA-DQA"):
                    dqa_alleles.append(allele)
                elif allele.startswith("HLA-DQB"):
                    dqb_alleles.append(allele)
        dp_alleles = ["-".join([x, y.replace("HLA-", "")]) for x in dpa_alleles for y in dpb_alleles]
        dq_alleles = ["-".join([x, y.replace("HLA-", "")]) for x in dqa_alleles for y in dqb_alleles]
        dp_dq_alleles = dp_alleles + dq_alleles
        for allele in dp_dq_alleles:
            if allele in set_available_mhc:
                allels_for_prediction.append(allele)
        logger.info(allels_for_prediction)
        return allels_for_prediction

    def mhcII_prediction(self, hla_alleles, set_available_mhc, tmpfasta, tmppred):
        ''' Performs netmhcIIpan prediction for desired hla alleles and writes result to temporary file.
        '''
        allels_for_prediction = self.generate_mhcII_alelles_combination_list(hla_alleles, set_available_mhc)
        hla_allele = ",".join(allels_for_prediction)
        tmp_folder = tempfile.mkdtemp(prefix="tmp_netmhcIIpan_")
        logger.debug(tmp_folder)
        lines, _ = self.runner.run_command([
            self.configuration.net_mhc2_pan,
            "-a", hla_allele,
            "-f", tmpfasta,
            "-tdir", tmp_folder,
            "-dirty"])
        logger.debug(lines)
        counter = 0
        # TODO: avoid writing a file here, just return some data structure no need to go to the file system
        with open(tmppred, "w") as f:
            for line in lines.splitlines():
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("#", "-", "Number", "Temporary")):
                        continue
                    if counter == 0 and line.startswith("Seq"):
                        counter += 1
                        line = line.split()
                        line = line[0:-1] if len(line) > 12 else line
                        f.write(";".join(line) + "\n")
                        continue
                    elif counter > 0 and line.startswith("Seq"):
                        continue
                    line = line.split()
                    line = line[0:-2] if len(line) > 11 else line
                    f.write(";".join(line) + "\n")

    def filter_binding_predictions(self, position_xmer_sequence, tmppred):
        """
        filters prediction file for predicted epitopes that cover mutations
        """
        header, data = data_import.import_dat_general(tmppred)
        dat_fil = []
        logger.debug(header)
        pos_epi = header.index("Seq")
        epi = header.index("Peptide")
        for ii, i in enumerate(data):
            if self.epitope_covers_mutation(position_xmer_sequence, i[pos_epi], len(i[epi])):
                dat_fil.append(data[ii])
        return header, dat_fil

    def minimal_binding_score(self, prediction_tuple, rank=True):
        '''reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used. if rank = False, Aff(nM) is used
        '''
        # TODO: generalize this method with netmhcpan_prediction.py + change input
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        if rank:
            mhc_sc = dat_head.index("%Rank")
        else:
            mhc_sc = dat_head.index("Affinity(nM)")
        epi = dat_head.index("Peptide")
        hla_allele = dat_head.index("Allele")
        max_score = float(1000000000)
        allele = "NA"
        epitope = "NA"
        row = []
        for ii, i in enumerate(dat):
            mhc_score = float(i[mhc_sc])
            if mhc_score < max_score:
                max_score = mhc_score
                row = i
        return dat_head, row

    def filter_for_wt_epitope_position(self, prediction_tuple, mut_seq, position_epi_xmer):
        '''returns wt epitope info for given mutated sequence. best wt that is allowed to bind to any allele of patient
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        seq_col = dat_head.index("Peptide")
        pos_col = dat_head.index("Seq")
        wt_epi = []
        for ii, i in enumerate(dat):
            wt_seq = i[seq_col]
            wt_pos = i[pos_col]
            if (len(wt_seq) == len(mut_seq)) & (wt_pos == position_epi_xmer):
                wt_epi.append(i)
        dt = (dat_head, wt_epi)
        min = self.minimal_binding_score(dt)
        return (min)