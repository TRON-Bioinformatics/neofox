#!/usr/bin/env python

from logzero import logger

from input.helpers import data_import
from input.helpers.epitope_helper import EpitopeHelper
from input.netmhcpan4.abstract_netmhcpan_predictor import AbstractNetMhcPanPredictor


class NetMhcPanPredictor(EpitopeHelper, AbstractNetMhcPanPredictor):

    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.mhc_score = "NA"
        self.epitope = "NA"
        self.allele = "NA"
        self.directed_to_TCR = "NA"
        self.affinity = "NA"
        self.affinity_epitope = "NA"
        self.affinity_allele = "NA"
        self.affinity_directed_to_TCR = "NA"
        self.mhcI_score_9mer = "NA"
        self.mhcI_score_allele_9mer = "NA"
        self.mhcI_score_epitope_9mer = "NA"
        self.mhcI_affinity_9mer = "NA"
        self.mhcI_affinity_allele_9mer = "NA"
        self.mhcI_affinity_epitope_9mer = "NA"

    def _mhc_allele_in_netmhcpan_available(self, allele, set_available_mhc):
        '''checks if mhc prediction is possible for given hla allele
        '''
        return allele in set_available_mhc

    def check_format_allele(self, allele):
        """
        sometimes genotyping may be too detailed. (e.g. HLA-DRB1*04:01:01 should be HLA-DRB1*04:01)
        :param allele: HLA-allele
        :return: HLA-allele in correct format
        """
        # TODO: was added to netMHCIIpan too --> combine
        if allele.count(":") > 1:
            allele_correct = ":".join(allele.split(":")[0:2])
        else:
            allele_correct = allele
        return allele_correct


    def mhc_prediction(self, hla_alleles, set_available_mhc, tmpfasta, tmppred):
        ''' Performs netmhcpan4 prediction for desired hla allele and writes result to temporary file.
        '''
        allels_for_prediction = []
        for allele in hla_alleles:
            allele = self.check_format_allele(allele)
            allele = allele.replace("*", "")
            if self._mhc_allele_in_netmhcpan_available(allele, set_available_mhc):
                allels_for_prediction.append(allele)
            else:
                logger.info(allele + "not available")
        hla_allele = ",".join(allels_for_prediction)
        cmd = [
            self.configuration.net_mhc_pan,
            "-a", hla_allele,
            "-f", tmpfasta,
            "-BA"]
        lines, _ = self.runner.run_command(cmd)
        counter = 0
        # TODO: avoid writing a file here, just return some data structure no need to go to the file system
        with open(tmppred, "w") as f:
            for line in lines.splitlines():
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("#", "-", "HLA", "Prot")):
                        continue
                    if counter == 0 and line.startswith("Pos"):
                        counter += 1
                        line = line.split()
                        line = line[0:-1] if len(line) > 14 else line
                        f.write(";".join(line) + "\n")
                        continue
                    elif counter > 0 and line.startswith("Pos"):
                        continue
                    line = line.split()
                    line = line[0:-2] if len(line) > 14 else line
                    line = ";".join(line)
                    f.write(line + "\n")

    def filter_binding_predictions(self, position_xmer, tmppred):
        '''filters prediction file for predicted epitopes that cover mutations
        '''
        dat_prediction = data_import.import_dat_general(tmppred)
        dat = dat_prediction[1]
        dat_head = dat_prediction[0]
        dat_fil = []
        pos_epi = dat_head.index("Pos")
        epi = dat_head.index("Peptide")
        for ii, i in enumerate(dat):
            if self.epitope_covers_mutation(position_xmer, i[pos_epi], len(i[epi])):
                dat_fil.append(dat[ii])
        return dat_head, dat_fil

    def minimal_binding_score(self, prediction_tuple, rank=True):
        '''reports best predicted epitope (over all alleles). indicate by rank = true if rank score should be used. if rank = False, Aff(nM) is used
        '''
        # TODO: generalize this method with netmhcIIpan_prediction.py + change input
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        if rank:
            mhc_sc = dat_head.index("%Rank")
        else:
            mhc_sc = dat_head.index("Aff(nM)")
        max_score = float(1000000000000)
        row = []
        for ii, i in enumerate(dat):
            mhc_score = float(i[mhc_sc])
            if mhc_score < max_score:
                max_score = mhc_score
                row = i
        return dat_head, row

    def mutation_in_loop(self, position_xmer_list, epitope_tuple):
        """
        returns if mutation is directed to TCR (yes or no)
        """
        dat_head = epitope_tuple[0]
        dat_epi = epitope_tuple[1]
        pos_epi = dat_head.index("Pos")
        del_pos = dat_head.index("Gp")
        del_len = dat_head.index("Gl")
        directed_to_tcr_list = [False]
        for position_mutation_xmer in position_xmer_list:
            if del_pos > 0:
                pos = int(dat_epi[pos_epi])
                start = pos + int(dat_epi[del_pos]) - 1
                end = start + int(dat_epi[del_len])
                if start < position_mutation_xmer <= end:
                    directed_to_tcr_list.append("yes")
        directed_to_tcr = any(directed_to_tcr_list)
        return directed_to_tcr


    def filter_for_9mers(self, prediction_tuple):
        '''returns only predicted 9mers
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        seq_col = dat_head.index("Peptide")
        dat_9mers = []
        for ii, i in enumerate(dat):
            seq = i[seq_col]
            if len(seq) == 9:
                dat_9mers.append(i)
        return dat_head, dat_9mers

    def filter_for_WT_epitope(self, prediction_tuple, mut_seq, mut_allele, number_snv):
        '''returns wt epitope info for given mutated sequence. best wt that is allowed to bind to any allele of patient
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        seq_col = dat_head.index("Peptide")
        allele_col = dat_head.index("HLA")
        wt_epi = []
        for ii, i in enumerate(dat):
            wt_seq = i[seq_col]
            wt_allele = i[allele_col]
            if (len(wt_seq) == len(mut_seq)):
                numb_mismatch = self.hamming_check_0_or_1(mut_seq, wt_seq)
                if numb_mismatch <= number_snv:
                    wt_epi.append(i)
        dt = (dat_head, wt_epi)
        min = self.minimal_binding_score(dt)
        return (min)

    def filter_for_WT_epitope_position(self, prediction_tuple, mut_seq, position_epi_xmer):
        '''returns wt epitope info for given mutated sequence. best wt that is allowed to bind to any allele of patient
        '''
        dat_head = prediction_tuple[0]
        dat = prediction_tuple[1]
        seq_col = dat_head.index("Peptide")
        pos_col = dat_head.index("Pos")
        wt_epi = []
        for ii, i in enumerate(dat):
            wt_seq = i[seq_col]
            wt_pos = i[pos_col]
            if (len(wt_seq) == len(mut_seq)) & (wt_pos == position_epi_xmer):
                wt_epi.append(i)
        dt = (dat_head, wt_epi)
        min = self.minimal_binding_score(dt)
        return (min)
