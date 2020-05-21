#!/usr/bin/env python

import tempfile

from logzero import logger

from input.helpers import properties_manager


class MixMHCpred:

    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.all_peptides = "NA"
        self.all_scores = "NA"
        self.all_ranks = "NA"
        self.all_alleles = "NA"
        self.best_peptide = "NA"
        self.best_score = "NA"
        self.best_rank = "NA"
        self.best_allele = "NA"
        self.best_peptide_wt = "NA"
        self.best_score_wt = "NA"
        self.best_rank_wt = "NA"
        self.difference_score_mut_wt = "NA"

    def mut_position_xmer_seq(self, props):
        '''returns position of mutation in xmer sequence
        '''
        xmer_wt = props["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        xmer_mut = props["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        p1 = -1
        if len(xmer_wt) == len(xmer_mut):
            p1 = -1
            for i, aa in enumerate(xmer_mut):
                if aa != xmer_wt[i]:
                    p1 = i + 1
        else:
            p1 = 0
            # in case sequences do not have same length
            for a1, a2 in zip(xmer_wt, xmer_mut):
                if a1 == a2:
                    p1 += 1
        return str(p1)

    def generate_nmers(self, props, list_lengths, mut=True):
        ''' generates peptides covering mutation of all lengts that are provided. Returns peptides as list
        '''
        xmer_wt = props["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        xmer_mut = props["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        list_peptides = []
        pos_mut = int(self.mut_position_xmer_seq(props))
        long_seq = xmer_mut if mut else xmer_wt
        for l in list_lengths:
            l = int(l)
            start_first = pos_mut - (l)
            starts = []
            for s in range(l):
                starts.append(int(start_first + s))
            ends = []
            [ends.append(int(s + (l))) for s in starts]
            for s, e in zip(starts, ends):
                list_peptides.append(long_seq[s:e])
        return list_peptides

    def generate_fasta(self, seqs, tmpfile):
        ''' Writes seqs given in seqs list into fasta file
        '''
        counter = 0
        with open(tmpfile, "w") as f:
            for seq in seqs:
                id = "".join([">seq", str(counter)])
                f.write(id + "\n")
                f.write(seq + "\n")
                counter += 1

    def mixmhcprediction(self, hla_alleles, tmpfasta, outtmp):
        ''' Performs MixMHCpred prediction for desired hla allele and writes result to temporary file.
        '''
        allels_for_prediction = []
        for allele in hla_alleles:
            allele = allele.replace("*", "")
            allele = allele.replace("HLA-", "")
            allels_for_prediction.append(allele)
        hla_allele = ",".join(allels_for_prediction)
        self.runner.run_command(cmd=[
            self.configuration.mix_mhc_pred,
            "-a", hla_allele,
            "-i", tmpfasta,
            "-o", outtmp])

    def read_mixmhcpred(self, outtmp):
        '''imports output of MixMHCpred prediction
        '''
        counter = 0
        header = []
        dat = []
        with open(outtmp) as f:
            for line in f:
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith("#"):
                        continue
                    if line.startswith("Peptide"):
                        counter += 1
                        header = line.split()
                        continue
                    line = line.split()
                    dat.append(line)
        return header, dat

    def extract_best_per_pep(self, pred_dat):
        '''extract info of best allele prediction for all potential ligands per muatation
        '''
        head = pred_dat[0]
        dat = pred_dat[1]
        peps = []
        scores = []
        alleles = []
        ranks = []
        pepcol = head.index("Peptide")
        scorecol = head.index("Score_bestAllele")
        allelecol = head.index("BestAllele")
        rankcol = head.index("%Rank_bestAllele")
        min_value = -1000000000000000000
        for ii, i in enumerate(dat):
            col_of_interest = [i[pepcol], i[scorecol], i[rankcol], i[allelecol]]
            # all potential peptides per mutation --> return ditionary
            peps.append(i[pepcol])
            scores.append(i[scorecol])
            ranks.append(i[rankcol])
            alleles.append(i[allelecol])
        return {"Peptide": peps, "Score_bestAllele": scores, "BestAllele": alleles, "%Rank_bestAllele": ranks}

    def extract_best_peptide_per_mutation(self, pred_dat):
        '''extract best predicted ligand per mutation
        '''
        head = pred_dat[0]
        dat = pred_dat[1]
        peps = []
        scores = []
        alleles = []
        ranks = []
        pepcol = head.index("Peptide")
        scorecol = head.index("Score_bestAllele")
        allelecol = head.index("BestAllele")
        rankcol = head.index("%Rank_bestAllele")
        min_value = -1000000000000000000
        for ii, i in enumerate(dat):
            col_of_interest = [str(i[pepcol]), str(i[scorecol]), str(i[rankcol]), str(i[allelecol])]
            # best ligand per mutation
            if float(i[scorecol]) > float(min_value):
                min_value = i[scorecol]
                min_pep = col_of_interest
        head_new = ["Peptide", "Score_bestAllele", "%Rank_bestAllele", "BestAllele"]
        return head_new, min_pep

    def add_best_epitope_info(self, epitope_tuple, column_name):
        '''returns desired information of prediction of best epitope from netmhcpan output;
        e.g. "%Rank": MHC I score, "HLA": HLA allele, "Icore": best epitope
        '''
        dat_head = epitope_tuple[0]
        dat = epitope_tuple[1]
        val = dat_head.index(column_name)
        try:
            return dat[val]
        except IndexError:
            return "NA"

    def extract_WT_for_best(self, props, best_mut_seq):
        '''extracts the corresponding WT epitope for best predicted mutated epitope
        '''
        xmer_wt = props["X.WT._..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        xmer_mut = props["X..13_AA_.SNV._._.15_AA_to_STOP_.INDEL."]
        start = xmer_mut.find(best_mut_seq)
        l = len(best_mut_seq)
        wt_epi = xmer_wt[start:(start + l)]
        return (wt_epi)

    def extract_WT_info(self, epitope_tuple, column_name):
        '''
        '''
        dat_head = epitope_tuple[0]
        dat = epitope_tuple[1][0]
        val = dat_head.index(column_name)
        try:
            return dat[val]
        except IndexError:
            return "NA"

    def difference_score(self, mut_score, wt_score):
        ''' calcualated difference in MixMHCpred scores between mutated and wt
        '''
        try:
            return str(float(mut_score) - float(wt_score))
        except ValueError:
            return "NA"

    def main(self, props_dict, dict_patient_hla):
        '''Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR
        '''
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix="tmp_sequence_", suffix=".fasta", delete=False)
        tmp_fasta = tmp_fasta_file.name
        tmp_prediction_file = tempfile.NamedTemporaryFile(prefix="mixmhcpred", suffix=".txt", delete=False)
        tmp_prediction = tmp_prediction_file.name
        seqs = self.generate_nmers(props_dict, [8, 9, 10, 11])
        self.generate_fasta(seqs, tmp_fasta)
        alleles = properties_manager.get_hla_allele(props_dict, dict_patient_hla)
        self.mixmhcprediction(alleles, tmp_fasta, tmp_prediction)
        pred = self.read_mixmhcpred(tmp_prediction)
        try:
            pred_all = self.extract_best_per_pep(pred)
        except ValueError:
            pred_all = {}
        if len(pred_all) > 0:
            pred_best = self.extract_best_peptide_per_mutation(pred)
            self.best_peptide = self.add_best_epitope_info(pred_best, "Peptide")
            self.best_score = self.add_best_epitope_info(pred_best, "Score_bestAllele")
            self.best_rank = self.add_best_epitope_info(pred_best, "%Rank_bestAllele")
            self.best_allele = self.add_best_epitope_info(pred_best, "BestAllele")
            self.best_allele = self.add_best_epitope_info(pred_best, "BestAllele")
            self.all_peptides = "|".join(pred_all["Peptide"])
            self.all_scores = "|".join(pred_all["Score_bestAllele"])
            self.all_ranks = "|".join(pred_all["%Rank_bestAllele"])
            self.all_alleles = "|".join(pred_all["BestAllele"])
            # prediction of for wt epitope that correspond to best epitope
            wt = self.extract_WT_for_best(props_dict, self.best_peptide)
            wt_list = [wt]
            tmp_fasta_file = tempfile.NamedTemporaryFile(prefix="tmp_sequence_wt_", suffix=".fasta", delete=False)
            tmp_fasta = tmp_fasta_file.name
            tmp_prediction_file = tempfile.NamedTemporaryFile(prefix="mixmhcpred_wt_", suffix=".txt", delete=False)
            tmp_prediction = tmp_prediction_file.name
            self.generate_fasta(wt_list, tmp_fasta)
            self.mixmhcprediction(alleles, tmp_fasta, tmp_prediction)
            pred_wt = self.read_mixmhcpred(tmp_prediction)
            logger.debug(pred_wt)
            self.best_peptide_wt = self.extract_WT_info(pred_wt, "Peptide")
            score_wt_of_interest = "_".join(["Score", self.best_allele])
            rank_wt_of_interest = "_".join(["%Rank", self.best_allele])
            self.best_score_wt = self.extract_WT_info(pred_wt, score_wt_of_interest)
            self.best_rank_wt = self.extract_WT_info(pred_wt, rank_wt_of_interest)
            # difference in scores between mut and wt
            self.difference_score_mut_wt = self.difference_score(self.best_score, self.best_score_wt)


# if __name__ == '__main__':
#
#     from input import predict_all_epitopes, epitope
#
#     # test with ott data set
#     file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
#     hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
#     # test inest data set
#     # file = "/flash/projects/WP3/AnFranziska/AnFranziska/head_seqs.txt"
#     # hla_file = "/flash/projects/WP3/AnFranziska/AnFranziska/alleles.csv"
#     dat = data_import.import_dat_icam(file, False)
#     if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
#         dat = data_import.change_col_names(dat)
#     # available MHC alleles
#     set_available_mhc = predict_all_epitopes.Bunchepitopes().load_available_hla_alleles()
#     # hla allele of patients
#     patient_hlaI = predict_all_epitopes.Bunchepitopes().load_patient_hla_I_allels(hla_file)
#     patient_hlaII = predict_all_epitopes.Bunchepitopes().load_patient_hla_II_allels(hla_file)
#
#     print(patient_hlaI)
#     print(patient_hlaII)
#
#     for ii, i in enumerate(dat[1]):
#         if ii < 10:
#             print(ii)
#             dict_epi = epitope.Epitope()
#             dict_epi.init_properties(dat[0], dat[1][ii])
#             prediction = MixMHCpred()
#             prediction.main(dict_epi.properties, patient_hlaI)
#             attrs = vars(prediction)
#             print(attrs)
