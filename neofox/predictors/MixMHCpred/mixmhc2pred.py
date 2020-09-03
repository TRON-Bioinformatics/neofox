#!/usr/bin/env python
from typing import List

from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.predictors.MixMHCpred.abstract_mixmhcpred import AbstractMixMHCpred
from neofox.helpers import intermediate_files


class MixMhc2Pred(AbstractMixMHCpred):

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self.available_alleles = self.load_available_allelles()
        self._initialise()

    def _initialise(self):
        self.all_peptides = None
        self.all_ranks = None
        self.all_alleles = None
        self.best_peptide = None
        self.best_rank = None
        self.best_allele = None
        self.best_peptide_wt = None
        self.best_rank_wt = None

    def load_available_allelles(self):
        """
        loads file with available hla alllels for MixMHC2pred prediction, returns set
        :return:
        """
        path_to_HLAII_file = self.configuration.mix_mhc2_pred_alleles_list
        avail_alleles = []
        with open(path_to_HLAII_file) as f:
            for line in f:
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("L", "A")):
                        continue
                    line1 = line.split()[0]
                    if line1 is not None:
                        avail_alleles.append(line1)
        return avail_alleles

    def prepare_dq_dp(self, list_alleles):
        ''' returns patient DQ/DP alleles that are relevant for prediction
        '''
        list_alleles_pairs = ["__".join([p1, p2]) for p1 in list_alleles for p2 in list_alleles if p1 != p2]
        list_alleles_triplets = ["__".join([p1, p2, p3]) for p1 in list_alleles for p2 in list_alleles for p3 in
                                 list_alleles if p1 != p2 and p1 != p3 and p2 != p3]
        list_alleles_all = list_alleles_pairs + list_alleles_triplets
        alleles4pred = [allele for allele in list_alleles_all if allele in self.available_alleles]
        return (alleles4pred)

    def hlaIIallels2prediction(self, hla_alleles):
        ''' prepares list of hla alleles for prediction
        '''
        allels_for_prediction = []
        alleles_dq = []
        alleles_dp = []
        # print hla_alleles
        for allele in hla_alleles:
            # print allele
            allele = allele.replace("*", "_").replace(":", "_").replace("HLA-", "")
            if allele.startswith("DR"):
                if allele in self.available_alleles:
                    allels_for_prediction.append(allele)
            elif allele.startswith("DP"):
                alleles_dp.append(allele)
            elif allele.startswith("DQ"):
                alleles_dq.append(allele)
        alleles_dp4pred = self.prepare_dq_dp(alleles_dp)
        alleles_dq4pred = self.prepare_dq_dp(alleles_dq)
        allels_for_prediction = allels_for_prediction + alleles_dq4pred + alleles_dp4pred
        hla_allele = " ".join(allels_for_prediction)
        # print hla_allele
        return hla_allele

    def mixmhc2prediction(self, hla_alleles, tmpfasta, outtmp, wt=False):
        ''' Performs MixMHC2pred prediction for desired hla allele and writes result to temporary file.
        '''
        if not wt:
            hla_allele = self.hlaIIallels2prediction(hla_alleles)
        elif wt:
            # use best allele from mutated seq prediction
            hla_allele = hla_alleles[0]
        cmd = [
            self.configuration.mix_mhc2_pred,
            "-a", hla_allele,
            "-i", tmpfasta,
            "-o", outtmp]
        self.runner.run_command(cmd)


    def extract_best_per_pep(self, pred_dat):
        '''extract info of best allele prediction for all potential ligands per muatation
        '''
        head = pred_dat[0]
        dat = pred_dat[1]
        peps = []
        alleles = []
        ranks = []
        pepcol = head.index("Peptide")
        allelecol = head.index("BestAllele")
        rankcol = head.index("%Rank")
        for entry in sorted(dat, key=lambda x: (float(x[rankcol]), x[allelecol])):
            # all potential peptides per mutation --> return ditionary
            peps.append(entry[pepcol])
            ranks.append(entry[rankcol])
            alleles.append(entry[allelecol])
        return {"Peptide": peps, "BestAllele": alleles, "%Rank": ranks}

    def extract_best_peptide_per_mutation(self, pred_dat):
        '''extract best predicted ligand per mutation
        '''
        head = pred_dat[0]
        dat = pred_dat[1]
        pepcol = head.index("Peptide")
        allelecol = head.index("BestAllele")
        rankcol = head.index("%Rank")
        min_value = 1000000000000000000
        for ii, i in enumerate(dat):
            col_of_interest = [str(i[pepcol]), str(i[rankcol]), str(i[allelecol])]
            # best ligand per mutation
            if float(i[rankcol]) < float(min_value):
                min_value = i[rankcol]
                min_pep = col_of_interest
        head_new = ["Peptide", "%Rank", "BestAllele"]
        return head_new, min_pep

    def import_available_HLAII_alleles(self, path_to_HLAII_file):
        """HLA II alleles for which MixMHC2pred predictions are possible"""
        avail_alleles = []
        with open(path_to_HLAII_file) as f:
            for line in f:
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("L", "A")):
                        continue
                    line1 = line.split()[0]
                    if line1 is not None:
                        avail_alleles.append(line1)
        return avail_alleles

    def run(self, alleles, xmer_wt, xmer_mut):
        '''Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR
        '''
        self._initialise()
        tmp_prediction = intermediate_files.create_temp_file(prefix="mixmhc2pred", suffix=".txt")
        # prediction for peptides of length 13 to 18 based on Suppl Fig. 6 a in Racle, J., et al.
        # Robust prediction of HLA class II epitopes by deep motif deconvolution of immunopeptidomes.
        # Nat. Biotech. (2019).
        seqs = self.generate_nmers(xmer_wt=xmer_wt, xmer_mut=xmer_mut, lengths=[13, 14, 15, 16, 17, 18])
        tmp_fasta = intermediate_files.create_temp_fasta(seqs, prefix="tmp_sequence_")
        # try except statement to prevent stop of neofox for mps shorter < 13aa
        # TODO: this needs to be fixed, we could filter the list of nmers by length
        try:
            self.mixmhc2prediction(alleles, tmp_fasta, tmp_prediction)
        except:
            pass
        # TODO: also all of this try-catch needs to be fixed, in general the risk here is that they hide errors
        try:
            pred = self.read_mixmhcpred(tmp_prediction)
        except:
            pass
        try:
            pred_all = self.extract_best_per_pep(pred)
        except ValueError:
            pred_all = {}
        if len(pred_all) > 0:
            pred_best = self.extract_best_peptide_per_mutation(pred)
            self.best_peptide = self.add_best_epitope_info(pred_best, "Peptide")
            # TODO: improve how data is fetched so types are maintained
            self.best_rank = float(self.add_best_epitope_info(pred_best, "%Rank"))
            self.best_allele = self.add_best_epitope_info(pred_best, "BestAllele")
            self.all_peptides = "|".join(pred_all["Peptide"])
            self.all_ranks = "|".join(pred_all["%Rank"])
            self.all_alleles = "|".join(pred_all["BestAllele"])
            # prediction of for wt epitope that correspond to best epitope
            wt = self.extract_WT_for_best(xmer_wt=xmer_wt, xmer_mut=xmer_mut, best_mut_seq=self.best_peptide)
            wt_list = [wt]
            tmp_prediction = intermediate_files.create_temp_file(prefix="mixmhc2pred_wt_", suffix=".txt")
            tmp_fasta = intermediate_files.create_temp_fasta(wt_list, prefix="tmp_sequence_wt_")
            self.mixmhc2prediction([self.best_allele], tmp_fasta, tmp_prediction, wt=True)
            pred_wt = self.read_mixmhcpred(tmp_prediction)
            self.best_peptide_wt = self.extract_WT_info(pred_wt, "Peptide")
            # TODO: improve how data is fetched so types are maintained
            self.best_rank_wt = float(self.extract_WT_info(pred_wt, "%Rank"))
            
    def get_annotations(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(value=self.all_peptides, name="MixMHC2pred_all_peptides"),
            AnnotationFactory.build_annotation(value=self.all_ranks, name="MixMHC2pred_all_ranks"),
            AnnotationFactory.build_annotation(value=self.all_alleles, name="MixMHC2pred_all_alleles"),
            AnnotationFactory.build_annotation(value=self.best_peptide, name="MixMHC2pred_best_peptide"),
            AnnotationFactory.build_annotation(value=self.best_rank, name="MixMHC2pred_best_rank"),
            AnnotationFactory.build_annotation(value=self.best_allele, name="MixMHC2pred_best_allele"),
            AnnotationFactory.build_annotation(value=self.best_peptide_wt, name="MixMHC2pred_best_peptide_wt"),
            AnnotationFactory.build_annotation(value=self.best_rank_wt, name="MixMHC2pred_best_rank_wt"),
            AnnotationFactory.build_annotation(
                value=self.best_rank - self.best_rank_wt, name="MixMHC2pred_difference_rank_mut_wt")
        ]
