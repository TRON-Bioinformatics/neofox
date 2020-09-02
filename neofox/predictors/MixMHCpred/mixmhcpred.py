#!/usr/bin/env python
from typing import List
from logzero import logger
from neofox.model.neoantigen import Annotation
from neofox.model.wrappers import AnnotationFactory
from neofox.predictors.MixMHCpred.abstract_mixmhcpred import AbstractMixMHCpred
from neofox.helpers import intermediate_files


class MixMHCpred(AbstractMixMHCpred):

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration
        self._initialise()

    def _initialise(self):
        self.all_peptides = None
        self.all_scores = None
        self.all_ranks = None
        self.all_alleles = None
        self.best_peptide = None
        self.best_score = None
        self.best_rank = None
        self.best_allele = None
        self.best_peptide_wt = None
        self.best_score_wt = None
        self.best_rank_wt = None

    def _mixmhcprediction(self, hla_alleles, tmpfasta, outtmp):
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

    def _extract_best_per_pep(self, pred_dat):
        '''extract info of best allele prediction for all potential ligands per muatation
        '''
        head = pred_dat[0]
        dat = pred_dat[1]
        peps = []
        scores = []
        alleles = []
        ranks = []
        result = {}
        try:
            pepcol = head.index("Peptide")
            scorecol = head.index("Score_bestAllele")
            allelecol = head.index("BestAllele")
            rankcol = head.index("%Rank_bestAllele")
            for entry in sorted(dat, key=lambda x: float(x[rankcol])):
                # all potential peptides per mutation --> return ditionary
                peps.append(entry[pepcol])
                scores.append(entry[scorecol])
                ranks.append(entry[rankcol])
                alleles.append(entry[allelecol])
            result = {"Peptide": peps, "Score_bestAllele": scores, "BestAllele": alleles, "%Rank_bestAllele": ranks}
        except ValueError:
            pass
        return result

    def _extract_best_peptide_per_mutation(self, pred_dat):
        '''extract best predicted ligand per mutation
        '''
        head = pred_dat[0]
        dat = pred_dat[1]
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

    def run(self, xmer_wt, xmer_mut, alleles):
        '''Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR
        '''
        self._initialise()
        tmp_prediction = intermediate_files.create_temp_file(prefix="mixmhcpred", suffix=".txt")
        seqs = self.generate_nmers(xmer_wt=xmer_wt, xmer_mut=xmer_mut, lengths=[8, 9, 10, 11])
        tmp_fasta = intermediate_files.create_temp_fasta(seqs, prefix="tmp_sequence_")
        self._mixmhcprediction(alleles, tmp_fasta, tmp_prediction)
        pred = self.read_mixmhcpred(tmp_prediction)
        pred_all = self._extract_best_per_pep(pred)
        if len(pred[1]) > 0:
            pred_best = self._extract_best_peptide_per_mutation(pred)
            self.best_peptide = self.add_best_epitope_info(pred_best, "Peptide")
            self.best_score = float(self.add_best_epitope_info(pred_best, "Score_bestAllele"))
            self.best_rank = self.add_best_epitope_info(pred_best, "%Rank_bestAllele")
            self.best_allele = self.add_best_epitope_info(pred_best, "BestAllele")
            self.all_peptides = "|".join(pred_all["Peptide"])
            self.all_scores = "|".join(pred_all["Score_bestAllele"])
            self.all_ranks = "|".join(pred_all["%Rank_bestAllele"])
            self.all_alleles = "|".join(pred_all["BestAllele"])
            # prediction of for wt epitope that correspond to best epitope
            wt = self.extract_WT_for_best(xmer_wt=xmer_wt, xmer_mut=xmer_mut, best_mut_seq=self.best_peptide)
            wt_list = [wt]
            tmp_prediction = intermediate_files.create_temp_file(prefix="mixmhcpred_wt_", suffix=".txt")
            tmp_fasta = intermediate_files.create_temp_fasta(wt_list, prefix="tmp_sequence_wt_")
            self._mixmhcprediction(alleles, tmp_fasta, tmp_prediction)
            pred_wt = self.read_mixmhcpred(tmp_prediction)
            logger.debug(pred_wt)
            self.best_peptide_wt = self.extract_WT_info(pred_wt, "Peptide")
            score_wt_of_interest = "_".join(["Score", self.best_allele])
            rank_wt_of_interest = "_".join(["%Rank", self.best_allele])
            self.best_score_wt = float(self.extract_WT_info(pred_wt, score_wt_of_interest))
            self.best_rank_wt = self.extract_WT_info(pred_wt, rank_wt_of_interest)

    def get_annotations(self) -> List[Annotation]:
        return [
            AnnotationFactory.build_annotation(value=self.all_peptides, name="MixMHCpred_all_peptides"),
            AnnotationFactory.build_annotation(value=self.all_scores, name="MixMHCpred_all_scores"),
            AnnotationFactory.build_annotation(value=self.all_ranks, name="MixMHCpred_all_ranks"),
            AnnotationFactory.build_annotation(value=self.all_alleles, name="MixMHCpred_all_alleles"),
            AnnotationFactory.build_annotation(value=self.best_peptide, name="MixMHCpred_best_peptide"),
            AnnotationFactory.build_annotation(value=self.best_score, name="MixMHCpred_best_score"),
            AnnotationFactory.build_annotation(value=self.best_rank, name="MixMHCpred_best_rank"),
            AnnotationFactory.build_annotation(value=self.best_allele, name="MixMHCpred_best_allele"),
            AnnotationFactory.build_annotation(value=self.best_peptide_wt, name="MixMHCpred_best_peptide_wt"),
            AnnotationFactory.build_annotation(value=self.best_score_wt, name="MixMHCpred_best_score_wt"),
            AnnotationFactory.build_annotation(value=self.best_rank_wt, name="MixMHCpred_best_rank_wt"),
            AnnotationFactory.build_annotation(
                value=self.best_score - self.best_score_wt, name="MixMHCpred_difference_score_mut_wt")
            ]
