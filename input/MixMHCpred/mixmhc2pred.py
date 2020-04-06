#!/usr/bin/env python

import subprocess
import sys
import tempfile
from input.helpers import data_import


class MixMHC2pred:
    def __init__(self):
        self.all_peptides = "NA"
        self.all_ranks = "NA"
        self.all_alleles = "NA"
        self.best_peptide = "NA"
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
            for i,aa in enumerate(xmer_mut):
                if aa != xmer_wt[i]:
                    p1 = i + 1
        else:
            p1 = 0
            # in case sequences do not have same length
            for a1, a2 in zip(xmer_wt, xmer_mut):
                if a1 == a2:
                    p1 += 1
        return str(p1)

    def generate_nmers(self, props, list_lengths, mut = True ):
        ''' generates peptides covering mutation of all lengths that are provided. Returns peptides as list
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
            for s,e in zip(starts, ends):
                list_peptides.append(long_seq[s:e])
        list_peptides_fil = []
        [list_peptides_fil.append(x) for x in list_peptides if not x == ""]
        return list_peptides_fil


    def generate_fasta(self, seqs, tmpfile):
        ''' Writes seqs given in seqs list into fasta file
        '''
        counter = 0
        with open(tmpfile,"w") as f:
            for seq in seqs:
                id = "".join([">seq", str(counter)])
                f.write(id + "\n")
                f.write(seq + "\n")
                counter += 1

    def get_hla_allels(self, props, hla_patient_dict):
        ''' returns hla allele of patients given in hla_file
        '''
        if "patient.id" in props:
            patientid = props["patient.id"]
        else:
            try:
                patientid = props["patient"]
            except KeyError:
                patientid = props["patient.x"]
        return hla_patient_dict[patientid]

    def prepare_dq_dp(self, list_alleles, avail_hlaII):
        ''' returns patient DQ/DP alleles that are relevant for prediction
        '''
        list_alleles_pairs = ["__".join([p1,p2]) for p1 in list_alleles for p2 in list_alleles if p1 != p2]
        list_alleles_triplets = ["__".join([p1,p2,p3]) for p1 in list_alleles for p2 in list_alleles for p3 in list_alleles if p1 != p2 and p1 != p3 and p2 != p3]
        list_alleles_all = list_alleles_pairs + list_alleles_triplets
        alleles4pred = [allele for allele in list_alleles_all if allele in avail_hlaII]
        return(alleles4pred)


    def hlaIIallels2prediction(self, hla_alleles, avail_hlaII):
        ''' prepares list of hla alleles for prediction
        '''
        allels_for_prediction = []
        alleles_dq = []
        alleles_dp = []
        #print hla_alleles
        for allele in hla_alleles:
            #print allele
            allele = allele.replace("*", "_").replace(":", "_").replace("HLA-", "")
            if allele.startswith("DR"):
                if allele in avail_hlaII:
                    allels_for_prediction.append(allele)
            elif allele.startswith("DP"):
                alleles_dp.append(allele)
            elif allele.startswith("DQ"):
                alleles_dq.append(allele)
        alleles_dp4pred = self.prepare_dq_dp(alleles_dp, avail_hlaII)
        alleles_dq4pred = self.prepare_dq_dp(alleles_dq, avail_hlaII)
        allels_for_prediction = allels_for_prediction + alleles_dq4pred + alleles_dp4pred
        hla_allele = " ".join(allels_for_prediction)
        #print hla_allele
        return hla_allele

    def mixmhc2prediction(self, hla_alleles, tmpfasta, outtmp, avail_hlaII, wt = False):
        ''' Performs MixMHC2pred prediction for desired hla allele and writes result to temporary file.
        '''
        if not wt:
            hla_allele = self.hlaIIallels2prediction(hla_alleles, avail_hlaII)
        elif wt:
            # use best allele from mutated seq prediction
            hla_allele = hla_alleles[0]
        #
        cmd = "MixMHC2pred -a " + hla_allele + " -i " + tmpfasta + " -o " + outtmp
        print(cmd, file=sys.stderr)
        p = subprocess.Popen(cmd.split(" "),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        p_return = p.communicate()

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
                    #print line
                    dat.append(line)
        return header, dat


    def extract_best_per_pep(self, pred_dat):
        '''extract info of best allele prediction for all potential ligands per muatation
        '''
        head = pred_dat[0]
        dat = pred_dat[1]
        peps = []
        alleles = []
        ranks = []
        pepcol = head.index("Peptide")
        #scorecol = head.index("Score_bestAllele")
        allelecol = head.index("BestAllele")
        rankcol = head.index("%Rank")
        min_value = -1000000000000000000
        for ii,i in enumerate(dat):
            col_of_interest = [i[pepcol], i[rankcol], i[allelecol]]
            # all potential peptides per mutation --> return ditionary
            peps.append(i[pepcol])
            ranks.append(i[rankcol])
            alleles.append(i[allelecol])
        return {"Peptide": peps, "BestAllele": alleles, "%Rank": ranks }

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
        allelecol = head.index("BestAllele")
        rankcol = head.index("%Rank")
        min_value = 1000000000000000000
        for ii,i in enumerate(dat):
            col_of_interest = [str(i[pepcol]), str(i[rankcol]), str(i[allelecol])]
            # best ligand per mutation
            if float(i[rankcol]) < float(min_value):
                min_value = i[rankcol]
                min_pep = col_of_interest
        head_new = ["Peptide",  "%Rank", "BestAllele" ]
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
        wt_epi = xmer_wt[start:(start+l)]
        return(wt_epi)

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
            return str(float(wt_score) - float(mut_score))
        except ValueError:
            return "NA"

    def import_available_HLAII_alleles(self, path_to_HLAII_file):
        '''HLA II alleles for which MixMHC2pred predictions are possible
        '''
        avail_alleles = []
        with open(path_to_HLAII_file) as f:
            for line in f:
                line = line.rstrip().lstrip()
                if line:
                    if line.startswith(("L", "A")):
                        continue
                    line1 = line.split()[0]
                    #print line
                    if line1 is not None:
                        avail_alleles.append(line1)
        return avail_alleles


    def main(self, props_dict, dict_patient_hlaII, list_avail_hlaII):
        '''Wrapper for MHC binding prediction, extraction of best epitope and check if mutation is directed to TCR
        '''
        #list_avail_hlaII = self.import_available_HLAII_alleles(path_to_HLAII_file)
        tmp_fasta_file = tempfile.NamedTemporaryFile(prefix ="tmp_sequence_", suffix = ".fasta", delete = False)
        tmp_fasta = tmp_fasta_file.name
        tmp_prediction_file = tempfile.NamedTemporaryFile(prefix ="mixmhc2pred", suffix = ".txt", delete = False)
        tmp_prediction = tmp_prediction_file.name
        # prediction for peptides of length 13 to 18 based on Suppl Fig. 6 a in Racle, J., et al. Robust prediction of HLA class II epitopes by deep motif deconvolution of immunopeptidomes. Nat. Biotech. (2019).
        seqs = self.generate_nmers(props_dict, [13,14,15,16,17,18])
        self.generate_fasta(seqs, tmp_fasta)
        alleles = self.get_hla_allels(props_dict, dict_patient_hlaII)
        self.mixmhc2prediction(alleles, tmp_fasta, tmp_prediction, list_avail_hlaII)
        pred = self.read_mixmhcpred(tmp_prediction)
        try:
            pred_all = self.extract_best_per_pep(pred)
        except ValueError:
            pred_all ={}
        if len(pred_all) > 0:
            pred_best = self.extract_best_peptide_per_mutation(pred)
            self.best_peptide = self.add_best_epitope_info(pred_best, "Peptide")
            self.best_rank = self.add_best_epitope_info(pred_best, "%Rank")
            self.best_allele = self.add_best_epitope_info(pred_best, "BestAllele")
            #print "best allele " + self.best_allele
            self.all_peptides = "|".join(pred_all["Peptide"])
            self.all_ranks = "|".join(pred_all["%Rank"])
            self.all_alleles = "|".join(pred_all["BestAllele"])
            # prediction of for wt epitope that correspond to best epitope
            wt = self.extract_WT_for_best(props_dict, self.best_peptide)
            wt_list = [wt]
            tmp_fasta_file = tempfile.NamedTemporaryFile(prefix ="tmp_sequence_wt_", suffix = ".fasta", delete = False)
            tmp_fasta = tmp_fasta_file.name
            tmp_prediction_file = tempfile.NamedTemporaryFile(prefix ="mixmhc2pred_wt_", suffix = ".txt", delete = False)
            tmp_prediction = tmp_prediction_file.name
            self.generate_fasta(wt_list, tmp_fasta)
            self.mixmhc2prediction([self.best_allele], tmp_fasta, tmp_prediction, list_avail_hlaII, wt = True)
            pred_wt = self.read_mixmhcpred(tmp_prediction)
            #print >> sys.stderr, pred_wt
            self.best_peptide_wt = self.extract_WT_info(pred_wt, "Peptide")
            #score_wt_of_interest = "_".join(["Score",self.best_allele])
            #rank_wt_of_interest = "_".join(["%Rank",self.best_allele])
            #self.best_score_wt = self.extract_WT_info(pred_wt, score_wt_of_interest)
            self.best_rank_wt = self.extract_WT_info(pred_wt, "%Rank")
            # difference in scores between mut and wt
            self.difference_score_mut_wt = self.difference_score(self.best_rank_wt, self.best_rank)




if __name__ == '__main__':

    from input import predict_all_epitopes, epitope

    # alleles available in MixMHC2pred
    path_to_HLAII_file = "/projects/SUMMIT/WP1.2/input/development/MixMHCpred/Alleles_list_pred2.txt"
    list_avail_hlaII = MixMHC2pred().import_available_HLAII_alleles(path_to_HLAII_file)
    print(list_avail_hlaII)

    # test with ott data set
    #file = "/projects/CM01_iVAC/immunogenicity_prediction/3rd_party_solutions/MHC_prediction_netmhcpan4/testdat_ott.txt"
    #hla_file ="/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/ott/icam_ott/alleles.csv"
    file = "/projects/SUMMIT/WP1.2/input/development/netmhcIIpan/PtCU9061.test.txt"
    hla_file = "/projects/SUMMIT/WP1.2/Literature_Cohorts/data_analysis/cohorts/rizvi/icam_rizvi/20190819_alleles_extended.csv"
    # test inest data set
    #file = "/flash/projects/WP3/AnFranziska/AnFranziska/head_seqs.txt"
    #hla_file = "/flash/projects/WP3/AnFranziska/AnFranziska/alleles.csv"
    dat = data_import.import_dat_icam(file, False)
    if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in dat[0]:
        dat = data_import.change_col_names(dat)
    if "patient.id" not in dat[0]:
        try:
            patient = file.split("/")[-3]
            if "Pt" not in patient:
                patient = file.split("/")[-1].split(".")[0]
        except IndexError:
            patient = file.split("/")[-1].split(".")[0]
        dat[0].append("patient.id")
        for ii,i in enumerate(dat[1]):
            dat[1][ii].append(str(patient))
    # available MHC alleles
    set_available_mhc = predict_all_epitopes.Bunchepitopes().add_available_hla_alleles()
    # hla allele of patients
    patient_hlaI = predict_all_epitopes.Bunchepitopes().add_patient_hla_I_allels(hla_file)
    patient_hlaII = predict_all_epitopes.Bunchepitopes().add_patient_hla_II_allels(hla_file)

    print(patient_hlaI)
    print(patient_hlaII)

    for ii,i in enumerate(dat[1]):
        if ii < 10:
            print(ii)
            dict_epi = epitope.Epitope()
            dict_epi.init_properties(dat[0], dat[1][ii])
            prediction = MixMHC2pred()
            #print ii
            #print dict_epi.properties

            prediction.main(dict_epi.properties,  patient_hlaII, list_avail_hlaII)
            attrs = vars(prediction)
            print(attrs)


        #def wrapper(func, *args, **kwargs):
        #    def wrapped():
        #        return func(*args, **kwargs)
        #    return wrapped
        #wrapped = wrapper(prediction.main, dict_epi.properties, set_available_mhc, patient_hlaI)
        #prediction.generate_fasta(dict_epi.properties)
        #wrapped = wrapper(prediction.mhc_prediction, set_available_mhc, patient_hlaI)
        #print timeit.timeit(wrapped, number=1)
        #print timeit.timeit(wrapped, number=3)
        #print prediction.mhc_score
        #print prediction.epitope
        #print prediction.allele
        #print prediction.directed_to_TCR
