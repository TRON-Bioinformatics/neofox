#!/usr/bin/env python

import os
import os.path
import tempfile

from input import MHC_I, MHC_II
from input.neoantigen_fitness.Aligner_modified import Aligner


class NeoantigenFitnessCalculator(object):

    def __init__(self, runner):
        """
        :type runner: input.helpers.runner.Runner
        """
        self.runner = runner

    def _calc_pathogensimilarity(self, fasta_file, n, iedb):
        '''
        This function determines the PATHOGENSIMILARITY of epitopes according to Balachandran et al. using a blast search against the IEDB pathogenepitope database
        '''
        outfile_file = tempfile.NamedTemporaryFile(prefix="tmp_iedb_", suffix=".xml", delete=False)
        outfile = outfile_file.name
        self.runner.run_command(cmd=[
            "/code/ncbi-blast/2.8.1+/bin/blastp",
            "-gapopen", "11",
            "-gapextend", "1",
            "-outfmt", "5",
            "-query", fasta_file,
            "-out", outfile,
            "-db", os.path.join(iedb, "iedb_blast_db"),
            "-evalue", "100000000"])
        a = Aligner()
        a.readAllBlastAlignments(outfile)
        a.computeR()
        kk = int(n.split("_")[1])
        x = a.Ri.get(kk)
        os.remove(fasta_file)
        os.remove(outfile)
        return x if x is not None else "NA"


    def wrap_pathogensimilarity(self, props, mhc, fastafile, iedb, affinity=False, nine_mer=False):
        if mhc == MHC_I:
            if nine_mer:
                mhc_mut = props["best_affinity_epitope_netmhcpan4_9mer"]
            elif affinity:
                mhc_mut = props["best_affinity_epitope_netmhcpan4"]
            else:
                mhc_mut = props["MHC_I_epitope_.best_prediction."]
        elif mhc == MHC_II:
            mhc_mut = props["MHC_II_epitope_.best_prediction."]
        with open(fastafile, "w") as f:
            id = ">M_1"
            f.write(id + "\n")
            f.write(mhc_mut + "\n")
        pathsim = self._calc_pathogensimilarity(fastafile, id, iedb)
        return str(pathsim) if pathsim != "NA" else "0"


    def calculate_amplitude_mhc(self, props, mhc, multiple_binding=False, affinity=False, netmhcscore=False, nine_mer=False):
        '''
        This function calculates the amplitude between mutated and wt epitope according to Balachandran et al.
        when affinity is used, use correction from Luksza et al. *1/(1+0.0003*aff_wt)
        '''
        if mhc == MHC_I:
            if multiple_binding:
                score_mutation = props["MB_score_top10_harmonic"].replace(",", ".")
                score_wild_type = props["MB_score_WT_top10_harmonic"].replace(",", ".")
            elif affinity:
                score_mutation = props["best_affinity_netmhcpan4"].replace(",", ".")
                score_wild_type = props["best_affinity_netmhcpan4_WT"].replace(",", ".")
            elif netmhcscore:
                score_mutation = props["best%Rank_netmhcpan4"].replace(",", ".")
                score_wild_type = props["best%Rank_netmhcpan4_WT"].replace(",", ".")
            elif nine_mer:
                score_mutation = props["best_affinity_netmhcpan4_9mer"].replace(",", ".")
                score_wild_type = props["best_affinity_netmhcpan4_9mer_WT"].replace(",", ".")
            else:
                score_mutation = props["MHC_I_score_.best_prediction."].replace(",", ".")
                score_wild_type = props["MHC_I_score_.WT."].replace(",", ".")
        elif mhc == MHC_II:
            if multiple_binding:
                score_mutation = props["MB_score_MHCII_top10_harmonic"].replace(",", ".")
                score_wild_type = props["MB_score_MHCII_top10_WT_harmonic"].replace(",", ".")
            elif affinity:
                score_mutation = props["best_affinity_netmhcIIpan"].replace(",", ".")
                score_wild_type = props["best_affinity_netmhcIIpan_WT"].replace(",", ".")
            elif netmhcscore:
                score_mutation = props["best%Rank_netmhcIIpan"].replace(",", ".")
                score_wild_type = props["best%Rank_netmhcIIpan_WT"].replace(",", ".")
            else:
                score_mutation = props["MHC_II_score_.best_prediction."].replace(",", ".")
                score_wild_type = props["MHC_II_score_.WT."].replace(",", ".")

        amplitude_mhc = "NA"
        try:
            candidate_amplitude_mhc = float(score_wild_type) / float(score_mutation)
            if nine_mer or affinity:
                amplitude_mhc = str(candidate_amplitude_mhc * (self._calculate_correction(score_wild_type)))
            else:
                amplitude_mhc = str(candidate_amplitude_mhc)
        except(ZeroDivisionError, ValueError) as e:
            pass
        return amplitude_mhc


    def _calculate_correction(self, score_wild_type):
        return 1 / (1 + 0.0003 * float(score_wild_type))


    def calculate_recognition_potential(self, props, mhc, affinity=False, netmhcscore=False, nine_mer=False):
        '''
        This function calculates the recognition potential, defined by the product of amplitude and pathogensimiliarity of an epitope according to Balachandran et al.
        F_alpha = - max (A_i x R_i)

        Returns (A_i x R_i) value only for nonanchor mutation and epitopes of length 9; only considered by Balachandran
        '''
        if mhc == MHC_I:
            if affinity:
                amplitude = props["Amplitude_mhcI_affinity"]
                pathogen_similarity = props["Pathogensimiliarity_mhcI_affinity_nmers"]
            elif netmhcscore:
                amplitude = props["Amplitude_mhcI_rank_netmhcpan4"]
                pathogen_similarity = props["Pathogensimiliarity_mhcI"]
            elif nine_mer:
                amplitude = props["Amplitude_mhcI_affinity_9mer_netmhcpan4"]
                pathogen_similarity = props["Pathogensimiliarity_mhcI_9mer"]
            else:
                amplitude = props["Amplitude_mhcI"]
                pathogen_similarity = props["Pathogensimiliarity_mhcI"]
        elif mhc == MHC_II:
            amplitude = props["Amplitude_mhcII"]
            pathogen_similarity = props["Pathogensimiliarity_mhcII"]

        recognition_potential = "NA"
        try:
            # mutation_in_anchor = props["Mutation_in_anchor"]
            candidate_recognition_potential = str(float(amplitude) * float(pathogen_similarity))
            if nine_mer:
                mutation_in_anchor = props["Mutation_in_anchor_netmhcpan_9mer"]
                mhc_affinity_mut = float(props["best_affinity_netmhcpan4_9mer"])
                if mutation_in_anchor == "0" and mhc_affinity_mut < 500:
                    recognition_potential = candidate_recognition_potential
            else:
                mutation_in_anchor = props["Mutation_in_anchor_netmhcpan"]
                if mutation_in_anchor == "0":
                    recognition_potential = candidate_recognition_potential
        except ValueError:
            pass
        return recognition_potential
