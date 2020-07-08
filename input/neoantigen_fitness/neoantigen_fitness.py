#!/usr/bin/env python

import os
import os.path
from logzero import logger

from input.helpers import intermediate_files
from input.helpers.blastp_runner import BlastpRunner


class NeoantigenFitnessCalculator(BlastpRunner):

    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        super().__init__(runner, configuration)

    def _calc_pathogen_similarity(self, fasta_file, iedb):
        """
        This function determines the PATHOGENSIMILARITY of epitopes according to Balachandran et al. using a blast
        search against the IEDB pathogenepitope database
        """
        outfile = self.run_blastp(fasta_file=fasta_file, database=os.path.join(iedb, "iedb_blast_db"))
        similarity = self.parse_blastp_output(blastp_output_file=outfile)
        os.remove(outfile)
        return similarity

    def wrap_pathogen_similarity(self, mutation, iedb):
        fastafile = intermediate_files.create_temp_fasta(sequences=[mutation], prefix="tmpseq", comment_prefix='M_')
        try:
            pathsim = self._calc_pathogen_similarity(fastafile, iedb)
        except Exception as ex:
            # TODO: do we need this at all? it should not fail and if it fails we probably want to just stop execution
            logger.exception(ex)
            pathsim = 0
        os.remove(fastafile)
        logger.info("Peptide {} has a pathogen similarity of {}".format(mutation, pathsim))
        return str(pathsim)

    def calculate_amplitude_mhc(self, score_mutation, score_wild_type, apply_correction=False):
        """
        This function calculates the amplitude between mutated and wt epitope according to Balachandran et al.
        when affinity is used, use correction from Luksza et al. *1/(1+0.0003*aff_wt)
        """
        amplitude_mhc = "NA"
        try:
            candidate_amplitude_mhc = float(score_wild_type) / float(score_mutation)
            if apply_correction:  #nine_mer or affinity:
                amplitude_mhc = str(candidate_amplitude_mhc * (self._calculate_correction(score_wild_type)))
            else:
                amplitude_mhc = str(candidate_amplitude_mhc)
        except(ZeroDivisionError, ValueError) as e:
            pass
        return amplitude_mhc

    def _calculate_correction(self, score_wild_type):
        return 1 / (1 + 0.0003 * float(score_wild_type))

    def calculate_recognition_potential(
            self, amplitude, pathogen_similarity, mutation_in_anchor, mhc_affinity_mut=None):
        """
        This function calculates the recognition potential, defined by the product of amplitude and pathogensimiliarity of an epitope according to Balachandran et al.
        F_alpha = - max (A_i x R_i)

        Returns (A_i x R_i) value only for nonanchor mutation and epitopes of length 9; only considered by Balachandran
        """
        recognition_potential = "NA"
        try:
            candidate_recognition_potential = str(float(amplitude) * float(pathogen_similarity))
            if mhc_affinity_mut:
                if mutation_in_anchor == "0" and float(mhc_affinity_mut) < 500.0:
                    recognition_potential = candidate_recognition_potential
            else:
                if mutation_in_anchor == "0":
                    recognition_potential = candidate_recognition_potential
        except ValueError:
            pass
        return recognition_potential
