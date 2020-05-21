#!/usr/bin/env python

import os
import os.path
import tempfile

from input.neoantigen_fitness.Aligner_modified import Aligner


class NeoantigenFitnessCalculator(object):

    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    def _calc_pathogensimilarity(self, fasta_file, n, iedb):
        '''
        This function determines the PATHOGENSIMILARITY of epitopes according to Balachandran et al. using a blast search against the IEDB pathogenepitope database
        '''
        outfile_file = tempfile.NamedTemporaryFile(prefix="tmp_iedb_", suffix=".xml", delete=False)
        outfile = outfile_file.name
        self.runner.run_command(cmd=[
            self.configuration.blastp,
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

    def wrap_pathogensimilarity(self, mutation, fastafile, iedb):
        with open(fastafile, "w") as f:
            id = ">M_1"
            f.write(id + "\n")
            f.write(mutation + "\n")
        try:
            pathsim = self._calc_pathogensimilarity(fastafile, id, iedb)
        except:
            pathsim = "NA"
        return str(pathsim) if pathsim != "NA" else "0"

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
                if mutation_in_anchor == "0" and mhc_affinity_mut < 500:
                    recognition_potential = candidate_recognition_potential
            else:
                if mutation_in_anchor == "0":
                    recognition_potential = candidate_recognition_potential
        except ValueError:
            pass
        return recognition_potential
