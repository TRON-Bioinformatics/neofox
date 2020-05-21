#!/usr/bin/env python

import os
import os.path

from input.helpers import intermediate_files
from input.neoantigen_fitness.Aligner_modified import Aligner


class DissimilarityCalculator(object):

    def __init__(self, runner, configuration):
        """
        :type runner: input.helpers.runner.Runner
        :type configuration: input.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    def _calc_dissimilarity(self, fasta_file, n, references):
        '''
        This function determines the dissimilarity to self-proteome of epitopes as described in Richman et al
        '''
        outfile = intermediate_files.create_temp_file(prefix="tmp_prot_", suffix=".xml")
        self.runner.run_command(cmd=[
            self.configuration.blastp,
            "-gapopen", "11",
            "-gapextend", "1",
            "-outfmt", "5",
            "-query", fasta_file,
            "-out", outfile,
            "-db", os.path.join(references.proteome_db, "homo_sapiens.mod"),
            "-evalue", "100000000"])
        aligner = Aligner()
        # set a to 32 for dissimilarity
        aligner.readAllBlastAlignments(outfile)
        aligner.computeR(a=32)
        kk = int(n.split("_")[1])
        x = aligner.Ri.get(kk)
        x_dis = "NA"
        if x is not None:
            x_dis = 1 - x
        os.remove(fasta_file)
        os.remove(outfile)
        return x_dis

    def calculate_dissimilarity(self, mhc_mutation, mhc_affinity, fastafile, references, filter_binder=False):
        '''wrapper for dissimilarity calculation
        '''
        with open(fastafile, "w") as f:
            id = ">M_1"
            f.write(id + "\n")
            f.write(mhc_mutation + "\n")
        dissim = self._calc_dissimilarity(fastafile, id, references)
        if filter_binder:
            if float(mhc_affinity) < 500:
                sc = str(dissim) if dissim != "NA" else "0"
            else:
                sc = "NA"
        else:
            sc = str(dissim) if dissim != "NA" else "0"
        return sc
