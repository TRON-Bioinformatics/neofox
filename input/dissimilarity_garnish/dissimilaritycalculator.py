#!/usr/bin/env python

import os
import os.path
import tempfile

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
        outfile_file = tempfile.NamedTemporaryFile(prefix="tmp_prot_", suffix=".xml", delete=False)
        outfile = outfile_file.name
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

    def calculate_dissimilarity(self, props, fastafile, references, filter_binder=False):
        '''wrapper for dissimilarity calculation
        '''
        mhc_mut = props["best_affinity_epitope_netmhcpan4"]
        mhc_aff = props["best_affinity_netmhcpan4"]
        with open(fastafile, "w") as f:
            id = ">M_1"
            f.write(id + "\n")
            f.write(mhc_mut + "\n")
        dissim = self._calc_dissimilarity(fastafile, id, references)
        if filter_binder:
            if float(mhc_aff) < 500:
                sc = str(dissim) if dissim != "NA" else "0"
            else:
                sc = "NA"
        else:
            sc = str(dissim) if dissim != "NA" else "0"
        return sc
