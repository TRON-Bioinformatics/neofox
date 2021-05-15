#
# Copyright (c) 2020-2030 Translational Oncology at the Medical Center of the Johannes Gutenberg-University Mainz gGmbH.
#
# This file is part of Neofox
# (see https://github.com/tron-bioinformatics/neofox).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.#
import os
from math import exp

import orjson as json
import subprocess

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

from neofox.helpers import intermediate_files

from neofox.helpers.runner import Runner
from neofox.published_features.neoantigen_fitness.aligner import Aligner
from neofox.references.references import DependenciesConfiguration


class BlastpRunner(object):

    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, proteome_db: str):
        self.runner = runner
        self.configuration = configuration
        self.database = proteome_db

    def calculate_similarity_database(self, peptide, a=26) -> float:
        """
        This function runs BLASTP on a given database and returns a score defining the similarity of the input sequence
        to best BLAST hit
        """
        cmd = [
            self.configuration.blastp,
            "-gapopen",
            "11",
            "-gapextend",
            "1",
            "-outfmt",
            "15",
            "-db",
            self.database,
            "-evalue",
            "100000000"
        ]

        hits = self._run_blastp(cmd=cmd, peptide=peptide)
        local_alignments = []
        for hit in hits:
            hsp = hit.get("hsps")[0]
            query = hsp.get("qseq")
            target = hsp.get("hseq")
            if "-" not in query and "-" not in target:
                al = BlastpRunner.align(query, target)
                if al and len(al) > 0:
                    local_alignments.append(al[0])
        similarity_score = self.computeR(alignments=local_alignments, a=a)
        return similarity_score

    def old_calculate_similarity_database(self, peptide, a=26) -> int:
        """
        This function runs BLASTP on a given database and returns a score defining the similarity of the input sequence
        to best BLAST hit
        """
        outfile = intermediate_files.create_temp_file(
            prefix="tmp_blastp_", suffix=".xml"
        )
        input_fasta = intermediate_files.create_temp_fasta(
            sequences=[peptide], prefix="tmp_dissimilarity_", comment_prefix="M_"
        )
        cmd = [
            self.configuration.blastp,
            "-gapopen",
            "11",
            "-gapextend",
            "1",
            "-outfmt",
            "5",
            "-out",
            outfile,
            "-query",
            input_fasta,
            "-db",
            self.database,
            "-evalue",
            "100000000"
        ]
        self.runner.run_command(cmd=cmd)
        os.remove(input_fasta)
        aligner = Aligner()
        aligner.readAllBlastAlignments(outfile)
        # TODO: return gene name related to wt peptide
        aligner.computeR(a=a)
        similarity_score = aligner.Ri.get(1, 0)
        return similarity_score

    def get_most_similar_wt_epitope(self, peptide):
        cmd = [
            self.configuration.blastp,
            "-gapopen",
            "11",
            "-gapextend",
            "1",
            "-outfmt",
            "15",
            "-db",
            self.database,
            "-evalue",
            "100000000",
            "-qcov_hsp_perc",
            "100",
            "-num_alignments",
            "1"
        ]

        hits = self._run_blastp(cmd=cmd, peptide=peptide)
        wt_peptide = None
        if hits is not None and len(hits) > 0:
            best_hit = hits[0]
            wt_peptide = best_hit.get("hsps")[0].get("hseq")
        return wt_peptide

    def _run_blastp(self, cmd, peptide):
        with subprocess.Popen(('echo', peptide), stdout=subprocess.PIPE) as echo:
            output, errors = self.runner.run_command(cmd=cmd, stdin=echo.stdout)
            echo.wait()
            echo.kill()
        results = json.loads(output)
        hits = results.get("BlastOutput2")[0].get("report").get("results").get("search").get("hits")
        return hits

    @staticmethod
    def align(seq1, seq2):
        """
        Smith-Waterman alignment with default parameters.
        """
        matrix = matlist.blosum62
        gap_open = -11
        gap_extend = -1
        aln = pairwise2.align.localds(
            seq1.upper(), seq2.upper(), matrix, gap_open, gap_extend
        )
        return aln

    @staticmethod
    def computeR(alignments, a=26, k=4.87) -> float:
        """
        Compute TCR-recognition probabilities for each neoantigen.
        """
        # energies of all bound states of neoantigen i
        bindingEnergies = [-k * (a - el[2]) for el in alignments]
        # partition function, over all bound states and an unbound state
        lZ = Aligner.logSum(bindingEnergies + [0])
        lGb = Aligner.logSum(bindingEnergies)
        R = exp(lGb - lZ)
        return R
