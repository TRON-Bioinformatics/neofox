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
from math import exp, log
import orjson as json
import subprocess
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from neofox.helpers.runner import Runner
from neofox.references.references import DependenciesConfiguration


class BlastpRunner(object):

    INF = float("inf")

    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, database: str):
        self.runner = runner
        self.configuration = configuration
        self.database = database
        self.cache_homologous_epitopes = {}

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
        try:
            for hit in hits:
                hsp = hit.get("hsps")[0]
                query = hsp.get("qseq")
                target = hsp.get("hseq")
                if "-" not in query and "-" not in target:
                    al = BlastpRunner.align(query, target)
                    if al and len(al) > 0:
                        local_alignments.append(al[0])
            similarity_score = self.computeR(alignments=local_alignments, a=a)
        except SystemError:
            # NOTE: some rarer aminoacids substitutions may not be present in the BLOSUM matrix and thus cause this to
            # fail, an example is U>Y
            similarity_score = None
        return similarity_score

    def get_most_similar_wt_epitope(self, peptide):
        if peptide not in self.cache_homologous_epitopes:
            cmd = [
                self.configuration.blastp,
                "-outfmt",
                "15",
                "-db",
                self.database,
                "-evalue",
                "100000000",
                "-qcov_hsp_perc",
                "100",
                "-comp_based_stats F",
                "-num_alignments 1",
                "-ungapped"
            ]

            hits = self._run_blastp(cmd=cmd, peptide=peptide, print_log=True)
            wt_peptide = None
            if hits is not None and len(hits) > 0:
                best_hit = hits[0]
                wt_peptide = best_hit.get("hsps")[0].get("hseq")
            self.cache_homologous_epitopes[peptide] = wt_peptide
        else:
            wt_peptide = self.cache_homologous_epitopes.get(peptide)
        return wt_peptide

    def _run_blastp(self, cmd, peptide, print_log=True):
        with subprocess.Popen(('echo', peptide), stdout=subprocess.PIPE) as echo:
            output, errors = self.runner.run_command(cmd=cmd, stdin=echo.stdout, print_log=print_log)
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
        lZ = BlastpRunner.logSum(bindingEnergies + [0])
        lGb = BlastpRunner.logSum(bindingEnergies)
        R = exp(lGb - lZ)
        return R

    @staticmethod
    def logSum(v):
        """
        compute the logarithm of a sum of exponentials
        """
        if len(v) == 0:
            return -BlastpRunner.INF
        ma = max(v)
        if ma == -BlastpRunner.INF:
            return -BlastpRunner.INF
        return log(sum([exp(x - ma) for x in v])) + ma
