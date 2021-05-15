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
import orjson as json
import subprocess
from neofox.helpers.runner import Runner
from neofox.references.references import DependenciesConfiguration


class BlastpRunner(object):

    def __init__(self, runner: Runner, configuration: DependenciesConfiguration, proteome_db: str):
        self.runner = runner
        self.configuration = configuration
        self.database = proteome_db

    def calculate_similarity_database(self, peptide, a=26) -> int:
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
            "100000000",
            "-num_alignments",
            "1"
        ]

        best_hit = self._run_blastp(cmd=cmd, peptide=peptide)
        similarity_score = None
        if best_hit is not None:
            similarity_score = best_hit.get("hsps")[0].get("score")
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

        best_hit = self._run_blastp(cmd=cmd, peptide=peptide)
        wt_peptide = None
        if best_hit is not None:
            wt_peptide = best_hit.get("hsps")[0].get("hseq")
        return wt_peptide

    def _run_blastp(self, cmd, peptide):
        with subprocess.Popen(('echo', peptide), stdout=subprocess.PIPE) as echo:
            output, errors = self.runner.run_command(cmd=cmd, stdin=echo.stdout)
            echo.wait()
            echo.kill()
        results = json.loads(output)
        hits = results.get("BlastOutput2")[0].get("report").get("results").get("search").get("hits")
        best_hit = None
        if hits is not None and len(hits) > 0:
            best_hit = hits[0]
        return best_hit
