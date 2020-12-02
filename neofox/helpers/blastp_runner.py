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
from neofox.helpers import intermediate_files
from neofox.published_features.neoantigen_fitness.aligner import Aligner
import os


class BlastpRunner(object):
    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    def run_blastp(self, peptide, database, a=26) -> int:
        """
        This function runs BLASTP on a given database
        """
        input_fasta = intermediate_files.create_temp_fasta(
            sequences=[peptide], prefix="tmp_dissimilarity_", comment_prefix="M_"
        )
        outfile = intermediate_files.create_temp_file(
            prefix="tmp_blastp_", suffix=".xml"
        )
        self.runner.run_command(
            cmd=[
                self.configuration.blastp,
                "-gapopen",
                "11",
                "-gapextend",
                "1",
                "-outfmt",
                "5",
                "-query",
                input_fasta,
                "-out",
                outfile,
                "-db",
                database,
                "-evalue",
                "100000000",
            ]
        )

        score = self._parse_blastp_output(outfile, a=a)
        os.remove(outfile)
        os.remove(input_fasta)
        return score

    def _parse_blastp_output(self, blastp_output_file, a) -> int:
        aligner = Aligner()
        # set a to 32 for dissimilarity
        try:
            aligner.readAllBlastAlignments(blastp_output_file)
            aligner.computeR(a=a)
            result = aligner.Ri.get(1, 0)  # NOTE: returns 0 when not present
        except SystemError:
            # NOTE: some rarer aminoacids substitutions may not be present in the BLOSUM matrix and thus cause this to
            # fail, an example is U>Y
            result = None
        return result
