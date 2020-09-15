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


class BlastpRunner(object):

    def __init__(self, runner, configuration):
        """
        :type runner: neofox.helpers.runner.Runner
        :type configuration: neofox.references.DependenciesConfiguration
        """
        self.runner = runner
        self.configuration = configuration

    def run_blastp(self, fasta_file, database):
        '''
        This function runs BLASTP on a given database
        '''
        outfile = intermediate_files.create_temp_file(prefix="tmp_blastp_", suffix=".xml")
        self.runner.run_command(cmd=[
            self.configuration.blastp,
            "-gapopen", "11",
            "-gapextend", "1",
            "-outfmt", "5",
            "-query", fasta_file,
            "-out", outfile,
            "-db", database,
            "-evalue", "100000000"])
        return outfile

    def parse_blastp_output(self, blastp_output_file, **kwargs) -> int:
        aligner = Aligner()
        # set a to 32 for dissimilarity
        aligner.readAllBlastAlignments(blastp_output_file)
        aligner.computeR(**kwargs)
        return aligner.Ri.get(1, 0)  # NOTE: returns 0 when not present