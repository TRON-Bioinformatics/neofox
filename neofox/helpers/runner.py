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
import subprocess
import time

from logzero import logger

from neofox.exceptions import NeofoxCommandException


class Runner(object):

    def __init__(self, verbose=True):
        self.verbose = verbose

    def run_command(self,  cmd, print_log=True, **kwargs):
        if print_log and self.verbose:
            logger.info("Starting command: {}".format(" ".join(cmd)))
        start = time.time()
        process = subprocess.Popen(
            self._preprocess_command(cmd),
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            **kwargs
        )
        output, errors = process.communicate()
        return_code = process.returncode
        end = time.time()
        if print_log and self.verbose:
            logger.info("Elapsed time {} seconds".format(round(end - start, 3)))
        if return_code == 0:
            if print_log and self.verbose:
                logger.info("Finished command correctly!")
        else:
            logger.error("Finished command with return code {}".format(return_code))
            logger.error(self._decode(output))
            logger.error(self._decode(errors))
            raise NeofoxCommandException(
                "Error running command '{}'".format(" ".join(cmd))
            )
        return self._decode(output), self._decode(errors)

    @staticmethod
    def _preprocess_command(cmd):
        """
        This makes sure that any parameter containing white spaces is passed appropriately
        """
        return " ".join(cmd).split(" ")

    def _decode(self, data):
        return data.decode("utf8")
