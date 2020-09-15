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
import sys
import time
from multiprocessing.pool import ThreadPool
from unittest import TestCase
import dask
from dask.distributed import Client
import numpy as np
from logzero import logger
from neofox.helpers.intermediate_files import create_temp_file

from neofox.model.neoantigen import Neoantigen, Patient

from neofox.tests.unit_tests.tools import get_random_neoantigen, get_random_patient

DATA = [(np.random.uniform(1, 100, 1000), np.random.uniform(1, 100, 1000)) for _ in range(10)]
NEOANTIGENS = [(get_random_neoantigen(), get_random_patient()) for _ in range(10)]


class TestDask(TestCase):

    @staticmethod
    def _lengthy_computation(d1, d2):
        logger.info("Zzzzzzz")
        time.sleep(10)
        logger.info("Waky waky!")
        return np.convolve(d1, d2)

    @staticmethod
    def _lengthy_computation_on_neoantigens(n: Neoantigen, p: Patient):
        assert isinstance(n, Neoantigen)
        assert isinstance(p, Patient)
        logger.info("Zzzzzzz")
        time.sleep(5)
        logger.info("Waky waky!")
        return n.to_json(indent=3) + p.identifier

    @staticmethod
    def _lengthy_computation_with_io(n: Neoantigen, p: Patient):
        assert isinstance(n, Neoantigen)
        assert isinstance(p, Patient)
        logger.info("Zzzzzzz")
        filename = create_temp_file()
        with open(filename, "w") as f:
            for _ in range(100000):
                f.write(n.to_json())
        logger.info("Waky waky!")
        return n.to_json(indent=3) + p.identifier

    # def test_delayed_with_thread_pool(self):
    #     dask.config.set(pool=ThreadPool(4))
    #     delayed = []
    #     for d1, d2 in DATA:
    #         delayed.append(dask.delayed(TestDask._lengthy_computation)(d1, d2))
    #
    #     result = dask.compute(*delayed)
    #     logger.info("Result: {}".format(result))
    #     self.assertTrue(result is not None)
    #     self.assertTrue(len(result) == 10)
    #
    # def test_delayed_with_thread_pool_on_neoantigens(self):
    #     dask.config.set(pool=ThreadPool(4))
    #     delayed = []
    #     for n, p in NEOANTIGENS:
    #         delayed.append(dask.delayed(TestDask._lengthy_computation_on_neoantigens)(n, p))
    #
    #     result = dask.compute(*delayed)
    #     logger.info("Result: {}".format(result))
    #     self.assertTrue(result is not None)
    #     self.assertTrue(len(result) == 10)

    # def test_futures_with_threads(self):
    #     client = Client(processes=False, n_workers=4)
    #     futures = []
    #     for d1, d2 in DATA:
    #         futures.append(client.submit(TestDask._lengthy_computation, d1, d2))
    #     result = client.gather(futures)
    #     logger.info("Result: {}".format(result))
    #     self.assertTrue(result is not None)
    #     self.assertTrue(len(result) == 10)
    #
    # def test_futures_with_threads_on_neoantigens(self):
    #     client = Client(processes=False, n_workers=4)
    #     futures = []
    #     for n, p in NEOANTIGENS:
    #         futures.append(client.submit(TestDask._lengthy_computation_on_neoantigens, n, p))
    #     result = client.gather(futures)
    #     logger.info("Result: {}".format(result))
    #     self.assertTrue(result is not None)
    #     self.assertTrue(len(result) == 10)
    #
    # def test_futures_with_processes(self):
    #     client = Client(processes=True, n_workers=4)
    #     futures = []
    #     for d1, d2 in DATA:
    #         futures.append(client.submit(TestDask._lengthy_computation, d1, d2))
    #     result = client.gather(futures)
    #     logger.info("Result: {}".format(result))
    #     self.assertTrue(result is not None)
    #     self.assertTrue(len(result) == 10)
    #
    # def test_futures_with_processes_on_neoantigens(self):
    #     client = Client(processes=True, n_workers=4)
    #     futures = []
    #     for n, p in NEOANTIGENS:
    #         futures.append(client.submit(TestDask._lengthy_computation_on_neoantigens, n, p))
    #     result = client.gather(futures)
    #     logger.info("Result: {}".format(result))
    #     self.assertTrue(result is not None)
    #     self.assertTrue(len(result) == 10)
    #
    # def test_futures_with_processes_with_io(self):
    #     client = Client(processes=True, n_workers=4, threads_per_worker=4)
    #     futures = []
    #     for n, p in NEOANTIGENS:
    #         futures.append(client.submit(TestDask._lengthy_computation_with_io, n, p))
    #     result = client.gather(futures)
    #     logger.info("Result: {}".format(result))
    #     self.assertTrue(result is not None)
    #     self.assertTrue(len(result) == 10)
