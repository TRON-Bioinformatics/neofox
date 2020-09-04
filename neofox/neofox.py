#!/usr/bin/env python
import logging
import os
import time
from typing import List
import logzero
from logzero import logger
import dask
from multiprocessing.pool import ThreadPool
from neofox import NEOFOX_LOG_FILE_ENV
from neofox.annotation_resources.gtex.gtex import GTEx
from neofox.annotator import NeoantigenAnnotator
from neofox.exceptions import NeofoxConfigurationException
from neofox.model.neoantigen import NeoantigenAnnotations, Neoantigen, Patient
from neofox.references.references import ReferenceFolder
from neofox.annotation_resources.uniprot.uniprot import Uniprot


class NeoFox:

    def __init__(self, neoantigens: List[Neoantigen], patient_id: str, patients: List[Patient], work_folder=None,
                 output_prefix = None, num_cpus: int):


        # initialise logs
        if work_folder and os.path.exists(work_folder):
            logfile = "".join([work_folder, "/", output_prefix, "_neofox.log"])
        else:
            logfile = os.environ.get(NEOFOX_LOG_FILE_ENV)
        if logfile is not None:
            logzero.logfile(logfile)
        # TODO: this does not work
        logzero.loglevel(logging.DEBUG)
        logger.info("Loading data...")

        # initialise dask
        dask.config.set(pool=ThreadPool(num_cpus))

        if neoantigens is None or patients is None:
            raise NeofoxConfigurationException("Missing input data to run Neofox")
        self.neoantigens = neoantigens
        self.patients = {patient.identifier: patient for patient in patients}
        # TODO: avoid overriding patient id parameter
        for n in self.neoantigens:
            if n.patient_identifier is None:
                n.patient_identifier = patient_id

        references = ReferenceFolder()

        # resources with long loading times
        self.uniprot = Uniprot(references.uniprot)
        self.gtex = GTEx()

        logger.info("Data loaded")

    def get_annotations(self) -> List[NeoantigenAnnotations]:
        """
        Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        logger.info("Starting NeoFox annotations...")
        # feature calculation for each epitope
        delayed_annotations = []
        for neoantigen in self.neoantigens:
            patient = self.patients.get(neoantigen.patient_identifier)
            logger.debug("Neoantigen: {}".format(neoantigen.to_json(indent=3)))
            logger.debug("Patient: {}".format(patient.to_json(indent=3)))
            delayed_annotation = dask.delayed(
                NeoantigenAnnotator(uniprot=self.uniprot, gtex=self.gtex).get_annotation)(neoantigen, patient)
            delayed_annotations.append(delayed_annotation)

        start = time.time()
        annotations = dask.compute(*delayed_annotations)
        end = time.time()
        logger.info("Elapsed time for annotating {} neoantigens {} seconds".format(
            len(self.neoantigens), int(end - start)))
        return annotations
