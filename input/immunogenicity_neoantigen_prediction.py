#!/usr/bin/env python

from logzero import logger

import input.aa_index.aa_index as aa_index
from input import MHC_I, MHC_II
from input.annotation_resources.nmer_frequency.nmer_frequency import AminoacidFrequency, FourmerFrequency
from input.epitope import Epitope
from input.annotation_resources.gtex.gtex import GTEx
from input.helpers import data_import
from input.helpers.properties_manager import PATIENT_ID
from input.helpers.runner import Runner
from input.new_features.conservation_scores import ProveanAnnotator
from input.references import ReferenceFolder, DependenciesConfiguration
from input.annotation_resources.uniprot.uniprot import Uniprot


class ImmunogenicityNeoantigenPredictionToolbox:

    def __init__(self, icam_file, patient_id, patients_file):

        self.patient_id = patient_id

        self.references = ReferenceFolder()
        self.configuration = DependenciesConfiguration()
        self.runner = Runner()
        self.gtex = GTEx()
        self.uniprot = Uniprot(self.references.uniprot)
        self.hla_available_alleles = self.references.load_available_hla_alleles(mhc=MHC_I)
        self.hlaII_available_alleles = self.references.load_available_hla_alleles(mhc=MHC_II)

        # TODO: remove this empty initialisations
        self.Allepit = {}

        # import epitope data
        self.header, self.rows = data_import.import_dat_icam(icam_file)
        patients = data_import.import_patients_data(patients_file)

        # TODO: remove once we are loading data into models
        if "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)" in self.header:
            self.header, self.rows = data_import.change_col_names(header=self.header, data=self.rows)

        # adds patient to the table
        self.header.append(PATIENT_ID)
        logger.debug(self.patient_id)

        # TODO: when we are moving away from the icam table we will have a patient id for each neoantigen and
        # TODO: we will be able to pass the whole list of patients below
        self.tissue = list(filter(lambda x: x.identifier == self.patient_id, patients))[0].tissue
        for row in self.rows:
            row.append(str(self.patient_id))

        self.aa_frequency = AminoacidFrequency()
        self.fourmer_frequency = FourmerFrequency()
        self.aa_index1_dict = aa_index.parse_aaindex1(self.references.aaindex1)
        self.aa_index2_dict = aa_index.parse_aaindex2(self.references.aaindex2)

        self.patient_hla_I_allels = {p.identifier: p.mhc_i_alleles for p in patients}
        self.patient_hla_II_allels = {p.identifier: p.mhc_i_i_alleles for p in patients}
        self.tumour_content = {p.identifier: p.estimated_tumor_content for p in patients}
        self.rna_avail = {p.identifier: p.is_rna_available for p in patients}
        self.provean_annotator = ProveanAnnotator(provean_file=self.references.prov_scores_mapped3,
                                                  header_epitopes=self.header, epitopes=self.rows)

    def write_to_file_sorted(self, d, header):
        """Transforms dictionary (property --> epitopes). To one unit (epitope) corresponding values are concentrated in one list
        and printed ';' separated."""
        features_names = []
        for key in d:
            if key not in header:
                features_names.append(key)
        features_names.sort()
        header.extend(features_names)
        print("\t".join(header))
        for i in range(len(d["mutation"])):
            z = [str(d[col][i]) for col in header]
            print("\t".join(z))

    def run(self):
        """ Loads epitope data (if file has been not imported to R; colnames need to be changed), adds data to class that are needed to calculate,
        calls epitope class --> determination of epitope properties,
        write to txt file
        """
        # feature calculation for each epitope
        for row in self.rows:
            # dict for each epitope
            epitope = Epitope(
                runner=self.runner, references=self.references, configuration=self.configuration,
                provean_annotator=self.provean_annotator, gtex=self.gtex, uniprot=self.uniprot)
            features = epitope.main(
                self.header, row, self.aa_frequency,
                self.fourmer_frequency, self.aa_index1_dict, self.aa_index2_dict,
                self.hla_available_alleles, self.hlaII_available_alleles, self.patient_hla_I_allels,
                self.patient_hla_II_allels, self.tumour_content, self.rna_avail, self.patient_id, self.tissue)
            for key in features:
                if key not in self.Allepit:
                    # keys are are feautres; values: list of feature values associated with mutated peptide sequence
                    self.Allepit[key] = [features[key]]
                else:
                    self.Allepit[key].append(features[key])
        self.write_to_file_sorted(self.Allepit, self.header)
