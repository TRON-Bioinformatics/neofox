from argparse import ArgumentParser

from logzero import logger

from input.immunogenicity_neoantigen_prediction import ImmunogenicityNeoantigenPredictionToolbox


def input_cli():
    parser = ArgumentParser(description='adds patient information given in sample file of a cohort to merged icam file')
    parser.add_argument('-i', '--icam-file', dest='icam_file', help='define iCaM file which should be annotated',
                        required=True)
    # TODO: once we support the input from the models this parameter will not be required
    parser.add_argument('-p', '--patient-id', dest='patient_id', help='the patient id for the iCaM file',
                        required=True)
    parser.add_argument('-d', '--patients-data', dest='patients_data',
                        help='file with data for patients with columns: identifier, estimated_tumor_content, '
                             'is_rna_available, mhc_i_alleles, mhc_i_i_alleles, tissue',
                        required=True)
    args = parser.parse_args()

    icam_file = args.icam_file
    patient_id = args.patient_id
    patients_data = args.patients_data

    logger.info("Loading data...")
    input = ImmunogenicityNeoantigenPredictionToolbox()
    logger.info("Starting INPuT...")
    input.run(icam_file, patient_id, patients_data)
    logger.info("Finished INPuT...")
