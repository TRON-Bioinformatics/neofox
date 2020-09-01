from argparse import ArgumentParser
from logzero import logger

from neofox.model.conversion import ModelConverter
from neofox.neofox import NeoFox


def neofox_cli():
    parser = ArgumentParser(description='adds patient information given in sample file of a cohort to merged icam file')
    parser.add_argument('-i', '--icam-file', dest='icam_file', help='define iCaM file which should be annotated',
                        required=True)
    # TODO: once we support the neofox from the models this parameter will not be required
    parser.add_argument('-p', '--patient-id', dest='patient_id', help='the patient id for the iCaM file',
                        required=True)
    parser.add_argument('-d', '--patients-data', dest='patients_data',
                        help='file with data for patients with columns: identifier, estimated_tumor_content, '
                             'is_rna_available, mhc_i_alleles, mhc_i_i_alleles, tissue',
                        required=True)
    parser.add_argument('-o', '--output-file', dest='output_file', help='output file', required=True)
    args = parser.parse_args()

    icam_file = args.icam_file
    patient_id = args.patient_id
    patients_data = args.patients_data
    output_file = args.output_file

    neoantigens = ModelConverter.parse_icam_file(icam_file)
    patients = ModelConverter.parse_patients_file(patients_data)
    annotations = NeoFox(neoantigens=neoantigens, patients=patients, patient_id=patient_id).get_annotations()
    ModelConverter.annotations2short_wide_table(annotations, neoantigens).to_csv(output_file, sep='\t', index=False)
    logger.info("Finished NeoFox")
