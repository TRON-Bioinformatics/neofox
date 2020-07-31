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
    input = ImmunogenicityNeoantigenPredictionToolbox(
        icam_file=icam_file, patients_file=patients_data, patient_id=patient_id)
    logger.info("Starting INPuT...")
    annotations, header = input.get_annotations()
    logger.info("Writing results...")
    write_to_file_sorted(annotations, header)
    logger.info("Finished INPuT...")


def write_to_file_sorted(annotations, header):
    """Transforms dictionary (property --> epitopes). To one unit (epitope) corresponding values are concentrated in one list
    and printed ';' separated."""

    transformed_annotations = {}
    for neoantigen in annotations:
        for key in neoantigen:
            if key not in transformed_annotations:
                # keys are are feautres; values: list of feature values associated with mutated peptide sequence
                transformed_annotations[key] = [neoantigen[key]]
            else:
                transformed_annotations[key].append(neoantigen[key])

    features_names = []
    for key in transformed_annotations:
        if key not in header:
            features_names.append(key)
    features_names.sort()
    header.extend(features_names)
    print("\t".join(header))
    for i in range(len(transformed_annotations["mutation"])):
        z = [str(transformed_annotations[col][i]) for col in header]
        print("\t".join(z))
