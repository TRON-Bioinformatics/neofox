from input import MHC_I, MHC_II
from input.exceptions import INPuTInputParametersException

PATIENT_ID3 = "patient.x"

PATIENT_ID2 = "patient"

PATIENT_ID = "patient.id"


def get_wild_type_and_mutations(properties, mhc):
    if mhc == MHC_I:
        mutation = properties["MHC_I_epitope_.best_prediction."]
        wild_type = properties["MHC_I_epitope_.WT."]
    elif mhc == MHC_II:
        mutation = properties["MHC_II_epitope_.best_prediction."]
        wild_type = properties["MHC_II_epitope_.WT."]
    else:
        raise INPuTInputParametersException("Bad MHC value: {}".format(mhc))
    return wild_type, mutation


def get_wild_type_and_mutation_from_netmhcpan4(properties, nine_mer=False):
    if nine_mer:
        mutation = properties["best_affinity_epitope_netmhcpan4_9mer"]
        wild_type = properties["best_epitope_netmhcpan4_9mer_WT"]
    else:
        mutation = properties["best_affinity_epitope_netmhcpan4"]
        wild_type = properties["best_affinity_epitope_netmhcpan4_WT"]
    return wild_type, mutation


def get_hla_allele(props, hla_patient_dict):
    ''' returns hla allele of patients given in hla_file
    '''
    patient_id = get_patient_id(props)
    return hla_patient_dict[patient_id]


def get_patient_id(props):
    if PATIENT_ID in props:
        patient_id = props[PATIENT_ID]
    elif PATIENT_ID2 in props:
        patient_id = props[PATIENT_ID2]
    else:
        patient_id = props[PATIENT_ID3]
    return patient_id