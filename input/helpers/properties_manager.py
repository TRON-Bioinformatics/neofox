from input import MHC_I, MHC_II
from input.exceptions import INPuTInputParametersException


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