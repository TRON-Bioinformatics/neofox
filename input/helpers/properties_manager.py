from input import MHC_I, MHC_II
from input.exceptions import INPuTInputParametersException

PATIENT_ID = "patient.id"


def get_gene(properties):
    if "gene.x" in properties:
        gene = properties["gene.x"]
    else:
        gene = properties["gene"]
    return gene


def get_substitution(properties):
    return properties["substitution"]


def get_expression(properties):
    expression = properties["transcript_expression"]
    expression = float(expression) if expression != "NA" else "NA"
    return expression


def get_mutation_aminoacid(properties):
    return properties["MUT_AA"]


def get_scores_multiple_binding(properties, mhc):
    if mhc == MHC_I:
        mutation = properties["MB_score_top10_harmonic"]
        wild_type = properties["MB_score_WT_top10_harmonic"]
    elif mhc == MHC_II:
        mutation = properties["MB_score_MHCII_top10_harmonic"]
        wild_type = properties["MB_score_MHCII_top10_WT_harmonic"]
    else:
        raise INPuTInputParametersException("Bad MHC value: {}".format(mhc))
    mutation = float(mutation) if mutation != "NA" else "NA"
    wild_type = float(wild_type) if wild_type != "NA" else "NA"
    return wild_type, mutation


def get_scores_netmhcpan4_affinity(properties, mhc):
    if mhc == MHC_I:
        mutation = properties["best_affinity_netmhcpan4"]
        wild_type = properties["best_affinity_netmhcpan4_WT"]
    elif mhc == MHC_II:
        mutation = properties["best_affinity_netmhcIIpan"]
        wild_type = properties["best_affinity_netmhcIIpan_WT"]
    else:
        raise INPuTInputParametersException("Bad MHC value: {}".format(mhc))
    mutation = float(mutation) if mutation != "NA" else "NA"
    wild_type = float(wild_type) if wild_type != "NA" else "NA"
    return wild_type, mutation


def get_scores_netmhcpan4_affinity_9mer(properties):
    mutation = properties["best_affinity_netmhcpan4_9mer"]
    wild_type = properties["best_affinity_netmhcpan4_9mer_WT"]
    mutation = float(mutation) if mutation != "NA" else "NA"
    wild_type = float(wild_type) if wild_type != "NA" else "NA"
    return wild_type, mutation


def get_scores_netmhcpan4_ranks(properties, mhc):
    if mhc == MHC_I:
        mutation = properties["best%Rank_netmhcpan4"]
        wild_type = properties["best%Rank_netmhcpan4_WT"]
    elif mhc == MHC_II:
        mutation = properties["best%Rank_netmhcIIpan"]
        wild_type = properties["best%Rank_netmhcIIpan_WT"]
    else:
        raise INPuTInputParametersException("Bad MHC value: {}".format(mhc))
    mutation = float(mutation) if mutation != "NA" else "NA"
    wild_type = float(wild_type) if wild_type != "NA" else "NA"
    return wild_type, mutation


def get_netmhcpan4_epitopes(properties, nine_mer=False):
    if nine_mer:
        mutation = properties["best_affinity_epitope_netmhcpan4_9mer"]
        wild_type = properties["best_affinity_epitope_netmhcpan4_9mer_WT"]
    else:
        mutation = properties["best_affinity_epitope_netmhcpan4"]
        wild_type = properties["best_affinity_epitope_netmhcpan4_WT"]
    return wild_type, mutation


def get_netmhcpan4_epitopes_rank(properties):
    mutation = properties["best_epitope_netmhcpan4"]
    wild_type = properties["best_epitope_netmhcpan4_WT"]
    return wild_type, mutation


def get_netmhciipan_epitopes(properties, affinity=False):
    if affinity:
        mutation = properties["best_affinity_epitope_netmhcIIpan"]
        wild_type = properties["best_affinity_epitope_netmhcIIpan_WT"]
    else:
        mutation = properties["best_epitope_netmhcIIpan"]
        wild_type = properties["best_epitope_netmhcIIpan_WT"]
    return wild_type, mutation


def get_patient_id(props):
    return props.get(PATIENT_ID)


def get_wt_mut_aa(substitution, mut_or_wt):
    """Returns wt and mut aa.
    """
    amino_acid = "NA"
    try:
        if mut_or_wt == "mut":
            amino_acid = substitution[-1]
        elif mut_or_wt == "wt":
            amino_acid = substitution[0]
    except ValueError:
        pass
    return amino_acid

