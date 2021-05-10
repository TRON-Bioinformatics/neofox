import pickle
import subprocess
import os
import pandas as pd
from Bio import SeqIO

from neofox import NEOFOX_HLA_DATABASE_ENV
from neofox.exceptions import NeofoxReferenceException
from neofox.helpers.runner import Runner
from neofox.references.references import (
    DependenciesConfigurationForInstaller,
    NETMHCPAN_AVAILABLE_ALLELES_FILE,
    NETMHC2PAN_AVAILABLE_ALLELES_FILE,
    IEDB_FOLDER,
    PROTEOME_DB_FOLDER,
    IEDB_BLAST_PREFIX,
    IEDB_FASTA,
    HOMO_SAPIENS_FASTA,
    PREFIX_HOMO_SAPIENS, HLA_DATABASE_AVAILABLE_ALLELES_FILE, HOMO_SAPIENS_PICKLE,
)
from logzero import logger


class NeofoxReferenceInstaller(object):
    def __init__(self, reference_folder, install_r_dependencies=False):
        self.config = DependenciesConfigurationForInstaller()
        self.runner = Runner()
        self.reference_folder = reference_folder
        self.install_r_dependencies = install_r_dependencies

    def install(self):
        # ensures the reference folder exists
        os.makedirs(self.reference_folder, exist_ok=True)

        self._set_netmhcpan_alleles()
        self._set_netmhc2pan_alleles()
        self._set_iedb()
        self._set_proteome()
        self._set_ipd_imgt_hla_database()
        if self.install_r_dependencies:
            self._install_r_dependencies()
        else:
            logger.warning("R dependencies will need to be installed manually")

    def _set_netmhcpan_alleles(self):
        # available MHC alleles netMHCpan
        # $NEOFOX_NETMHCPAN -listMHC | grep "HLA-" > "$NEOFOX_REFERENCE_FOLDER"/MHC_available.csv
        logger.info("Fetching available alleles from NetMHCpan")
        available_alleles_file = os.path.join(
            self.reference_folder, NETMHCPAN_AVAILABLE_ALLELES_FILE
        )
        cmd = '{netmhcpan} -listMHC | grep "HLA-" | grep -v "#" > {available_alleles_file}'.format(
            netmhcpan=self.config.net_mhc_pan,
            available_alleles_file=available_alleles_file,
        )
        self._run_command(cmd)

    def _set_netmhc2pan_alleles(self):
        # available MHCII alleles netMHCIIpan
        # $NEOFOX_NETMHC2PAN -list  > "$NEOFOX_REFERENCE_FOLDER"/avail_mhcII.txt
        logger.info("Fetching available alleles from NetMHC2pan")
        available_alleles_file = os.path.join(
            self.reference_folder, NETMHC2PAN_AVAILABLE_ALLELES_FILE
        )
        cmd = "{netmhc2pan} -list > {available_alleles_file}".format(
            netmhc2pan=self.config.net_mhc2_pan,
            available_alleles_file=available_alleles_file,
        )
        self._run_command(cmd)

    def _set_iedb(self):
        # build IEDB blast database
        # mkdir "$NEOFOX_REFERENCE_FOLDER"/iedb
        # wget "http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip" -O "$NEOFOX_REFERENCE_FOLDER"/iedb/Iedb.zip
        # unzip "$NEOFOX_REFERENCE_FOLDER"/iedb/Iedb.zip -d "$NEOFOX_REFERENCE_FOLDER"/iedb/
        # python $DIR/build_IEDB_db.py "$NEOFOX_REFERENCE_FOLDER"/iedb/tcell_full_v3.csv "$NEOFOX_REFERENCE_FOLDER"/iedb/IEDB.fasta
        # $NEOFOX_MAKEBLASTDB -in "$NEOFOX_REFERENCE_FOLDER"/iedb/IEDB.fasta -dbtype prot -parse_seqids -out "$NEOFOX_REFERENCE_FOLDER"/iedb/iedb

        logger.info("Configuring IEDB")

        os.makedirs(os.path.join(self.reference_folder, IEDB_FOLDER), exist_ok=True)

        # download IEDB
        iedb_zip = os.path.join(self.reference_folder, IEDB_FOLDER, "Iedb.zip")
        cmd = 'wget "http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip" -O {}'.format(
            iedb_zip
        )
        self._run_command(cmd)

        # unzip IEDB
        path_to_iedb_folder = os.path.join(self.reference_folder, IEDB_FOLDER)
        cmd = "unzip -o {iedb_zip} -d {iedb_folder}".format(
            iedb_zip=iedb_zip, iedb_folder=path_to_iedb_folder
        )
        self._run_command(cmd)

        # transforms IEDB into fasta
        iedb_tcell_csv = os.path.join(
            self.reference_folder, IEDB_FOLDER, "tcell_full_v3.csv"
        )
        iedb_fasta = os.path.join(self.reference_folder, IEDB_FOLDER, IEDB_FASTA)
        IedbFastaBuilder(iedb_tcell_csv, iedb_fasta).build_fasta()

        # run makeblastdb on iedb fasta
        cmd = "{makeblastdb} -in {iedb_fasta} -dbtype prot -out {iedb_folder}".format(
            makeblastdb=self.config.make_blastdb,
            iedb_fasta=iedb_fasta,
            iedb_folder=os.path.join(path_to_iedb_folder, IEDB_BLAST_PREFIX),
        )
        self._run_command(cmd)

    def _set_proteome(self):
        # human proteome database
        # mkdir "$NEOFOX_REFERENCE_FOLDER"/proteome_db
        # wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O "$NEOFOX_REFERENCE_FOLDER"/proteome_db/Homo_sapiens.fa.gz
        # gunzip "$NEOFOX_REFERENCE_FOLDER"/proteome_db/Homo_sapiens.fa.gz
        # $NEOFOX_MAKEBLASTDB -in "$NEOFOX_REFERENCE_FOLDER"/proteome_db/Homo_sapiens.fa -dbtype prot -parse_seqids -out "$NEOFOX_REFERENCE_FOLDER"/proteome_db/homo_sapiens

        logger.info("Configuring the proteome DB")

        os.makedirs(os.path.join(self.reference_folder, PROTEOME_DB_FOLDER))

        # download proteome
        proteome_compressed_file = os.path.join(
            self.reference_folder, PROTEOME_DB_FOLDER, "%s.gz" % HOMO_SAPIENS_FASTA
        )
        # ftp_url = "ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
        ftp_url = "ftp://ftp.ensembl.org/pub/grch37/release-101/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz"
        cmd = "wget {ftp_url} -O {proteome_file}".format(
            ftp_url=ftp_url, proteome_file=proteome_compressed_file
        )
        self._run_command(cmd)

        cmd = "gunzip -f {proteome_file}".format(proteome_file=proteome_compressed_file)
        self._run_command(cmd)

        proteome_file = os.path.join(
            self.reference_folder, PROTEOME_DB_FOLDER, HOMO_SAPIENS_FASTA
        )
        output_folder = os.path.join(
            self.reference_folder, PROTEOME_DB_FOLDER, PREFIX_HOMO_SAPIENS
        )
        cmd = "{makeblastdb} -in {proteome_file} -dbtype prot -parse_seqids -out {output_folder}".format(
            makeblastdb=self.config.make_blastdb,
            proteome_file=proteome_file,
            output_folder=output_folder,
        )
        self._run_command(cmd)

        # builds proteome in pickle for querying
        prepared_proteome = []
        for record in SeqIO.parse(proteome_file, "fasta"):
            prepared_proteome.append(str(record.seq))

        proteome_pickle = os.path.join(
            self.reference_folder, PROTEOME_DB_FOLDER, HOMO_SAPIENS_PICKLE
        )
        outfile = open(proteome_pickle, 'wb')
        pickle.dump("\n".join(prepared_proteome), outfile)
        outfile.close()

    def _set_ipd_imgt_hla_database(self):
        logger.info("Downloading the IPD-IMGT/HLA database")
        allele_list = os.path.join(self.reference_folder, HLA_DATABASE_AVAILABLE_ALLELES_FILE)
        url = os.environ.get(
            NEOFOX_HLA_DATABASE_ENV, "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.3430.txt")
        cmd = "wget {url} -O {allele_list}".format(url=url, allele_list=allele_list)
        self._run_command(cmd)

    def _install_r_dependencies(self):
        logger.info("Installing R dependencies...")
        cmd = "{rscript} {dependencies_file}".format(
            rscript=self.config.rscript,
            dependencies_file=os.path.join(
                os.path.abspath(os.path.dirname(__file__)), "install_r_dependencies.R"
            ),
        )
        self._run_command(cmd)

    def _run_command(self, cmd):
        logger.info(cmd)
        process = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        output, errors = process.communicate()
        if process.returncode != 0:
            logger.error(errors)
            raise NeofoxReferenceException("Error running command {}".format(cmd))
        logger.info("Success")


class IedbFastaBuilder:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def build_fasta(self):
        # read IEDB input file
        iedb = pd.read_csv(self.input_file, skiprows=1)

        # filter entries
        filtered_iedb = iedb[
            (
                iedb["Name"].isin(
                    [
                        "Homo sapiens",
                        "Homo sapiens (human)",
                        "Homo sapiens Caucasian",
                        "Homo sapiens Black",
                    ]
                )
            )
            & (iedb["Object Type"] == "Linear peptide")
            & (iedb["Process Type"] == "Occurrence of infectious disease")
            & (iedb["Qualitative Measure"] == "Positive")
            & (iedb["Class"] == "I")
        ]

        # sets values for identifiers and sequences
        filtered_iedb.loc[:, "seq"] = filtered_iedb.loc[:, "Description"].transform(
            lambda x: x.strip()
        )
        # build fasta header: 449|FL-160-2 protein - Trypanosoma cruzi|JH0823|Trypanosoma cruzi|5693
        # epitope id|Antigen Name|antigen_id|Organism Name|organism_id
        filtered_iedb.loc[:, "epitope_id"] = filtered_iedb.loc[
            :, "Epitope IRI"
        ].transform(lambda x: x.replace("http://www.iedb.org/epitope/", "", regex=True))
        filtered_iedb.loc[:, "antigen_id"] = filtered_iedb.loc[
            :, "Antigen IRI"
        ].transform(
            lambda x: x.replace(
                "http://www.ncbi.nlm.nih.gov/protein/", "", regex=True
            ).replace("https://ontology.iedb.org/ontology/", "", regex=True)
        )
        filtered_iedb.loc[:, "organism_id"] = filtered_iedb.loc[
            :, "Organism IRI"
        ].transform(
            lambda x: x.replace(
                "http://purl.obolibrary.org/obo/NCBITaxon_", "", regex=True
            )
        )
        filtered_iedb.loc[:, "fasta_header"] = filtered_iedb.apply(
            lambda row: ">{epitope_id}|{antigen_name}|{antigen_id}|{organism_name}|{organism_id}".format(
                epitope_id=str(row["epitope_id"]),
                antigen_name=row["Antigen Name"],
                antigen_id=str(row["antigen_id"]),
                organism_name=row["Organism Name"],
                organism_id=str(row["organism_id"]),
            ),
            axis=1,
        )
        filtered_iedb.drop_duplicates(subset="seq", keep="last", inplace=True)

        # writes output FASTA file
        with open(self.output_file, "w") as fasta:
            for index, row in filtered_iedb.iterrows():
                fasta.write(">{header}\n".format(header=str(row["fasta_header"])))
                fasta.write("{sequence}\n".format(sequence=str(row["seq"])))
