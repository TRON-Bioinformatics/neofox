import json
import pickle
import subprocess
import os
from shutil import copyfile
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
import hashlib
from neofox import NEOFOX_HLA_DATABASE_ENV
from neofox.exceptions import NeofoxReferenceException
from neofox.helpers.runner import Runner
from neofox.model.neoantigen import Resource
from neofox.references.references import (
    DependenciesConfigurationForInstaller,
    NETMHCPAN_AVAILABLE_ALLELES_FILE,
    NETMHC2PAN_AVAILABLE_ALLELES_FILE,
    IEDB_FOLDER,
    PROTEOME_DB_FOLDER,
    IEDB_FASTA_HOMO_SAPIENS,
    HOMO_SAPIENS_FASTA,
    PREFIX_HOMO_SAPIENS, HLA_DATABASE_AVAILABLE_ALLELES_FILE, HOMO_SAPIENS_PICKLE,
    NETMHCPAN_AVAILABLE_ALLELES_MICE_FILE, NETMHC2PAN_AVAILABLE_ALLELES_MICE_FILE, MUS_MUSCULUS_FASTA,
    PREFIX_MUS_MUSCULUS, MUS_MUSCULUS_PICKLE, IEDB_FASTA_MUS_MUSCULUS, IEDB_BLAST_PREFIX_HOMO_SAPIENS,
    IEDB_BLAST_PREFIX_MUS_MUSCULUS, H2_DATABASE_AVAILABLE_ALLELES_FILE, RESOURCES_VERSIONS,
)
from logzero import logger

IMGT_HLA_DB_URL = "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt"
MOUSE_PROTEOME_VERSION = "GRCm39"
HUMAN_PROTEOME_VERSION = "GRCh37"
MOUSE_PROTEOME = "ftp://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz"
HUMAN_PROTEOME = "ftp://ftp.ensembl.org/pub/grch37/release-101/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz"
IEDB_URL = 'http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip'


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
        self._set_netmhcpan_alleles_mouse()
        self._set_netmhc2pan_alleles_mouse()
        iedb_download_timestamp, iedb_hash = self._set_iedb()
        hash_human, hash_mouse = self._set_proteome()
        hla_url, hla_version, hla_hash = self._set_ipd_imgt_hla_database()
        self._set_h2_resource()
        if self.install_r_dependencies:
            self._install_r_dependencies()
        else:
            logger.warning("R dependencies will need to be installed manually")
        self._save_resources_versions(
            iedb_timestamp=iedb_download_timestamp,
            iedb_hash=iedb_hash,
            hla_url=hla_url,
            hla_version=hla_version,
            hla_hash=hla_hash,
            human_proteome_hash=hash_human,
            mouse_proteome_hash=hash_mouse,
        )

    def _save_resources_versions(
            self, iedb_timestamp, iedb_hash, hla_url, hla_version, hla_hash, human_proteome_hash, mouse_proteome_hash):
        resources_version_file = os.path.join(self.reference_folder, RESOURCES_VERSIONS)
        resources_version = [
            Resource(name="netMHCpan", version="4.1"),
            Resource(name="netMHCIIpan", version="4.0"),
            Resource(name="mixMHCpred", version="2.1"),
            Resource(name="mixMHC2pred", version="1.2"),
            Resource(name="IEDB", version=iedb_timestamp, url=IEDB_URL, hash=iedb_hash),
            Resource(name="Human Ensembl proteome", version=HUMAN_PROTEOME_VERSION, url=HUMAN_PROTEOME,
                     hash=human_proteome_hash),
            Resource(name="Mouse Ensembl proteome", version=MOUSE_PROTEOME_VERSION, url=MOUSE_PROTEOME,
                     hash=mouse_proteome_hash),
            Resource(name="IMGT/HLA database", version=hla_version, url=hla_url, hash=hla_hash)
        ]
        json.dump([r.to_dict() for r in resources_version], open(resources_version_file, "w"), indent=4)

    def _set_netmhcpan_alleles(self):
        # available MHC alleles netMHCpan
        # $NEOFOX_NETMHCPAN -listMHC | grep "HLA-" > "$NEOFOX_REFERENCE_FOLDER"/MHC_available.csv
        logger.info("Fetching available MHC alleles from NetMHCpan for human")
        available_alleles_file = os.path.join(
            self.reference_folder, NETMHCPAN_AVAILABLE_ALLELES_FILE
        )
        cmd = '{netmhcpan} -listMHC | grep "HLA-" | grep -v "#" > {available_alleles_file}'.format(
            netmhcpan=self.config.net_mhc_pan,
            available_alleles_file=available_alleles_file,
        )
        self._run_command(cmd)

    def _set_netmhcpan_alleles_mouse(self):
        logger.info("Fetching available MHC alleles from NetMHCpan for mouse")
        available_alleles_file = os.path.join(
            self.reference_folder, NETMHCPAN_AVAILABLE_ALLELES_MICE_FILE
        )
        cmd = '{netmhcpan} -listMHC | grep -e "^H-2-" | grep -v "#" > {available_alleles_file}'.format(
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

    def _set_netmhc2pan_alleles_mouse(self):
        logger.info("Fetching available alleles from NetMHC2pan")
        available_alleles_file = os.path.join(
            self.reference_folder, NETMHC2PAN_AVAILABLE_ALLELES_MICE_FILE
        )
        cmd = '{netmhc2pan} -list | grep -e "^H-2-" | grep -v "#" > {available_alleles_file}'.format(
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
        cmd = 'wget "{}" -O {}'.format(IEDB_URL, iedb_zip)
        self._run_command(cmd)
        download_timestamp = datetime.today().strftime('%Y%m%d%H%M%S')

        # unzip IEDB
        path_to_iedb_folder = os.path.join(self.reference_folder, IEDB_FOLDER)
        cmd = "unzip -o {iedb_zip} -d {iedb_folder}".format(
            iedb_zip=iedb_zip, iedb_folder=path_to_iedb_folder
        )
        self._run_command(cmd)

        # transforms IEDB into fasta
        tcell_full_iedb_file = os.path.join(self.reference_folder, IEDB_FOLDER, "tcell_full_v3.csv")
        hash = self._get_md5_hash(tcell_full_iedb_file)
        iedb_builder = IedbFastaBuilder(tcell_full_iedb_file)

        iedb_fasta_homo_sapiens = os.path.join(self.reference_folder, IEDB_FOLDER, IEDB_FASTA_HOMO_SAPIENS)
        iedb_builder.build_fasta(
            organism="Homo sapiens",
            process_type="Occurrence of infectious disease",
            output_file=iedb_fasta_homo_sapiens)

        # run makeblastdb on iedb fasta
        cmd = "{makeblastdb} -in {iedb_fasta} -dbtype prot -out {iedb_folder}".format(
            makeblastdb=self.config.make_blastdb,
            iedb_fasta=iedb_fasta_homo_sapiens,
            iedb_folder=os.path.join(path_to_iedb_folder, IEDB_BLAST_PREFIX_HOMO_SAPIENS),
        )
        self._run_command(cmd)

        # TODO: clarify which value of "Process Type" we need
        # values: ['Administration in vivo', 'No immunization', 'Administration in vivo to prevent or reduce disease',
        # 'Administration in vivo to cause disease', 'Unknown', 'Occurrence of disease',
        # 'Occurrence of autoimmune disease', 'Occurrence of cancer', 'Exposure without evidence for disease',
        # 'Transplant/transfusion', 'Environmental exposure to endemic/ubiquitous agent without evidence for disease']
        iedb_fasta_mus_musculus = os.path.join(self.reference_folder, IEDB_FOLDER, IEDB_FASTA_MUS_MUSCULUS)
        iedb_builder.build_fasta(
            organism="Mus musculus",
            process_type="Occurrence of disease",
            output_file=iedb_fasta_mus_musculus)

        # run makeblastdb on iedb fasta
        cmd = "{makeblastdb} -in {iedb_fasta} -dbtype prot -out {iedb_folder}".format(
            makeblastdb=self.config.make_blastdb,
            iedb_fasta=iedb_fasta_mus_musculus,
            iedb_folder=os.path.join(path_to_iedb_folder, IEDB_BLAST_PREFIX_MUS_MUSCULUS),
        )
        self._run_command(cmd)
        return download_timestamp, hash

    def _get_md5_hash(self, filepath):
        file_hash = hashlib.md5()
        with open(filepath, "rb") as f:
            while True:
                chunk = f.read(8192)
                if not chunk:
                    break
                file_hash.update(chunk)
        return file_hash.hexdigest()

    def _set_proteome(self):
        # human proteome database
        # mkdir "$NEOFOX_REFERENCE_FOLDER"/proteome_db
        # wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O "$NEOFOX_REFERENCE_FOLDER"/proteome_db/Homo_sapiens.fa.gz
        # gunzip "$NEOFOX_REFERENCE_FOLDER"/proteome_db/Homo_sapiens.fa.gz
        # $NEOFOX_MAKEBLASTDB -in "$NEOFOX_REFERENCE_FOLDER"/proteome_db/Homo_sapiens.fa -dbtype prot -parse_seqids -out "$NEOFOX_REFERENCE_FOLDER"/proteome_db/homo_sapiens

        logger.info("Configuring the proteome DB")

        os.makedirs(os.path.join(self.reference_folder, PROTEOME_DB_FOLDER))

        # installs Homo sapiens proteome
        # url_human = "ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
        url_human = HUMAN_PROTEOME
        hash_human = self._prepare_proteome(
            url=url_human,
            proteome_file_name=HOMO_SAPIENS_FASTA,
            proteome_prefix=PREFIX_HOMO_SAPIENS,
            proteome_pickle_file_name=HOMO_SAPIENS_PICKLE)

        # installs Mus musculus proteome
        url_mouse = MOUSE_PROTEOME
        hash_mouse = self._prepare_proteome(
            url=url_mouse,
            proteome_file_name=MUS_MUSCULUS_FASTA,
            proteome_prefix=PREFIX_MUS_MUSCULUS,
            proteome_pickle_file_name=MUS_MUSCULUS_PICKLE)

        return hash_human, hash_mouse

    def _prepare_proteome(self, url, proteome_file_name, proteome_prefix, proteome_pickle_file_name):
        # download proteome
        proteome_compressed_file = os.path.join(
            self.reference_folder, PROTEOME_DB_FOLDER, "%s.gz" % proteome_file_name
        )

        cmd = "wget {url} -O {proteome_file}".format(
            url=url, proteome_file=proteome_compressed_file
        )
        self._run_command(cmd)
        hash = self._get_md5_hash(proteome_compressed_file)

        cmd = "gunzip -f {proteome_file}".format(proteome_file=proteome_compressed_file)
        self._run_command(cmd)
        proteome_file = os.path.join(
            self.reference_folder, PROTEOME_DB_FOLDER, proteome_file_name
        )
        output_folder = os.path.join(
            self.reference_folder, PROTEOME_DB_FOLDER, proteome_prefix
        )

        # prepare proteome for blast database and for pickle
        prepared_proteome = []
        proteome_database = []
        for record in SeqIO.parse(proteome_file, "fasta"):
            seq = str(record.seq).replace("*", "")
            record.seq = Seq(seq)
            proteome_database.append(record)
            prepared_proteome.append(seq)
        SeqIO.write(proteome_database, proteome_file, "fasta")

        cmd = "{makeblastdb} -in {proteome_file} -dbtype prot -parse_seqids -out {output_folder}".format(
            makeblastdb=self.config.make_blastdb,
            proteome_file=proteome_file,
            output_folder=output_folder,
        )
        self._run_command(cmd)

        proteome_pickle = os.path.join(
            self.reference_folder, PROTEOME_DB_FOLDER, proteome_pickle_file_name
        )
        outfile = open(proteome_pickle, 'wb')
        pickle.dump("\n".join(prepared_proteome), outfile)
        outfile.close()

        return hash

    def _set_ipd_imgt_hla_database(self):
        logger.info("Adding the available human HLA allele database")
        allele_list = os.path.join(self.reference_folder, HLA_DATABASE_AVAILABLE_ALLELES_FILE)
        url = os.environ.get(
            NEOFOX_HLA_DATABASE_ENV, IMGT_HLA_DB_URL)
        cmd = "wget {url} -O {allele_list}".format(url=url, allele_list=allele_list)
        self._run_command(cmd)
        hash = self._get_md5_hash(allele_list)

        version = None
        with open(allele_list, 'r') as fd:
            while True:
                line = fd.readline()
                if "version" in line:
                    version = line.split(" ")[-1]
                    break

        return url, version, hash

    def _set_h2_resource(self):
        logger.info("Adding the available mouse MHC alleles resource")
        source_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), H2_DATABASE_AVAILABLE_ALLELES_FILE)
        target_file = os.path.join(self.reference_folder, H2_DATABASE_AVAILABLE_ALLELES_FILE)
        copyfile(source_file, target_file)

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
    def __init__(self, input_file):
        self.input_file = input_file

    def build_fasta(self, organism, process_type, output_file):
        # read IEDB input file
        iedb = pd.read_csv(self.input_file, skiprows=1)

        # filter entries
        filtered_iedb = iedb[
            (iedb["Name"].str.contains(organism))
            & (iedb["Object Type"] == "Linear peptide")
            & (iedb["Process Type"] == process_type)
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
        with open(output_file, "w") as fasta:
            for index, row in filtered_iedb.iterrows():
                fasta.write(">{header}\n".format(header=str(row["fasta_header"])))
                fasta.write("{sequence}\n".format(sequence=str(row["seq"])))
