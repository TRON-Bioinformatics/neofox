# Installation

This guide contains two alternatives to install NeoFox:
- Building a docker image that automates the installation into a container (**NOTE**: the docker recipe is not supported in neofox-v1.1.0. Please use an older version (<v1.1.0) for building the docker image at the moment.)
- A set of detailed step by step installation instructions without docker

The first approach has the lowest entry barrier to use NeoFox as a command line tool.
While the second provides access to the command line tool and allows the integration of the NeoFox API.

## Build and run the docker image

Clone the repository: `git clone git@github.com:TRON-Bioinformatics/neofox.git`

Move into the neofox folder: `cd neofox`

To download NetMHCpan and NetMHCIIpan, each user must explicitly accept the software license. Thus, we cannot
distribute the software or provide a direct URL to download it. Please make sure to download the right version from 
the sites indicated below.

- NetMHCpan-4.1: https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1 (`netMHCpan-4.1b.Linux.tar.gz`)
- NetMHCIIpan-4.0: https://services.healthtech.dtu.dk/software.php (`netMHCIIpan-4.0.Linux.tar.gz`)

Please also check the licences of the other third party dependencies ( e.g. listed in the docker recipe `Dockerfile`). 

Store these in the root folder of the repository, next to the `Dockerfile`. Do not rename the installer files.

Build the docker image: `docker build --platform linux/amd64 --tag neofox-docker .`

Run NeoFox: `docker run neofox-docker neofox --help`

See the usage guide [here](03_03_usage.md) for further details.


## Step by step guide without docker

These installation instructions were tested on Ubuntu 18.04.

Python 3.7 or 3.8 should be preinstalled.

The libz compression development library is required. This can be installed in Ubuntu as follows:
```
apt-get install libbz2-dev
```

### Install NeoFox

Install from PyPI:
```
pip install neofox
```

or install from bioconda:
```
conda install bioconda::neofox
```

### Install third-party dependencies

**NOTE**: Please, check the licences of third-party dependencies. 

#### Install R

R 3.6.0 is required.

Optionally set the environment variable pointing to `Rscript`, otherwise neofox will look for it in the path.
```
export NEOFOX_RSCRIPT=`which Rscript`
```

**NOTE**: when installing from conda this dependency is already installed.

#### Install BLASTP

The version of BLASTP that was tested is 2.10.1, other versions may work but that is untested.
```
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz
tar -xvf ncbi-blast-2.10.1+-x64-linux.tar.gz
```

Optionally set the environment variable pointing to `blastp`, otherwise neofox will look for it in the path.
```
export NEOFOX_BLASTP=/path/to/ncbi-blast-2.10.1+/bin/blastp
```

**NOTE**: when installing from conda this dependency is already installed.

#### Install NetMHCpan-4.1

NetMHCpan-4.1 can be downloaded by academic users from https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1

```
tar -xvf netMHCpan-4.1b.Linux.tar.gz
cd netMHCpan-4.1
wget https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz
tar -xvf data.tar.gz
```

Optionally set the environment variable pointing to `netMHCpan`, otherwise neofox will look for it in the path.
```
export NEOFOX_NETMHCPAN=/path/to/netMHCpan-4.1/netMHCpan
```

Configure NetMHCpan as explained in the file `netMHCpan-4.1/netMHCpan-4.1.readme`


#### Install NetMHCIIpan-4.0

NetMHCIIpan-4.0 can be downloaded by academic users from https://services.healthtech.dtu.dk/software.php

```
tar -xvf netMHCIIpan-4.0.Linux.tar.gz
cd netMHCIIpan-4.0
# download the data
wget http://www.cbs.dtu.dk/services/NetMHCIIpan-4.0/data.Linux.tar.gz
tar -xvf data.Linux.tar.gz
# install tcsh shell interpreter if not available yet
sudo apt-get install tcsh
```

Optionally set the environment variable pointing to `netMHCIIpan`, otherwise neofox will look for it in the path.
```
export NEOFOX_NETMHC2PAN=/path/to/netMHCIIpan-4.0/netMHCIIpan
```

Configure NetMHCIIpan-4.0 as explained in the file `netMHCIIpan-4.0/netMHCIIpan-4.0.readme`
         

#### Install MixMHCpred-2.2 (recommended but optional)

```
wget https://github.com/GfellerLab/MixMHCpred/archive/refs/tags/v2.2.tar.gz
tar -xvf v2.2.tar.gz
```

Set the environment variable pointing to `MixMHCpred`, there will be no search in the path as the installation folder
is also needed to determine the supported alleles.
```
export NEOFOX_MIXMHCPRED=/path/to/MixMHCpred-2.2/MixMHCpred
```

Configure MixMHCpred-2.2 as explained in the file `MixMHCpred-2.2/README`



#### Install MixMHC2pred-2.0.2 (recommended but optional)

```
wget https://github.com/GfellerLab/MixMHC2pred/archive/refs/tags/v2.0.2.2.tar.gz
tar -xvf v2.0.2.2.tar.gz
```

Set the environment variable pointing to `MixMHC2pred_unix`, there will be no search in the path as the installation 
folder is also needed to determine the supported alleles.
```
export NEOFOX_MIXMHC2PRED=`pwd`/MixMHC2pred-2.0.2/MixMHC2pred_unix
```

#### PRIME-1.0 (recommended but optional)

```
wget https://github.com/GfellerLab/PRIME/archive/refs/tags/v2.0.tar.gz
tar -xvf v2.0.tar.gz
```

Set the environment variable pointing to `PRIME`, there will be no search in the path as the installation folder
is also needed to determine the supported alleles.
```
export NEOFOX_PRIME==`pwd`/PRIME-master/PRIME
```

Configure PRIME as explained in the file `PRIME-master/README`

### Configuration of the reference folder 

To configure the reference folder, set the environment variables for `makeblastdb`, NetMHCpan, NetMHCIIpan and Rscript,
 or alternatively rely on these being fetched from the path:

```
export NEOFOX_MAKEBLASTDB=`pwd`/ncbi-blast-2.10.1+/bin/makeblastdb
export NEOFOX_RSCRIPT=`which Rscript`
export NEOFOX_NETMHCPAN=`pwd`/netMHCpan-4.1/netMHCpan
export NEOFOX_NETMHC2PAN=`pwd`/netMHCIIpan-4.0/netMHCIIpan
```

Furthermore, a list of available MHC alleles is required. Optionally, you can provide the URL to the IPD-IMGT/HLA database CSV table, see releases here https://www.ebi.ac.uk/ipd/imgt/hla/docs/release.html. 
If not provided the default value is the latest version at the time of this writing https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.3430.txt

```
export NEOFOX_HLA_DATABASE=https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.3430.txt
```

Run the following to configure the NeoFox reference folder:
```
neofox-configure --reference-folder /your/neofox/folder [--install-r-dependencies --install_mouse_mixmhc2pred]
```

**NOTE**: when installing from conda `--install-r-dependencies` is not needed. 

The above command will download and transform several resources and store in the annotations metadata their version, MD5 checksum and 
download timestamp. 


To run NeoFox on data from mouse with MixMHC2pred, mouse-specific PMWs are required. For such use cases the reference folder needs to be configured with `--install_mouse_mixmhc2pred` (see also )

Depending on your use case please check the licences of these third-party resources (see urls in neofox/references/installer.py). 

Unless indicated to the installer by flag `--install-r-dependencies` you will need to install manually some R packages. These packages are the following:
```
lattice
ggplot2
caret
Peptides
doParallel
gbm
Biostrings
```

Add the reference folder to the Path
```
export NEOFOX_REFERENCE_FOLDER=path/to/reference/folder
```

## Test installation   

The user can test if all the installations have been successful by testing NeoFox with some test data. 
The test data can be downloaded here:

* [test_data.tsv](_static/test_data.tsv)
* [test_patients.tsv](_static/test_patients.tsv)

````commandline
neofox --input-file /path/to/test_data.txt --patient-data /path/to/test_patients.txt --output-folder  /path/to/outputfolder --output-prefix test
````

The resulting output files can be compared to the following test output files:

* [test_neoantigen_candidates_annotated.tsv](_static/test_neoantigen_candidates_annotated.tsv)
* [test_neoantigen_candidates_annotated.json](_static/test_neoantigen_candidates_annotated.json)

