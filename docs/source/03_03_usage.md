# Usage

There are two ways to use NeoFox for annotation of neoantigen candidates with neoantigen features: directly from the [command line](#command-line) or [programmatically](#api). 

## Command line

To call NeoFox from the command line, use the following command. Make sure that the requirements have been added to PATH as described [here](02_installation.md) or add a config file as described below:  

````commandline
neofox --candidate-file/--json-file neoantigens_candidates.tab/neoantigens_candidates.json --patient-data/--patient-data-json patient_data.txt/patient_data.json --output-folder /path/to/out --output-prefix out_prefix [--patient-id] [--with-table] [--with-json] [--num_cpus] [--affinity-threshold]
````

where:
- `--candidate-file`: tab-separated values table with neoantigen candidates represented by long mutated peptide sequences as described [here](03_01_input_data.md#tabular-file-format)
- `--json-file`: JSON file neoantigens in NeoFox model format as  described [here](03_01_input_data.md#json-file-format)
- `--patient-id`: patient identifier (*optional*, this will be used if the patient id the column `patient` is missing the candidate input file)
- `--patient-data`: a table of tab separated values containing metadata on the patient as  described [here](03_01_input_data.md#file-with-patient-information)
- `--output-folder`: path to the folder to which the output files should be written 
- `--output-prefix`: prefix for the output files (*optional*)
- `--with-table`: output file in [short-wide](03_02_output_data.md#tabular-format) tabular format (*default*, *optional*)
- `--with-json`: output file in [JSON](03_02_output_data.md#json-format) format (*optional*)
- `--num_cpus`: number of CPUs to use (*optional*)
- `--config`: a config file with the paths to dependencies as shown below  (*optional*)
- `--affinity-threshold`: a affinity value (*optional*) neoantigen candidates with a best predicted affinity greater than or equal than this threshold will be not annotated with features that specifically model
                        neoepitope recognition. A threshold that is commonly used is 500 nM. 
- `--organism`: the organism to which the data corresponds. Possible values: [human, mouse]. Default value: human


**PLEASE NOTE THE FOLLOWING HINTS**:   
- provide the neoantigen candidate file either as `--candidate-file` or `--json-file` 
- if no specific output format is selected, the output will be written in [short-wide](03_02_output_data.md#short-wide-format) format
- if all expression values related to a patient are NA, imputated expression will be used for the relevant features

**EXAMPLE**  
This is an example to call NeoFox with a candidate-file and obtaining the annotated neoantigen candidates in [short-wide](03_02_output_data.md#short-wide-format) format:  

````commandline
neofox --candidate-file neoantigens_candidates.tab --patient-id Ptx --patient-data patient_data.tab --output-folder /path/to/out --output-prefix test
````

The optional **config** file with the paths to the dependencies can look like this:  
````commandline
NEOFOX_REFERENCE_FOLDER=path/to/reference/folder
NEOFOX_RSCRIPT=`which Rscript`
NEOFOX_BLASTP=path/to/ncbi-blast-2.10.1+/bin/blastp
NEOFOX_NETMHCPAN=path/to/netMHCpan-4.1/netMHCpan
NEOFOX_NETMHC2PAN=path/to/netMHCIIpan-4.0/netMHCIIpan
NEOFOX_MIXMHCPRED=path/to/MixMHCpred-2.1/MixMHCpred
NEOFOX_MIXMHC2PRED=path/to/MixMHC2pred-1.2/MixMHC2pred_unix
NEOFOX_MAKEBLASTDB=path/to/ncbi-blast-2.8.1+/bin/makeblastdb
NEOFOX_PRIME=/path/to/PRIME/PRIME
````

### Running from docker

In order to run the command line in a docker image, all of the above applies but
some additional steps are required.

If the docker image is named `neofox-docker`, run as follows: `docker run neofox-docker neofox --help`

In order to copy the NeoFox input and output data to and from the docker container, a docker volume needs to be created for 
mapping a folder in the host to a folder in the container using the `-v VOLUME_NAME:ABSOLUTE_FOLDER_IN_CONTAINER` argument.

First create a volume:
```
docker volume create neofox-volume
```

Identify the folder where the volume is mounted in the host:
```
$ docker volume inspect neofox-volume
[
    {
        "CreatedAt": "2021-03-25T21:40:23+01:00",
        "Driver": "local",
        "Labels": {},
        "Mountpoint": "/var/snap/docker/common/var-lib-docker/volumes/neofox-volume/_data",
        "Name": "neofox-volume",
        "Options": {},
        "Scope": "local"
    }
]
```

In the case above the folder is `/var/snap/docker/common/var-lib-docker/volumes/neofox-volume/_data`.
Copy the input data into that folder.

Now, NeoFox can be run as follows mounting the volume as indicated. 
Note that the output folder needs to be specified within the volume, if the output from NeoFox should be recovered.
```
docker run -v neofox-volume:/app/data neofox-docker \
neofox --candidate-file /app/data/test_model_file.txt \
--patient-data /app/data/test_patient_info.txt \
--output-folder /app/data/output
```

## API

NeoFox can be used programmatically and by that integrated into existing tools. 
Here, we will explain the use of NeoFox by API in short with the help of a dummy example that includes building models 
from scratch (For a detailed description, please refer to [this notebook](notebooks/api_usage.ipynb)). 
The models can be created based on files too. 
In this case, ignore step 2-5 and refer to the note on the bottom of this paragraph.

NeoFox can be executed as follows.
    
```python
annotations = NeoFox(neoantigens=[neoantigen], patients=[patient], num_cpus=2).get_annotations()
```

Follow the next detailed steps.

**NOTE**: the environment variables need to be loaded. Hint: use `dotenv.load_dotenv()`.

### Create a neoantigen object  

Create a neoantigen candidate using `NeoantigenFactory`.
The data will be internally validated.
Additional annotations with custom names are supported.

```python
from neofox.model.factories import NeoantigenFactory

# create a neoantigen candidate using the factory
neoantigen = NeoantigenFactory.build_neoantigen(
    mutated_xmer="AAAAAAAAAAAAARAAAAAAAAAAAAA",
    wild_type_xmer="AAAAAAAAAAAAAMAAAAAAAAAAAAA",
    patient_identifier="Ptx", 
    rna_expression=0.52, 
    rna_variant_allele_frequency=0.88, 
    dna_variant_allele_frequency=0.29,
    my_custom_annotation="add any custom annotation as additional fields with any name"
)
```   

### Create a patient model

In order to parse MHC alleles and being able to normalize them into a standard nomenclature, load the following resources.
```python
from neofox.references.references import ReferenceFolder

reference_folder = ReferenceFolder(organism='human')
```

To parse mouse H-2 alleles use `ReferenceFolder(organism='mouse')`.

Create a patient model based on models for MHC I and MHC II alleles.
```python
from neofox.model.factories import PatientFactory

patient = PatientFactory.build_patient(
    identifier="Ptx",
    mhc_alleles=["HLA-A*01:01:02:03N", "HLA-A*01:02:02:03N", "HLA-B*01:01:02:03N", "HLA-B*01:01:02:04N", "HLA-C*01:01"],
    mhc2_alleles=["HLA-DPA1*01:01", "HLA-DPA1*01:02", "HLA-DPB1*01:01", "HLA-DPB1*01:01", "HLA-DRB1*01:01", "HLA-DRB1*01:01"],
    mhc_database=reference_folder.get_mhc_database()
)
```

 **WARNING**: alleles in homozygous state need to be provided twice, otherwise they are considered as hemizygous. 
 For instance `["HLA-A*01:01"]` would be interpreted as hemizygous and `["HLA-A*01:01", "HLA-A*01:01"]` as homozygous.

   
### Run NeoFox  

Run NeoFox by passing the neoantigen and patients object to get the neoantigen features. 
The output is a list of type `NeoantigenAnnotations`:  

```python
from neofox.neofox import NeoFox

annotated_neoantigens = NeoFox(neoantigens=[neoantigen], patients=[patient], num_cpus=2).get_annotations()
```  

where:
       - `neoantigens`: a list of neoantigen objects  
       - `patients`: a list of patient objects  
       - `num_cpus`: number of CPUs to use (*optional*)
       
### Transformation of the output   
    
Depending on the use case, the user can transform the resulting neoantigen feature annotations into the formats described [here](03_02_output_data.md).

```python
from neofox.model.conversion import ModelConverter

# Pandas data frame
annotations_table = ModelConverter.annotations2table(neoantigens=[neoantigen])

# JSON 
neoantigen_json = ModelConverter.objects2json(model_objects=[neoantigen])
```
   
   
**PLEASE NOTE THE FOLLOWING HINTS**:   

- process multiple neoantigens by passing a list of neoantigens and a list of patients to `NeoFox().get_annotations()`.
- only the transformation of the annotation with `ModelConverter.annotations2short_wide_table()` will keep both information of neoantigen candidates and annotated features values. with Neoantigen objects can be transformed into other formats, too when `ModelConverter.annotations2tall_skinny_table()` or `ModelConverter.objects2json()` were used. In case of our example:  

```python
# convert neoantigens into data frame
neoantigens_df = ModelConverter.objects2dataframe(model_objects=[neoantigen])

# convert neoantigens into JSON format 
neoantigens_json = ModelConverter.objects2json(model_objects=[neoantigen])
```   
- instead of creating neoantigen or patient models, tabular or json files containing this information can be passed:  
  The neoantigen candidates can be provided in [candidate-file format](03_01_input_data.md#tabular-file-format)

```python
model_file = "/path/to/neoantigen_candidates.tab"
neoantigens = ModelConverter.parse_neoantigens_file(neoantigens_file=model_file)
```
  
 
or in [JSON format](03_01_input_data.md#json-file-format). 

```python
json_file = "/path/to/neoantigen_candidates.json"
neoantigens = ModelConverter.parse_neoantigens_json_file(json_file=json_file)  
```  

The patient information should be provided in [tabular format](03_01_input_data.md#file-with-patient-information)

```python
patient_file = "/path/to/patients.tab"
patients = ModelConverter.parse_patients_file(patients_data)
```  
  
Then, run NeoFox as explained before by calling:

```python
neoantigens_annotated = NeoFox(neoantigens=neoantigens, patients=patients, num_cpus=2).get_annotations()
```

## Performance

As indicated above NeoFox can run in parallel using the parameter `--num-cpus`. 
Each CPU will process one neoantigen candidate at a time, thus NeoFox uses only as many CPUs as candidats are to be processed.

We processed several simulated datasets with 10, 100, 1000 and 10000 neoantigen candidates on 1, 5, 10 and 50 CPUs. We obtained 
that the average time to process a single candidate in a single CPU takes 37.516 seconds, with a standard deviation of 
6.739 seconds. No significant overhead due to parallelization was observed. 
In terms of memory the application uses less than 0.5 GB for up to 1000 neoantigen candidates irrespective of the number of CPUs used. 
The memory use grows to around 2.5 GB when processing 10000 candidates. 

![Neofox model](../figures/performance_1.jpg)

If either MHC I or II alleles are not provided at all for a given patient the computation will be lighter as no 
annotations run for the missing MHC. Likewise, if the optional tools are unset performance improves.


## A note on the support for mouse

NeoFox was originally conceived for Homo sapiens data and later extended to Mus musculus.
There are some fundamental differences to take into account when processing mouse data.
First, not all annotations are available for mouse. 
MixMHCpred, MixMHC2pred, PRIME and HEX are not available for mouse. 
At the time of this writing the expression imputation is not available either.
The MHC for Mus musculus, H-2, is not described with the same level of detail as HLA is. 
Furthermore, the Mus musculus strains used in laboratory experimentation are different from the wild type Mus musculus 
and their variability is much more limited. 
NetMHCpan and netMHC2pan, and by extension NeoFox, support a subset of the H-2 alleles found in laboratory mice which
is again a small subset of the wild type.
A consequence of this is that MHC II is highly simplified, the genes for chains alpha and beta are considered to be 
part of the same haplotype always. Furthermore, only homozygosity is considered. 
Thus there is only one possible MHC II isoform for each pair of genes, as opposed to four in human.

H-2 alleles are not standardized as HLA alleles are through the WHO Nomenclature Committee for Factors of the HLA System.
In NeoFox we represent H-2 alleles as follows.
For MHC I genes K, D and L: H2K, H2D and H2L
For MHC II genes A and E: H2A and H2E
Then a last small case single letter (eg: d, k, p) with an optional number (eg: d1, p2) represent a given allele.

This is an example of H-2 alleles: H2Kd, H2Dd, H2Lp 


