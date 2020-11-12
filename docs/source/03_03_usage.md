# Usage

## General information
There are two ways to use NeoFox for annotation of neoantigen candidates with neoantigen features: directly from the command line [command line](#command-line) or [programmatically](#api). 

## Command line

````commandline
neofox --model-file/--candidate-file neoantigens.txt --patient-id Ptx --patient-data patient_data.txt --output-folder /path/to/out --output-prefix out_prefix [--with-short-wide-table] [--with-tall-skinny-table] [--with-json] [--num_cpus]
````

where:
- `--candidate-file`: tab-separated values table with neoantigen candidates represented by long mutated peptide sequences as described [here](03_01_input_data.md#tabular-format)
- `--model-file`: tab-separated values table with neoantigens in Neofox model format described in [protobuf model](neofox/model/neoantigen.proto)
- `--patient-id`: patient identifier (**optional**, this will be used as the patient id the column `patient` is missing the candidate input file)
- `--patient data`: a table of tab separated values containing metadata on the patient
## API
