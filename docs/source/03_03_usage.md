# Usage

## General information
There are two ways to use NeoFox for annotation of neoantigen candidates with neoantigen features: directly from the command line [command line](#command-line) or [programmatically](#api). 

## Command line
To call NeoFox from the command line, please use the following command:  

````commandline
neofox --model-file/--candidate-file/--json-file neoantigens_candidates.tab/neoantigens_candidates.json --patient-id Ptx --patient-data/--patient-data-json patient_data.txt/patient_data.json --output-folder /path/to/out --output-prefix out_prefix [--with-short-wide-table] [--with-tall-skinny-table] [--with-json] [--num_cpus]
````

where:
- `--candidate-file`: tab-separated values table with neoantigen candidates represented by long mutated peptide sequences as described [here](03_01_input_data.md#tabular-format)
- `--model-file`: tab-separated values table with neoantigens in NeoFox model format as  described [here](03_01_input_data.md#tabular-format)
- `--json-file`: JSON file neoantigens in NeoFox model format as  described [here](03_01_input_data.md#json-format)
- `--patient-id`: patient identifier (*optional*, this will be used as the patient id the column `patient` is missing the candidate input file)
- `--patient-data`: a table of tab separated values containing metadata on the patient as  described [here](03_01_input_data.md#file-with-patient-information)
- `--patient-data-json`: a table patient models as described [here](03_01_input_data.md#patient-file-in-json-format)
- `--output-folder`: path to the folder to which the output files should be written 
- `--output-prefix`: prefix for the output files (*optional*)
- `--with-short-wide-table`: output file in [short-wide](03_02_output_data.md#short-wide-format) format (*optional*)
- `--with-tall-skinny-table`: output file in [tall-skinny](03_02_output_data.md#tall-skinny-format) format (*optional*)
- `--with-json`: output file in [JSON](03_02_output_data.md#json-format) format (*optional*)
- `--num_cpus`: number of CPUs to use (*optional*)

**PLEASE NOTE THE FOLLOWING HINTS**:   
- provide the neoantigen candidate file either as `--candidate-file`, `--model-file` or `--json-file` 
- provide the patient data in tabular format format (`--patient-data`) if neoantigen candidates are provided with `--candidate-file` or `--model-file`
- provide the patient data in JSON format format (`--patient-data-json`) if neoantigen candidates are provided with `--json-file` 
- if no specific output format is selected, the output will be written in [short-wide](03_02_output_data.md#short-wide-format) format
- indicate in the `isRnaAvailable` column of the [patient file](03_01_input_data.md#file-with-patient-information) if expression should be imputed for neoantigen candidates of the respective patient  

**EXAMPLE**  
This is an example to call NeoFox with a model-file and obtaining the annotated neoantigen candidates in [short-wide](03_02_output_data.md#short-wide-format) format:  

````commandline
neofox --model-file neoantigens_candidates.tab --patient-id Ptx --patient-data patient_data.tab --output-folder /path/to/out --output-prefix test
````

## API
