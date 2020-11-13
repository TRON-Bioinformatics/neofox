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
NeoFox can be used programmatically and by that integrated into existing tools. Here, we will explain step by step the use of NeoFox by API using an example.  

1. **Import requirements**   
    ````python
    from neofox.model.conversion import ModelConverter
    from neofox.model.conversion import ModelValidator
    from neofox.model.neoantigen import Neoantigen, Transcript, Mutation, Patient
    from neofox.neofox import NeoFox
   ````    
2. **Create a neoantigen model**  
    Create a neoantigen candidate model based on Transcript and Mutation model. Initialise each of these models by passing the required information. The following shows a dummy example:
    ````python
    # model the transcript related to the neoantigen candidate
    transcript = Transcript(assembly="hg19", gene="VCAN", identifier="uc003kii.3")
    # model the mutation related to the neoantigen candidate
    mutation = Mutation(position=1007, wild_type_aminoacid="I", mutated_aminoacid="T", left_flanking_region="DEVLGEPSQDILV", right_flanking_region="DQTRLEATISPET")
    # create a neoantigen candidate model using the transcript and mutation model
    neoantigen = Neoantigen(transcript=transcript, mutation=mutation, patient_identifier="Ptx", rna_expression=0.519506894, rna_variant_allele_frequency=0.857142857, dna_variant_allele_frequency=0.294573643)
    ````   
   where:  
       - transcript: Transcript model (explanation of the parameters is provided [here](05_models.md#transcript))  
       - mutation: Mutation model (explanation of the parameters is provided [here](05_models.md#mutation))  
       - neoantigen: Neoantigen candidate model (explanation of the parameters is provided [here](05_models.md#neoantigen))

3. **Validate the neoantigen model**  
    Check for validity of the entered parameters into the neoantigen models and the validity of the full neoantigen model:   
    ````python
    validated_neoantigen = ModelValidator.validate_neoantigen(neoantigen=neoantigen)
   ````

4. **Create a patient model**  
    Create a patient model based on models for MHC I and MHC II alleles. Initialise each of these models by passing the required information. The following shows a dummy example:
    ````python
    # create MHC I model for the patient
    mhc1 = ModelConverter.parse_mhc1_alleles(["HLA-A*01:01:02:03N", "HLA-A*01:02:02:03N", "HLA-B*01:01:02:03N", "HLA-B*01:01:02:04N", "HLA-C*01:01"])
    # create MHC II model for the patient
    mhc2 = ModelConverter.parse_mhc2_alleles(["HLA-DPA1*01:01", "HLA-DPA1*01:02", "HLA-DPB1*01:01", "HLA-DPB1*01:01", "HLA-DRB1*01:01", "HLA-DRB1*01:01"])
    patient = Patient(identifier="P123", is_rna_available=True, mhc1=mhc1, mhc2=mhc2)
   ````
      where:  
       - mhc1: Model of MHC class I alleles (explanation is provided [here](05_models.md#mhc1)). Single alleles should be provided with *at least 4digits* but more digits are allowed.  
       - mhc2: Mutation model (explanation is provided [here](05_models.md#mhc2))  Single alleles should be provided with *at least 4digits* but more digits are allowed  
       - patient: Patient model (explanation of the parameters is provided [here](05_models.md#patient))
       
5. **Validate the patient model**  
    Check for validity of the entered parameters into the patient models and the validity of the full patient model: 
    ````python
    validated_patient = ModelValidator.validate_patient(patient)
   ````
   
6. **Run NeoFox**  
    Run NeoFox by passing the validated neoantigen object and the validated patient object to get the neoantigen features. The output is a list of type `NeoantigenAnnotations`:  
    ````python
    annotations = NeoFox(neoantigens=[validated_neoantigen], patients=[validated_patient], num_cpus=2).get_annotations()
   ````  
      where:  
       - `anotations`: list of type `NeoantigenAnnotations`, i.e. a list of neoantigen features and there values for a given neoantigen candidate (further explanation is provided [here](05_models.md#neoantigenannotations))  
       - `neoantigens`: a list of validated neoantigen objects  
       - `patients`: a list of validated patient objects  
       - `num_cpus`: number of CPUs to use (*optional*)
       
7. **Transformation of output**   
    Depending on the use case, the user can transform the resulting neoantigen feature annotations into the formats described [here](03_02_output_data.md).
    ````python
   # short-wide 
   annotations_sw = ModelConverter.annotations2short_wide_table(neoantigen_annotations=annotations, neoantigens = [validated_neoantigen])
   # tall-skinny
   annotations_ts = ModelConverter.annotations2tall_skinny_table(neoantigen_annotations=annotations)
   # JSON 
   annotations_json = ModelConverter.objects2json(model_objects=annotations)
   ````
   
   Neoantigen obejcts can be transformed into other formats, too. In case, of our example:  
   ````python
    # convert neoantigens into data frame
   neoantigens_df = ModelConverter.objects2dataframe(model_objects=[validated_neoantigen])
   # convert neoantigens into JSON format 
   neoantiges_json = ModelConverter.objects2json(model_objects=[validated_neoantigen]
   ```` 