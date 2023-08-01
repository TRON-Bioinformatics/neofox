# Data models
<a name="top"></a>

Protocol Buffers is employed to model Neofox's input and output data: neoantigens, Major Histocompatibility Complex (MHC) alleles, patients and output annotations.

![Neofox model](../figures/neofox_model.png)

## Table of Contents

- [neoantigen.proto](#neoantigen.proto)
    - [Annotation](#neoantigen.Annotation)
    - [Annotations](#neoantigen.Annotations)
    - [Mhc1](#neoantigen.Mhc1)
    - [Mhc2](#neoantigen.Mhc2)
    - [Mhc2Gene](#neoantigen.Mhc2Gene)
    - [Mhc2Isoform](#neoantigen.Mhc2Isoform)
    - [MhcAllele](#neoantigen.MhcAllele)
    - [Neoantigen](#neoantigen.Neoantigen)
    - [Patient](#neoantigen.Patient)
    - [PredictedEpitope](#neoantigen.PredictedEpitope)
    - [Resource](#neoantigen.Resource)
  
    - [Mhc1Name](#neoantigen.Mhc1Name)
    - [Mhc2GeneName](#neoantigen.Mhc2GeneName)
    - [Mhc2Name](#neoantigen.Mhc2Name)
    - [Zygosity](#neoantigen.Zygosity)
  
- [Scalar Value Types](#scalar-value-types)



<a name="neoantigen.proto"></a>
<p align="right"><a href="#top">Top</a></p>

## neoantigen.proto



<a name="neoantigen.Annotation"></a>

### Annotation
This is a generic class to hold annotations from Neofox


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| name | [string](#string) |  | The name of the annotation |
| value | [string](#string) |  | The value of the annotation |






<a name="neoantigen.Annotations"></a>

### Annotations
A set of annotations for a neoantigen candidate


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| annotations | [Annotation](#neoantigen.Annotation) | repeated | List of annotations |
| annotator | [string](#string) |  | The annotator |
| annotatorVersion | [string](#string) |  | The version of the annotator |
| timestamp | [string](#string) |  | A timestamp determined when the annotation was created |
| resources | [Resource](#neoantigen.Resource) | repeated | List of resources |






<a name="neoantigen.Mhc1"></a>

### Mhc1
Models MHC I alleles related to the same MHC I gene, i.e. 2 alleles/2 isoforms per gene


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| name | [Mhc1Name](#neoantigen.Mhc1Name) |  | MHC I gene name |
| zygosity | [Zygosity](#neoantigen.Zygosity) |  | Zygosity of the gene |
| alleles | [MhcAllele](#neoantigen.MhcAllele) | repeated | The alleles of the gene (0, 1 or 2) |






<a name="neoantigen.Mhc2"></a>

### Mhc2
Models MHC II alleles related to the same MHC II protein, i.e. 4 isoforms related to 2 genes with 2 alleles each


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| name | [Mhc2Name](#neoantigen.Mhc2Name) |  | MHC II molecule name |
| genes | [Mhc2Gene](#neoantigen.Mhc2Gene) | repeated | List of MHC II genes |
| isoforms | [Mhc2Isoform](#neoantigen.Mhc2Isoform) | repeated | Different combinations of MHC II alleles building different isoforms |






<a name="neoantigen.Mhc2Gene"></a>

### Mhc2Gene
MHC II gene


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| name | [Mhc2GeneName](#neoantigen.Mhc2GeneName) |  | MHC II gene name |
| zygosity | [Zygosity](#neoantigen.Zygosity) |  | Zygosity of the gene |
| alleles | [MhcAllele](#neoantigen.MhcAllele) | repeated | The alleles of the gene (0, 1 or 2) |






<a name="neoantigen.Mhc2Isoform"></a>

### Mhc2Isoform
MHC II isoform


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| name | [string](#string) |  | Name to refer to the MHC II isoform |
| alphaChain | [MhcAllele](#neoantigen.MhcAllele) |  | The alpha chain of the isoform |
| betaChain | [MhcAllele](#neoantigen.MhcAllele) |  | The beta chain of the isoform |






<a name="neoantigen.MhcAllele"></a>

### MhcAllele
MHC allele representation. It does not include non synonymous changes to the sequence, changes in the non coding region
or changes in expression. See http://hla.alleles.org/nomenclature/naming.html for details


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| fullName | [string](#string) |  | HLA full name as provided by the user (e.g.: HLA-DRB1*13:01:02:03N). This will be parsed into name, gene and group. Any digit format is allowed for this field (ie: 4, 6 or 8 digits), 2 digits names are not specific enough for our purpose and thus invalid |
| name | [string](#string) |  | A specific HLA protein (e.g. HLA-DRB1*13:01). Alleles whose numbers differ in group and protein must differ in one or more nucleotide substitutions that change the amino acid sequence of the encoded protein. This name is normalized to avoid different representations of the same allele. For instance both HLA-DRB113:01 and HLA-DRB1*13:01:02:03N will be transformed into their normalised version HLA-DRB1*13:01. This name is also truncated to 4 digits. 2 digits names are not specific enough for our purpose and thus invalid |
| gene | [string](#string) |  | The gene from either MHC I or II (e.g. DRB1, A) (this information is redundant with the Mhc1Gene.name and Mhc2Gene.name but it is convenient to have this at this level too, code will check for data coherence) |
| group | [string](#string) |  | A group of alleles defined by a common serotype ie: Serological antigen carried by an allotype (e.g. 13 from HLA-DRB1*13) |
| protein | [string](#string) |  | A specific protein (e.g.: 02 from HLA-DRB1*13:02) |






<a name="neoantigen.Neoantigen"></a>

### Neoantigen
A neoantigen minimal definition


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| patientIdentifier | [string](#string) |  | Patient identifier |
| gene | [string](#string) |  | The HGNC gene symbol or gene identifier |
| position | [int32](#int32) | repeated | The amino acid position within the neoantigen candidate sequence. 1-based, starting in the N-terminus |
| wildTypeXmer | [string](#string) |  | Amino acid sequence of the WT corresponding to the neoantigen candidate sequence (IUPAC 1 letter codes) |
| mutatedXmer | [string](#string) |  | Amino acid sequence of the neoantigen candidate (IUPAC 1 letter codes) |
| rnaExpression | [float](#float) |  | Expression value of the transcript from RNA data. Range [0, +inf]. |
| imputedGeneExpression | [float](#float) |  | Expression value of the transcript from TCGA data. Range [0, +inf]. |
| dnaVariantAlleleFrequency | [float](#float) |  | Variant allele frequency from the DNA. Range [0.0, 1.0] |
| rnaVariantAlleleFrequency | [float](#float) |  | Variant allele frequency from the RNA. Range [0.0, 1.0] |
| neofoxAnnotations | [Annotations](#neoantigen.Annotations) |  | The NeoFox neoantigen annotations |
| externalAnnotations | [Annotation](#neoantigen.Annotation) | repeated | List of external annotations |
| neoepitopesMhcI | [PredictedEpitope](#neoantigen.PredictedEpitope) | repeated | List of predicted neoepitopes for MHC-I with feature annotation (optional) |
| neoepitopesMhcII | [PredictedEpitope](#neoantigen.PredictedEpitope) | repeated | List of predicted neoepitopes for MHC-II with feature annotation (optional) |






<a name="neoantigen.Patient"></a>

### Patient
The metadata required for analysis for a given patient + its patient identifier


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| identifier | [string](#string) |  | Patient identifier |
| tumorType | [string](#string) |  | Tumor entity in TCGA study abbrevation style as described here: https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations |
| mhc1 | [Mhc1](#neoantigen.Mhc1) | repeated | MHC I classic molecules |
| mhc2 | [Mhc2](#neoantigen.Mhc2) | repeated | MHC II classic molecules |






<a name="neoantigen.PredictedEpitope"></a>

### PredictedEpitope



| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| position | [int32](#int32) |  | Not sure that we need this... this is in the old PredictedEpitope model |
| mutatedPeptide | [string](#string) |  | The mutated peptide |
| wildTypePeptide | [string](#string) |  | Closest wild type peptide |
| alleleMhcI | [MhcAllele](#neoantigen.MhcAllele) |  | MHC I allele |
| isoformMhcII | [Mhc2Isoform](#neoantigen.Mhc2Isoform) |  | MHC II isoform |
| core | [string](#string) |  | MHCII core part of the peptide ligand that primarily interacts with the MHC binding groove, predicted by NetMHCpan/NetMHCIIpan |
| affinityMutated | [float](#float) |  | MHC binding affinity for the mutated peptide. This value is estimated with NetMHCpan in case of MHC-I peptides and NetMHCIIpan in cas of MHC-II peptides |
| rankMutated | [float](#float) |  | MHC binding rank for the mutated peptide. This value is estimated with NetMHCpan in case of MHC-I peptides and NetMHCIIpan in cas of MHC-II peptides |
| affinityWildType | [float](#float) |  | MHC binding affinity for the wild type peptide. This value is estimated with NetMHCpan in case of MHC-I peptides and NetMHCIIpan in cas of MHC-II peptides |
| rankWildType | [float](#float) |  | MHC binding rank for the wild type peptide. This value is estimated with NetMHCpan in case of MHC-I peptides and NetMHCIIpan in cas of MHC-II peptides |
| neofoxAnnotations | [Annotations](#neoantigen.Annotations) |  | The NeoFox neoantigen annotations |
| patientIdentifier | [string](#string) |  | Patient identifier |
| gene | [string](#string) |  | The HGNC gene symbol or gene identifier |
| rnaExpression | [float](#float) |  | Expression value of the transcript from RNA data. Range [0, +inf]. |
| imputedGeneExpression | [float](#float) |  | Expression value of the transcript from TCGA data. Range [0, +inf]. |
| dnaVariantAlleleFrequency | [float](#float) |  | Variant allele frequency from the DNA. Range [0.0, 1.0] |
| rnaVariantAlleleFrequency | [float](#float) |  | Variant allele frequency from the RNA. Range [0.0, 1.0] |
| externalAnnotations | [Annotation](#neoantigen.Annotation) | repeated | External annotations for neoepitope mode. |






<a name="neoantigen.Resource"></a>

### Resource
This is a class to track the version of an annotation resource


| Field | Type | Label | Description |
| ----- | ---- | ----- | ----------- |
| name | [string](#string) |  | The name of the resource |
| version | [string](#string) |  | The version of the resource |
| url | [string](#string) |  | The URL of the resource if applicable |
| hash | [string](#string) |  | The MD5 hash of the resource if applicable. This may be used when version is not available |
| download_timestamp | [string](#string) |  | The timestamp when the download happened |





 <!-- end messages -->


<a name="neoantigen.Mhc1Name"></a>

### Mhc1Name
Valid names for MHC I classic genes
Mus musculus gene names are preceded by the prefix H2 to avoid naming collisions.

| Name | Number | Description |
| ---- | ------ | ----------- |
| A | 0 | Homo sapiens |
| B | 1 |  |
| C | 2 |  |
| H2K | 3 | Mus musculus |
| H2D | 4 |  |
| H2L | 5 |  |



<a name="neoantigen.Mhc2GeneName"></a>

### Mhc2GeneName
Valid names for MHC II classic genes.
DRA is not included in this list as it does not have much variability in the population and for our purpose is
considered constant.
For Mus musculus we do not represent alpha and beta chains as they are homozygotes at all their MHC loci.
Hence, they can be treated as a single gene, like DR is for HLA.
See http://www.imgt.org/IMGTrepertoireMH/Polymorphism/haplotypes/mouse/MHC/Mu_haplotypes.html
Mus musculus gene names are preceded by the prefix H2 to avoid naming collisions.

| Name | Number | Description |
| ---- | ------ | ----------- |
| DRB1 | 0 | Homo sapiens |
| DPA1 | 1 |  |
| DPB1 | 2 |  |
| DQA1 | 3 |  |
| DQB1 | 4 |  |
| H2A | 5 | Mus musculus |
| H2E | 6 |  |



<a name="neoantigen.Mhc2Name"></a>

### Mhc2Name
Valid names for MHC II classic molecules

| Name | Number | Description |
| ---- | ------ | ----------- |
| DR | 0 |  |
| DP | 1 |  |
| DQ | 2 |  |
| H2A_molecule | 3 |  |
| H2E_molecule | 4 |  |



<a name="neoantigen.Zygosity"></a>

### Zygosity
The zygosity of a given gene

| Name | Number | Description |
| ---- | ------ | ----------- |
| HOMOZYGOUS | 0 | Two equal copies of the gene |
| HETEROZYGOUS | 1 | Two different copies of the gene |
| HEMIZYGOUS | 2 | Only one copy of the gene |
| LOSS | 3 | No copy of the gene |


 <!-- end enums -->

 <!-- end HasExtensions -->

 <!-- end services -->



## Scalar Value Types

| .proto Type | Notes | C++ | Java | Python | Go | C# | PHP | Ruby |
| ----------- | ----- | --- | ---- | ------ | -- | -- | --- | ---- |
| <a name="double" /> double |  | double | double | float | float64 | double | float | Float |
| <a name="float" /> float |  | float | float | float | float32 | float | float | Float |
| <a name="int32" /> int32 | Uses variable-length encoding. Inefficient for encoding negative numbers – if your field is likely to have negative values, use sint32 instead. | int32 | int | int | int32 | int | integer | Bignum or Fixnum (as required) |
| <a name="int64" /> int64 | Uses variable-length encoding. Inefficient for encoding negative numbers – if your field is likely to have negative values, use sint64 instead. | int64 | long | int/long | int64 | long | integer/string | Bignum |
| <a name="uint32" /> uint32 | Uses variable-length encoding. | uint32 | int | int/long | uint32 | uint | integer | Bignum or Fixnum (as required) |
| <a name="uint64" /> uint64 | Uses variable-length encoding. | uint64 | long | int/long | uint64 | ulong | integer/string | Bignum or Fixnum (as required) |
| <a name="sint32" /> sint32 | Uses variable-length encoding. Signed int value. These more efficiently encode negative numbers than regular int32s. | int32 | int | int | int32 | int | integer | Bignum or Fixnum (as required) |
| <a name="sint64" /> sint64 | Uses variable-length encoding. Signed int value. These more efficiently encode negative numbers than regular int64s. | int64 | long | int/long | int64 | long | integer/string | Bignum |
| <a name="fixed32" /> fixed32 | Always four bytes. More efficient than uint32 if values are often greater than 2^28. | uint32 | int | int | uint32 | uint | integer | Bignum or Fixnum (as required) |
| <a name="fixed64" /> fixed64 | Always eight bytes. More efficient than uint64 if values are often greater than 2^56. | uint64 | long | int/long | uint64 | ulong | integer/string | Bignum |
| <a name="sfixed32" /> sfixed32 | Always four bytes. | int32 | int | int | int32 | int | integer | Bignum or Fixnum (as required) |
| <a name="sfixed64" /> sfixed64 | Always eight bytes. | int64 | long | int/long | int64 | long | integer/string | Bignum |
| <a name="bool" /> bool |  | bool | boolean | boolean | bool | bool | boolean | TrueClass/FalseClass |
| <a name="string" /> string | A string must always contain UTF-8 encoded or 7-bit ASCII text. | string | String | str/unicode | string | string | string | String (UTF-8) |
| <a name="bytes" /> bytes | May contain any arbitrary sequence of bytes. | string | ByteString | str | []byte | ByteString | string | String (ASCII-8BIT) |
