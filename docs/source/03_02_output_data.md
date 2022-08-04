# Output data

## Neoantigens

NeoFox returns the neoantigen candidates and their annotated features as output. 
Two output formats are supported: tabular format or JSON format.
The user gets the neoantigen annotations in all formats.
Despite different structures, both formats provide the same content with the exception of the metadata on 
the annotations which is only present in the JSON format.

The following table describes each of the annotations in the output:  
  
**TABLE 1** 

| Column   Name                             |  Description                                                                                                                                                                                                                                                                                                    |  Feature group/ Paper             |
|-------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------|
| dnaVariantAlleleFrequency                 |  the variant allele frequency calculated from the DNA                                                                                                                                                                                                                                                           |  -                                |
| mutatedXmer                               |  the long mutated amino acid sequence                                                                                                                                                                                                                                                                           |  -                                |
| wildTypeXmer                              |  the long non-mutated amino acid sequence (when not provided in the input this will contain the Blastp closest sequence in the proteome)                                                                                                                                                                                                          |  -                                |
| patientIdentifier                         |  the patient identifier                                                                                                                                                                                                                                                                                         |  -                                |
| rnaExpression                             |  the RNA expression. If expression was imputed, this will will be `imputedGeneExpression`                                                                                                                                                                                                                       |  expression                       |
| imputedGeneExpression                     |  median gene expression in the TCGA cohort of the tumor entity provided in the patient file.                                                                                                                                                                                                                    |  expression                       |
| rnaVariantAlleleFrequency                 |  the variant allele frequency calculated from the RNA                                                                                                                                                                                                                                                           |  -                                |
| gene                                      |  the HGNC gene symbol                                                                                                                                                                                                                                                                                           |  -                                |
| Expression_mutated_transcript             |  transcript expression normalized by the variant allele frequency of the mutation                                                                                                                                                                                                                               |  expression                       |
| mutation_not_found_in_proteome            | indicates if mutated amino acid sequence was not found in the WT proteome by exact search                                                                                                                                                                                                                       |  Priority score                   |
| Best_rank_MHCI_score                      |  minimal MHC I binding rank score over all neoepitope candidates (8-11mers) and MHC I alleles                                                                                                                                                                                                                   |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_score_epitope              |  neoepitope candidate sequence with minimal MHC I binding rank score                                                                                                                                                                                                                                            |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_score_allele               |  the MHC I allele related to ` Best_rank_MHCI_score_epitope`                                                                                                                                                                                                                                                    |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_score                  |  minimal MHC I binding affinity  over all neoepitope candidates (8-11mers) and MHC I alleles                                                                                                                                                                                                                    |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_epitope                |  neoepitope candidate sequence with minimal MHC I binding affinity                                                                                                                                                                                                                                              |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_allele                 |  the MHC I allele related to ` Best_affinity_MHCI_epitope`                                                                                                                                                                                                                                                      |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_9mer_score                 |  minimal MHC I binding rank score over all neoepitope candidates (9mers only) and MHC I alleles                                                                                                                                                                                                                 |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_9mer_epitope               |  neoepitope candidate sequence (9mer) with minimal MHC I binding rank score                                                                                                                                                                                                                                     |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_9mer_allele                |  the MHC I allele related to ` Best_rank_MHCI_9mer_epitope`                                                                                                                                                                                                                                                     |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_9mer_score             |  minimal MHC I binding affinity over all neoepitope candidates (9mers) and MHC I alleles                                                                                                                                                                                                                        |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_9mer_allele            |  the MHC I allele related to ` Best_affinity_MHCI_9mer_epitope           `                                                                                                                                                                                                                                      |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_9mer_epitope           |  neoepitope candidate sequence (9mer) with minimal MHC I binding affinity                                                                                                                                                                                                                                       |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_score_WT               |  MHC I binding affinity  of `Best_affinity_MHCI_epitope_WT`                                                                                                                                                                                                                                                     |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_epitope_WT             |  WT epitope that corresponds to ` Best_affinity_MHCI_epitope`                                                                                                                                                                                                                                                   |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_allele_WT              |  the MHC I allele related to `Best_affinity_MHCI_epitope_WT`                                                                                                                                                                                                                                                    |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_score_WT                   | MHC I binding rank score of `Best_rank_MHCI_score_epitope_WT`                                                                                                                                                                                                                                                   |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_score_epitope_WT           |  WT epitope that corresponds to `Best_rank_MHCI_score_epitope`                                                                                                                                                                                                                                                  |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_score_allele_WT            |  the MHC I allele related to `Best_rank_MHCI_score_epitope_WT`                                                                                                                                                                                                                                                  |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_9mer_score_WT              | MHC I binding rank score of `Best_rank_MHCI_9mer_epitope_WT `                                                                                                                                                                                                                                                   |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_9mer_epitope_WT            |  WT epitope that corresponds to `Best_rank_MHCI_9mer_epitope`                                                                                                                                                                                                                                                   |  MHC I binding with netMHCpan     |
| Best_rank_MHCI_9mer_allele_WT             |  the MHC I allele related to `Best_rank_MHCI_9mer_epitope_WT `                                                                                                                                                                                                                                                  |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_9mer_score_WT          | MHC I binding affinity of ` Best_affinity_MHCI_9mer_allele_WT `                                                                                                                                                                                                                                                 |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_9mer_allele_WT         |  the MHC I allele related to ` Best_affinity_MHCI_9mer_epitope_WT`                                                                                                                                                                                                                                              |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_9mer_epitope_WT        |  WT epitope that corresponds to `Best_affinity_MHCI_9mer_epitope`                                                                                                                                                                                                                                               |  MHC I binding with netMHCpan     |
| PHBR-I                                    |  harmonic mean of minimal MHC I binding rank scores of all MHC I alleles of a patient                                                                                                                                                                                                                           |  PHBR-I                           |
| Best_affinity_MHCI_9mer_position_mutation |  indicates position of the mutation in ` Best_affinity_MHCI_9mer_epitope`                                                                                                                                                                                                                                       |  MHC I binding with netMHCpan     |
| Best_affinity_MHCI_9mer_anchor_mutated    |  mutation in ` Best_affinity_MHCI_9mer_epitope` in an anchor position (i.e. position 2 or 9)                                                                                                                                                                                                                    |  anchor/non-anchor                |
| Best_rank_MHCII_score                     |  minimal MHC II binding rank score over all neoepitope candidates (15mers) and all MHC II alleles                                                                                                                                                                                                               |  MHC II binding with netMHCIIpan  |
| Best_rank_MHCII_score_epitope             |  neoepitope candidate sequence with minimal MHC II binding rank score                                                                                                                                                                                                                                           |  MHC II binding with netMHCIIpan  |
| Best_rank_MHCII_score_allele              |  the MHC II isoform related to ` Best_rank_MHCII_score_epitope`                                                                                                                                                                                                                                                 |  MHC II binding with netMHCIIpan  |
| Best_affinity_MHCII_score                 |  minimal MHC II binding affinity  over all neoepitope candidates (15mers) and all MHC II alleles                                                                                                                                                                                                                |  MHC II binding with netMHCIIpan  |
| Best_affinity_MHCII_epitope               |  neoepitope candidate sequence with minimal MHC II binding affinity                                                                                                                                                                                                                                             |  MHC II binding with netMHCIIpan  |
| Best_affinity_MHCII_allele                |  the MHC II isoform related to ` Best_affinity_MHCII_epitope `                                                                                                                                                                                                                                                  |  MHC II binding with netMHCIIpan  |
| Best_rank_MHCII_score_WT                  |  minimal MHC II binding rank of  `Best_rank_MHCII_score_epitope_WT  `                                                                                                                                                                                                                                           |  MHC II binding with netMHCIIpan  |
| Best_rank_MHCII_score_epitope_WT          |  WT epitope sequence (15mer) that corresponds to ` Best_rank_MHCII_score_epitope `                                                                                                                                                                                                                              |  MHC II binding with netMHCIIpan  |
| Best_rank_MHCII_score_allele_WT           |  the MHC II isoform related to ` Best_rank_MHCII_score_epitope_WT`                                                                                                                                                                                                                                              |  MHC II binding with netMHCIIpan  |
| Best_affinity_MHCII_score_WT              |  minimal MHC II binding rank of  `Best_affinity_MHCII_epitope_WT`                                                                                                                                                                                                                                               |  MHC II binding with netMHCIIpan  |
| Best_affinity_MHCII_epitope_WT            |  WT epitope sequence (15mer) that corresponds to `  Best_affinity_MHCII_epitope`                                                                                                                                                                                                                                |  MHC II binding with netMHCIIpan  |
| Best_affinity_MHCII_allele_WT             |  the MHC II isoform related to ` Best_affinity_MHCII_epitope_WT`                                                                                                                                                                                                                                                |  MHC II binding with netMHCIIpan  |
| PHBR-II                                   |  harmonic mean of minimal MHC II binding rank scores of all MHC II alleles of a patient                                                                                                                                                                                                                         |  PHBR-II                          |
| Amplitude_MHCI_affinity_9mer              |  ratio of  `Best_affinity_MHCI_9mer_score_WT` and   `Best_affinity_MHCI_9mer_score`                                                                                                                                                                                                                             |  Recognition Potential            |
| Amplitude_MHCI_affinity                   |  ratio of   `Best_affinity_MHCI_score_WT` and `Best_affinity_MHCI_score`                                                                                                                                                                                                                                        |  Generator rate                   |
| Amplitude_MHCII_rank                      |  ratio of   `Best_rank_MHCII_score_WT` and `Best_rank_MHCII_score` and                                                                                                                                                                                                                                          |  Generator rate                   |
| Pathogensimiliarity_MHCI_9mer             |  score representing the   similarity of    `Best_affinity_MHCI_9mer_epitope` to pathogen sequences in IEDB   database                                                                                                                                                                                           |  Recognition Potential            |
| Pathogensimiliarity_MHCII                 |  score representing the   similarity of    `Best_affinity_MHCII_epitope` to pathogen sequences in IEDB   database                                                                                                                                                                                               |  Recognition Potential            |
| Recognition_Potential_MHCI_9mer           |  product of   `Amplitude_MHCI_affinity_9mer` and `Pathogensimiliarity_MHCI_affinity_9mer`                                                                                                                                                                                                                       |  Recognition Potential            |
| DAI_MHCI_affinity                         |  difference of   `Best_affinity_MHCI_score_WT` and `Best_affinity_MHCI_score`                                                                                                                                                                                                                                   |  DAI                              |
| CDN_MHCI                                  |  `Best_affinity_MHCI_score` <   50 nM                                                                                                                                                                                                                                                                           |  Generator rate                   |
| ADN_MHCI                                  |  `Best_affinity_MHCI_score` <   5000 nM and `Amplitude_MHCI_affinity` > 10                                                                                                                                                                                                                                      |  Generator rate                   |
| CDN_MHCII                                 |  `Best_rank_MHCII_score` < 1                                                                                                                                                                                                                                                                                    |  Generator rate                   |
| ADN_MHCII                                 |  `Best_rank_MHCII_score` < 4   and `Amplitude_MHCII_rank` < 2                                                                                                                                                                                                                                                   |  Generator rate                   |
| Generator_rate_CDN_MHCI                   |  number of neoepitope candidates   with MHC I binding affinity < 50 nM per neoantigen canidate                                                                                                                                                                                                                  |  Generator rate                   |
| Generator_rate_ADN_MHCI                   |  number of neoepitope candidates  with MHC I binding affinity < 5000 nM per neoantigen canidate 10x better affinity in comparison to corresponding WT peptide                                                                                                                                                   |  Generator rate                   |
| Generator_rate_MHCI                       | sum of `Generator_rate_CDN_MHCI` and `Generator_rate_ADN_MHCI`                                                                                                                                                                                                                                                  |  Generator rate                   |
| Generator_rate_CDN_MHCII                  |  number of neoepitope candidates   with MHC II binding rank score < 1 per neoantigen canidate                                                                                                                                                                                                                   |  Generator rate                   |
| Generator_rate_ADN_MHCII                  |  number of neoepitope candidates  with MHC II binding rank score < 4 per neoantigen candidate 4x better rank in comparison to corresponding WT peptide                                                                                                                                                          |  Generator rate                   |
| Generator_rate_MHCII                      | sum of `Generator_rate_CDN_MHCII` and `Generator_rate_ADN_MHCII`                                                                                                                                                                                                                                                |  Generator rate                   |
| Tcell_predictor_score                     |  output score of T cell predictor   model                                                                                                                                                                                                                                                                       |  Tcell predictor                  |
| Improved_Binder_MHCI                      |  ratio of   `Best_rank_MHCI_score_WT` and `Best_rank_MHCI_score` > 1.2                                                                                                                                                                                                                                          |  self-similarity                  |
| Selfsimilarity_MHCI_conserved_binder      |  score representing the   similarity between `Best_rank_MHCI_score_epitope` and   `Best_affinity_MHCI_epitope_WT`   For conservered binder only                                                                                                                                                                 |  self-similarity                  |
| Selfsimilarity_MHCI                       |  score representing the   similarity between `Best_rank_MHCI_score_epitope` and   `Best_affinity_MHCI_epitope_WT`                                                                                                                                                                                               |  self-similarity                  |
| Selfsimilarity_MHCII                      |  score representing the   similarity between `Best_affinity_MHCII_epitope` and    Best_affinity_MHCII_epitope_WT`                                                                                                                                                                                               |  self-similarity                  |
| Number_of_mismatches_MCHI                 |  number of amino acids that do no   match between `Best_rank_MHCI_score_epitope` and   `Best_rank_MHCI_score_epitope_WT`                                                                                                                                                                                        |  Priority score                   |
| Priority_score                            |  combinatorial score of several   features such as MHC binding, expression and VAF                                                                                                                                                                                                                              |  Priority score                   |
| Neoag_immunogenicity                      |  output score of neoag model                                                                                                                                                                                                                                                                                    |  neoag                            |
| IEDB_Immunogenicity_MHCI                  |  IEDB Immunogenicity score  for ` Best_affinity_MHCI_epitope `                                                                                                                                                                                                                                                  |  IEDB Immunogenicity              |
| IEDB_Immunogenicity_MHCII                 |  IEDB Immunogenicity score   for ` Best_affinity_MHCII_epitope`                                                                                                                                                                                                                                                 |  IEDB Immunogenicity              |
| MixMHCpred_best_peptide                   |  MHC class I neoepitope candidate   sequence with maximum MixMHCpred score over all neoepitope canidates   (8-11mers) and MHC I alleles                                                                                                                                                                         |  MixMHCpred                       |
| MixMHCpred_best_score                     |  maximum MixMHCpred score over   all neoepitope canidates (8-11mers) and MHC I alleles                                                                                                                                                                                                                          |  MixMHCpred                       |
| MixMHCpred_best_rank                      |  rank that corresponds to   `MixMHCpred_best_score`                                                                                                                                                                                                                                                             |  MixMHCpred                       |
| MixMHCpred_best_allele                    |  the allele with maximum   MixMHCpred score                                                                                                                                                                                                                                                                     |  MixMHCpred                       |
| MixMHC2pred_best_peptide                  |  MHC class II neoepitope   candidate sequence with minimal MixMHC2pred score over all neoepitope   canidates (13-18mers) and MHC II alleles                                                                                                                                                                     |  MixMHC2pred                      |
| MixMHC2pred_best_rank                     |  minimal MixMHC2pred score over   all neoepitope canidates (13-18mers) and MHC II alleles                                                                                                                                                                                                                       |  MixMHC2pred                      |
| MixMHC2pred_best_allele                   |  the MHC II isoform with minimum   MixMHC2pred rank score                                                                                                                                                                                                                                                       |  MixMHC2pred                      |
| Dissimilarity_MHCI                        |  score reflecting the   dissimilarity of `Best_affinity_MHCI_epitope` to the self-proteome                                                                                                                                                                                                                      |  dissimilarity                    |
| Dissimilarity_MHCII                       |  score reflecting the   dissimilarity of `Best_affinity_MHCII_epitope` to the self-proteome                                                                                                                                                                                                                     |  dissimilarity                    |
| vaxrank_binding_score                     |  total binding score of vaxrank                                                                                                                                                                                                                                                                                 |  vaxrank                          |
| vaxrank_total_score                       |  product of total binding score   and expression score. Originally, the root of the number of reads   supporting the mutation are used in the original implementation. To simplify,   the expression normalised to VAF is used.                                                                                 |  vaxrank                          |
| PRIME_best_allele                         | best predicted MHC allele by PRIME model                                                                                                                                                                                                                                                                        | PRIME                             |
| PRIME_best_peptide                        | best predicted neoepitope candidate by PRIME model                                                                                                                                                                                                                                                              | PRIME                             |
| PRIME_best_rank                           | output rank score of PRIME model                                                                                                                                                                                                                                                                                | PRIME                             |
| PRIME_best_score                          | output score of PRIME model                                                                                                                                                                                                                                                                                     | PRIME                             |
| Hex_alignment_score_MHCI                  | the alignment score by HEX for ` Best_affinity_MHCI_epitope `                                                                                                                                                                                                                                                   | HEX                               |
| Hex_alignment_score_MHCII                 | the alignment score by HEX for ` Best_affinity_MHCII_epitope`                                                                                                                                                                                                                                                   | HEX                               |

### Tabular format

An output table with the suffix "*_neoantigen_candidates_annotated.tsv*" is created.  
This table contains the neoantigen candidates information, the neoantigen annotations and if some user-specific additional columns  
were provided in the input table, these external annotations.

This is a dummy example:  

| dnaVariantAlleleFrequency | gene  | imputedGeneExpression  | mutatedXmer                 | position  | wildTypeXmer                | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency | +-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) | ADN_MHCI | ADN_MHCII | Amplitude_MHCII_rank | Amplitude_MHCI_affinity | Amplitude_MHCI_affinity_9mer | Best_affinity_MHCII_allele | Best_affinity_MHCII_allele_WT | Best_affinity_MHCII_epitope | Best_affinity_MHCII_epitope_WT | Best_affinity_MHCII_score | Best_affinity_MHCII_score_WT | Best_affinity_MHCI_9mer_allele | Best_affinity_MHCI_9mer_allele_WT | Best_affinity_MHCI_9mer_anchor_mutated | Best_affinity_MHCI_9mer_epitope | Best_affinity_MHCI_9mer_epitope_WT | Best_affinity_MHCI_9mer_position_mutation | Best_affinity_MHCI_9mer_score | Best_affinity_MHCI_9mer_score_WT | Best_affinity_MHCI_allele | Best_affinity_MHCI_allele_WT | Best_affinity_MHCI_epitope | Best_affinity_MHCI_epitope_WT | Best_affinity_MHCI_score | Best_affinity_MHCI_score_WT | Best_rank_MHCII_score | Best_rank_MHCII_score_WT | Best_rank_MHCII_score_allele | Best_rank_MHCII_score_allele_WT | Best_rank_MHCII_score_epitope | Best_rank_MHCII_score_epitope_WT | Best_rank_MHCI_9mer_allele | Best_rank_MHCI_9mer_allele_WT | Best_rank_MHCI_9mer_epitope | Best_rank_MHCI_9mer_epitope_WT | Best_rank_MHCI_9mer_score | Best_rank_MHCI_9mer_score_WT | Best_rank_MHCI_score | Best_rank_MHCI_score_WT | Best_rank_MHCI_score_allele | Best_rank_MHCI_score_allele_WT | Best_rank_MHCI_score_epitope | Best_rank_MHCI_score_epitope_WT | CDN_MHCI | CDN_MHCII | DAI_MHCI_affinity_cutoff500nM | Dissimilarity_MHCI_cutoff500nM | Expression_mutated_transcript | Generator_rate | IEDB_Immunogenicity_MHCI_cutoff500nM | Improved_Binder_MHCI | MixMHC2pred_best_allele | MixMHC2pred_best_peptide | MixMHC2pred_best_rank | MixMHCpred_best_allele | MixMHCpred_best_peptide | MixMHCpred_best_rank | MixMHCpred_best_score | Neoag_immunogenicity | Number_of_mismatches_MCHI | PHBR-I  | PHBR-II | Pathogensimiliarity_MHCI_affinity_9mer | Priority_score | Recognition_Potential_MHCI_affinity_9mer | Selfsimilarity_MHCI_conserved_binder | Tcell_predictor_score_cutoff500nM | VAF_in_RNA | VAF_in_tumor | [WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) | mutation_not_found_in_proteome | patient | substitution | transcript_expression | vaxrank_binding_score | vaxrank_total_score |
|---------------------------|-------|------------------------|-----------------------------|-----------|-----------------------------|-------------------|---------------|---------------------------|----------------------------------------|----------|-----------|----------------------|-------------------------|------------------------------|----------------------------|-------------------------------|-----------------------------|--------------------------------|---------------------------|------------------------------|--------------------------------|-----------------------------------|----------------------------------------|---------------------------------|------------------------------------|-------------------------------------------|-------------------------------|----------------------------------|---------------------------|------------------------------|----------------------------|-------------------------------|--------------------------|-----------------------------|-----------------------|--------------------------|------------------------------|---------------------------------|-------------------------------|----------------------------------|----------------------------|-------------------------------|-----------------------------|--------------------------------|---------------------------|------------------------------|----------------------|-------------------------|-----------------------------|--------------------------------|------------------------------|---------------------------------|----------|-----------|-------------------------------|--------------------------------|-------------------------------|----------------|--------------------------------------|----------------------|-------------------------|--------------------------|-----------------------|------------------------|-------------------------|----------------------|-----------------------|----------------------|---------------------------|---------|---------|----------------------------------------|----------------|------------------------------------------|--------------------------------------|-----------------------------------|------------|--------------|---------------------------------------------|--------------------------------|---------|--------------|-----------------------|-----------------------|---------------------|
 | 0.294                     | BRCA2 | 0.5                    | AAAAAAAAAAAAAFAAAAAAAAAAAAA | 14        | AAAAAAAAAAAAALAAAAAAAAAAAAA | Ptx               | 0.51950689    | 0.857                     | AAAAAAAAAAAAAFAAAAAAAAAAAAA            | 0        | 1         | 28                   | 0.88723                 | 0.88723                      | HLA-DQA10401-DQB10402      | HLA-DQA10401-DQB10402         | AAAAFAAAAAAAAAA             | AAAALAAAAAAAAAA                | 251.77                    | 513.02                       | HLA-C*16:01                    | HLA-C*16:01                       | 1                                      | AAAAAAAAF                       | AAAAAAAAL                          | 9                                         | 24.3                          | 21.7                             | HLA-C*16:01               | HLA-C*16:01                  | AAAAAAAAF                  | AAAAAAAAL                     | 24.3                     | 21.7                        | 0.05                  | 1.4                      | HLA-DQA10301-DQB10402        | HLA-DQA10301-DQB10402           | AAAAFAAAAAAAAAA               | AAAALAAAAAAAAAA                  | HLA-C*16:01                | HLA-C*16:01                   | AAAAAAAAF                   | AAAAAAAAL                      | 0.0592                    | 0.0493                       | 0.0592               | 0.0493                  | HLA-C*16:01                 | HLA-C*16:01                    | AAAAAAAAF                    | AAAAAAAAL                       | 1        | 1         | -2.6                          | 1                              | 0.44522                       | 1              | 0.18288                              | 0                    | DPA1_01_03__DPB1_04_01  | AAAAFAAAAAAAAAAA         | 0.997                 | B0702                  | AAAAAAAAF               | 0.1                  | 0.50487               | 13.16998             | 1                         | 0.31193 | 0.21892 | 0                                      | 0.07017        | 0                                        | 0.99178271                           | 0.40327581                        | 0.857      | 0.294        | AAAAAAAAAAAAALAAAAAAAAAAAAA                 | 1                              | Ptx     | I547T        | 0.51950689            | 3.7689                | 1.678               |
| 0.173                     | BRCA2 | 0.5                    | AAAAAAAAAAAAAMAAAAAAAAAAAAA | 14        | AAAAAAAAAAAAARAAAAAAAAAAAAA | Ptx               | 0.71575659    | 0.556                     | AAAAAAAAAAAAAMAAAAAAAAAAAAA            | 1        | 1         | 10                   | 90.685                  | 90.685                       | HLA-DQA10401-DQB10402      | HLA-DQA10401-DQB10402         | AAAAAAAAAMAAAAA             | AAAAAAAAARAAAAA                | 421.53                    | 554.92                       | HLA-C*16:01                    | HLA-C*16:01                       | 1                                      | AAAAAAAAM                       | AAAAAAAAR                          | 9                                         | 24.1                          | 6346.9                           | HLA-C*16:01               | HLA-C*16:01                  | AAAAAAAAM                  | AAAAAAAAR                     | 24.1                     | 6346.9                      | 0.25                  | 2.5                      | HLA-DQA10401-DQB10302        | HLA-DQA10401-DQB10302           | AAAAAAAAAAMAAAA               | AAAAAAAAAARAAAA                  | HLA-C*16:01                | HLA-C*16:01                   | AAAAAAAAM                   | AAAAAAAAR                      | 0.0587                    | 8.9317                       | 0.0587               | 8.9317                  | HLA-C*16:01                 | HLA-C*16:01                    | AAAAAAAAM                    | AAAAAAAAR                       | 1        | 1         | 6322.8                        | 1                              | 0.39796                       | 1              | 0.18288                              | 1                    | DPA1_01_03__DPB1_04_01  | AAAAMAAAAAAAAAAA         | 2.44                  | B0702                  | AAAAAAAAM               | 0.07                 | 0.5444                | 39.51379             | 1                         | 0.29303 | 1.5594  | 0                                      | 0.10626        | 0                                        | NA                                   | 0.46452844                        | 0.556      | 0.173        | AAAAAAAAAAAAARAAAAAAAAAAAAA                 | 1                              | Ptx     | E135S        | 0.71575659            | 3.8741                | 1.5417              |


### JSON format

An output file with the suffix "*_neoantigen_candidates_annotated.json*" is created.  
This file contains neoantigen candidates information in JSON format.  
Furthermore, a second file with the suffix *"_neoantigen_features.json"* is created.  
This file contains the annotated neoantigen features in JSON format.  
The names within the models are described in **TABLE 1**.

This is a dummy example of a "*_neoantigen_candidates.json*" file.  
This file contains a list of neoantigen candidate models (for further information, please see [here](05_models.md).  
To simplify, only one full neoantigen candidate model is shown:
```json
[{
    "patient_identifier": "Ptx",
    "gene": "BRCA2",
    "position": [14], 
    "wild_type_xmer": "AAAAAAAAAAAAALAAAAAAAAAAAAA",
    "mutated_xmer": "AAAAAAAAAAAAAFAAAAAAAAAAAAA",
    "rna_expression": 0.5195068939999999,
    "imputed_gene_expression": 0.5,
    "dna_variant_allele_frequency": 0.294,
    "rna_variant_allele_frequency": 0.857,
    "neofox_annotations": [...],
    "external_annotations": [...]
}, {
    "patient_identifier": "Ptx",
    "gene": "BRCA2",
    "position": [14], 
    "wild_type_xmer": "AAAAAAAAAAAAARAAAAAAAAAAAAA",
    "mutated_xmer": "AAAAAAAAAAAAAMAAAAAAAAAAAAA",
    "rna_expression": 0.715756594,
    "imputed_gene_expression": 0.5,
    "dna_variant_allele_frequency": 0.17300000000000001,
    "rna_variant_allele_frequency": 0.556,
    "neofox_annotations": [ ... ],
    "external_annotations": [ ... ]
}]
```

Notice that for simplicity purposes both fields `neofox_annotations` and `external_annotations` are not shown above. 
For further information, please see [here](05_models.md).

This is a dummy example of the field `neofox_annotations`.

```json
{
    "annotations": [{
        "name": "Best_rank_MHCI_score",
        "value": "0.0592"
    }, {
        "name": "Best_rank_MHCI_score_epitope",
        "value": "AAAAAAAAF"
    }, {
        "name": "Best_rank_MHCI_score_allele",
        "value": "HLA-C*16:01"
    }, {
        "name": "Best_affinity_MHCI_score",
        "value": "24.3"
    }, {
        "name": "Best_affinity_MHCI_epitope",
        "value": "AAAAAAAAF"
    }, {
        "name": "Best_affinity_MHCI_allele",
        "value": "HLA-C*16:01"
    }, {
        "name": "Best_rank_MHCI_9mer_score",
        "value": "0.0592"
    }, {
        "name": "Best_rank_MHCI_9mer_epitope",
        "value": "AAAAAAAAF"
    }, {
        "name": "Best_rank_MHCI_9mer_allele",
        "value": "HLA-C*16:01"
    }, {
        "name": "Best_affinity_MHCI_9mer_score",
        "value": "24.3"
    }, {
        "name": "Best_affinity_MHCI_9mer_allele",
        "value": "HLA-C*16:01"
    }, {
        "name": "Best_affinity_MHCI_9mer_epitope",
        "value": "AAAAAAAAF"
    }, {
        "name": "Best_affinity_MHCI_score_WT",
        "value": "21.7"
    }, {
        "name": "Best_affinity_MHCI_epitope_WT",
        "value": "AAAAAAAAL"
    }, {
        "name": "Best_affinity_MHCI_allele_WT",
        "value": "HLA-C*16:01"
    }, {
        "name": "Best_rank_MHCI_score_WT",
        "value": "0.0493"
    }, {
        "name": "Best_rank_MHCI_score_epitope_WT",
        "value": "AAAAAAAAL"
    }, {
        "name": "Best_rank_MHCI_score_allele_WT",
        "value": "HLA-C*16:01"
    }, {
        "name": "Best_rank_MHCI_9mer_score_WT",
        "value": "0.0493"
    }, {
        "name": "Best_rank_MHCI_9mer_epitope_WT",
        "value": "AAAAAAAAL"
    }, {
        "name": "Best_rank_MHCI_9mer_allele_WT",
        "value": "HLA-C*16:01"
    }, {
        "name": "Best_affinity_MHCI_9mer_score_WT",
        "value": "21.7"
    }, {
        "name": "Best_affinity_MHCI_9mer_allele_WT",
        "value": "HLA-C*16:01"
    }, {
        "name": "Best_affinity_MHCI_9mer_epitope_WT",
        "value": "AAAAAAAAL"
    }, {
        "name": "Generator_rate",
        "value": "1"
    }, {
        "name": "PHBR-I",
        "value": "0.31193"
    }, {
        "name": "Best_affinity_MHCI_9mer_position_mutation",
        "value": "9"
    }, {
        "name": "Best_affinity_MHCI_9mer_anchor_mutated",
        "value": "1"
    }, {
        "name": "Best_rank_MHCII_score",
        "value": "0.05"
    }, {
        "name": "Best_rank_MHCII_score_epitope",
        "value": "AAAAFAAAAAAAAAA"
    }, {
        "name": "Best_rank_MHCII_score_allele",
        "value": "HLA-DQA10301-DQB10402"
    }, {
        "name": "Best_affinity_MHCII_score",
        "value": "251.77"
    }, {
        "name": "Best_affinity_MHCII_epitope",
        "value": "AAAAFAAAAAAAAAA"
    }, {
        "name": "Best_affinity_MHCII_allele",
        "value": "HLA-DQA10401-DQB10402"
    }, {
        "name": "Best_rank_MHCII_score_WT",
        "value": "1.4"
    }, {
        "name": "Best_rank_MHCII_score_epitope_WT",
        "value": "AAAALAAAAAAAAAA"
    }, {
        "name": "Best_rank_MHCII_score_allele_WT",
        "value": "HLA-DQA10301-DQB10402"
    }, {
        "name": "Best_affinity_MHCII_score_WT",
        "value": "513.02"
    }, {
        "name": "Best_affinity_MHCII_epitope_WT",
        "value": "AAAALAAAAAAAAAA"
    }, {
        "name": "Best_affinity_MHCII_allele_WT",
        "value": "HLA-DQA10401-DQB10402"
    }, {
        "name": "PHBR-II",
        "value": "0.21892"
    }, {
        "name": "MixMHCpred_best_peptide",
        "value": "AAAAAAAAF"
    }, {
        "name": "MixMHCpred_best_score",
        "value": "0.50487"
    }, {
        "name": "MixMHCpred_best_rank",
        "value": "0.1"
    }, {
        "name": "MixMHCpred_best_allele",
        "value": "B0702"
    }, {
        "name": "MixMHC2pred_best_peptide",
        "value": "AAAAFAAAAAAAAAAA"
    }, {
        "name": "MixMHC2pred_best_rank",
        "value": "0.997"
    }, {
        "name": "MixMHC2pred_best_allele",
        "value": "DPA1_01_03__DPB1_04_01"
    }, {
        "name": "Expression_mutated_transcript",
        "value": "0.44522"
    }, {
        "name": "mutation_not_found_in_proteome",
        "value": "1"
    }, {
        "name": "Amplitude_MHCI_affinity_9mer",
        "value": "0.88723"
    }, {
        "name": "Amplitude_MHCI_affinity",
        "value": "0.88723"
    }, {
        "name": "Amplitude_MHCII_rank",
        "value": "28"
    }, {
        "name": "Pathogensimiliarity_MHCI_affinity_9mer",
        "value": "0"
    }, {
        "name": "Recognition_Potential_MHCI_affinity_9mer",
        "value": "0"
    }, {
        "name": "DAI_MHCI_affinity_cutoff500nM",
        "value": "-2.6"
    }, {
        "name": "CDN_MHCI",
        "value": "1"
    }, {
        "name": "ADN_MHCI",
        "value": "0"
    }, {
        "name": "CDN_MHCII",
        "value": "1"
    }, {
        "name": "ADN_MHCII",
        "value": "1"
    }, {
        "name": "Tcell_predictor_score_cutoff500nM",
        "value": "0.4032758100297953"
    }, {
        "name": "Improved_Binder_MHCI",
        "value": "0"
    }, {
        "name": "Selfsimilarity_MHCI_conserved_binder",
        "value": "0.9917827053614943"
    }, {
        "name": "Number_of_mismatches_MCHI",
        "value": "1"
    }, {
        "name": "Priority_score",
        "value": "0.07017"
    }, {
        "name": "Neoag_immunogenicity",
        "value": "13.16998"
    }, {
        "name": "IEDB_Immunogenicity_MHCI_cutoff500nM",
        "value": "0.18288"
    }, {
        "name": "Dissimilarity_MHCI_cutoff500nM",
        "value": "1"
    }, {
        "name": "vaxrank_binding_score",
        "value": "3.7689"
    }, {
        "name": "vaxrank_total_score",
        "value": "1.678"
    }, {
        "name": "patient",
        "value": "Ptx"
    }, {
        "name": "substitution",
        "value": "I547T"
    }, {
        "name": "transcript_expression",
        "value": "0.5195068939999999"
    }, {
        "name": "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)",
        "value": "AAAAAAAAAAAAAFAAAAAAAAAAAAA"
    }, {
        "name": "[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)",
        "value": "AAAAAAAAAAAAALAAAAAAAAAAAAA"
    }, {
        "name": "VAF_in_tumor",
        "value": "0.294"
    }, {
        "name": "VAF_in_RNA",
        "value": "0.857"
    }, {
        "name": "Unnamed: 8",
        "value": null
    }],
    "annotator": "Neofox",
    "annotator_version": "0.4.0",
    "timestamp": "20201211115212061465"
}
```

And this is a dummy example of the field `external_annotations`:
```json
[
  {
        "name": "external_annotation_1",
        "value": "0.857"
    }, {
        "name": "external_annotation_2",
        "value": "that"
    }
]
```

The metadata on the annotations will look as follows:
```json
{
 annotator: "NeoFox",
 annotator_version: "0.5.4",
 timestamp: "20211104154037935536",
 resources: [
  {
   name: "netMHCpan",
   version: "4.1"
  },
  {
   name: "netMHCIIpan",
   version: "4.0"
  },
  {
   name: "mixMHCpred",
   version: "2.1"
  },
  {
   name: "mixMHC2pred",
   version: "1.2"
  },
  {
   name: "IEDB",
   url: "http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip",
   hash: "d225ab671ef375400d387354a5f450ff",
   download_timestamp: "20211103221051"
  },
  {
   name: "Human Ensembl proteome",
   version: "2021_03",
   url: "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz",
   hash: "a41303fd38380ca0321cf8a3d9beb4bc",
   download_timestamp: "20211103221051"
  },
  {
   name: "Mouse Ensembl proteome",
   version: "2021_03",
   url: "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000589/UP000000589_10090.fasta.gz",
   hash: "27a5de8c1eca42eebaf56400945cf7cb",
   download_timestamp: "20211103221051"
  },
  {
   name: "IMGT/HLA database",
   version: "3.46.0",
   url: "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist.txt",
   hash: "5a7618819498b88d0790bf4d58975d13",
   download_timestamp: "20211103221051"
  }
 ]
}
```


## Neoepitopes

NeoFox returns the neoepitopes candidates and their annotated features as output when using  
the flag `--with-all-neoepitopes` or when explicitly annotating neoepitopes with the command `neofox-epitope`.
Two output formats are supported: tabular format or JSON format.
The user gets the neoepitope annotations in all formats.
Despite different structures, both formats provide the same content with the exception of the metadata on
the annotations which is only present in the JSON format.
The tabular format is split into two tables:  
a first one for the MHC-I neoepitope candidates and a second one for  
the MHC-II neoepitope candidates.

### Tabular format

Two output files with the suffix "*_mhcI_epitope_candidates_annotated.tsv" and ""*_mhcII_epitope_candidates_annotated.tsv"" are created.

The following table describes each of the annotations in the output:

**TABLE 2**

| Column   Name                   | Description                                                                                                                        | Feature group/ Paper                                    |
|---------------------------------|------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| dnaVariantAlleleFrequency       | the variant allele frequency calculated from the DNA                                                                               | -                                                       |
| mutatedSequence                 | the mutated amino acid sequence                                                                                                    | -                                                       |
| wildTypeSequence                | the non-mutated amino acid sequence (when not provided in the input this will contain the Blastp closest sequence in the proteome) | -                                                       |
| alleleMhcI / isoformMhcII       | either the MHC-I allele for MHC-I neoepitopes or the MHC-II isoform for MHC-II neoepitopes                                         | -                                                       |
| patientIdentifier               | the patient identifier (optional)                                                                                                  | -                                                       |
| rnaExpression                   | the RNA expression. If expression was imputed, this will will be `imputedGeneExpression`                                           | expression                                              |
| imputedGeneExpression           | median gene expression in the TCGA cohort of the tumor entity provided in the patient file.                                        | expression                                              |
| rnaVariantAlleleFrequency       | the variant allele frequency calculated from the RNA                                                                               | -                                                       |
| gene                            | the HGNC gene symbol                                                                                                               | -                                                       |
| affinityMutated                 | NetMHCpan / NetMHCIIpan affinity score of the mutated peptide for MHC-I or MHC-II neoepitopes respectively                         | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan     |
| affinityWildType                | NetMHCpan / NetMHCIIpan affinity score of the wild type peptide for MHC-I or MHC-II neoepitopes respectively                       | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan     |
| rankMutated                     | NetMHCpan / NetMHCIIpan rank of the mutated peptide for MHC-I or MHC-II neoepitopes respectively                                   | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan     |
| rankWildType                    | NetMHCpan / NetMHCIIpan rank of the wild type peptide for MHC-I or MHC-II neoepitopes respectively                                 | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan     |
| MixMHCpred_affinity_score       | MixMHCpred score of the mutated peptide for MHC-I neoepitopes                                                                      | MHC-I binding with mixMHCpred                           |
| MixMHCpred_rank                 | MixMHCpred rank of the mutated peptide for MHC-I neoepitopes                                                                       | MHC-I binding with mixMHCpred                           |
| MixMHCpred_WT_affinity_score    | MixMHCpred score of the wild type peptide for MHC-I neoepitopes                                                                    | MHC-I binding with mixMHCpred                           |
| MixMHCpred_WT_rank              | MixMHCpred rank of the wild type peptide for MHC-I neoepitopes                                                                     | MHC-I binding with mixMHCpred                           |
| MixMHC2pred_affinity_score      | MixMHC2pred score of the mutated peptide for MHC-II neoepitopes                                                                    | MHC-II binding with mixMHC2pred                         |
| MixMHC2pred_rank                | MixMHC2pred rank of the mutated peptide for MHC-II neoepitopes                                                                     | MHC-II binding with mixMHC2pred                         |
| MixMHC2pred_WT_affinity_score   | MixMHC2pred score of the wild type peptide for MHC-II neoepitopes                                                                  | MHC-II binding with mixMHC2pred                         |
| MixMHC2pred_WT_rank             | MixMHC2pred rank of the wild type peptide for MHC-II neoepitopes                                                                   | MHC-II binding with mixMHC2pred                         |
| PRIME_affinity_score            | PRIME score of the mutated peptide for MHC-I neoepitopes                                                                           | MHC-I binding with PRIME                                |
| PRIME_rank                      | PRIME rank of the mutated peptide for MHC-I neoepitopes                                                                            | MHC-I binding with PRIME                                |
| PRIME_WT_affinity_score         | PRIME score of the wild type peptide for MHC-I neoepitopes                                                                         | MHC-I binding with PRIME                                |
| PRIME_WT_rank                   | PRIME rank of the wild type peptide for MHC-I neoepitopes                                                                          | MHC-I binding with PRIME                                |
| DAI                             | difference of `affinityWildType` and `affinityMutated`                                                                             | DAI (only availble for MHC-I)                           |
| IEDB_Immunogenicity             | IEDB Immunogenicity score for `affinityMutated`                                                                                    | IEDB immunogenicity                                     |
| Improved_Binder_MHCI            | ratio of `affinityWildType` and `affinityMutated` > 1.2                                                                            | self-similarity (only available for MHC-I)              |
| Priority_score                  | combinatorial score of several features such as MHC binding, expression and VAF                                                    | Priority score                                          |
| mutation_not_found_in_proteome  | indicates if mutated amino acid sequence was not found in the WT proteome by exact search                                          | Priority score                                          |
| Selfsimilarity                  | score representing the similarity between `rankMutated` and `rankWildType`                                                         | self-similarity                                         |
| Selfsimilarity_conserved_binder | score representing the similarity between `rankMutated` and `rankWildType` for conserved binder only                               | self-similarity (only available for MHC-I)              |
| dissimilarity_score             | score reflecting the dissimilarity of `affinityMutated` to the self-proteome                                                       | dissimilarity                                           |
| Tcell_predictor_score           | output score of T cell predictor model                                                                                             | Tcell predictor (only available for MHC-I)              |
| amplitude                       | ratio of `affinityWildType` and `affinityMutated` for MHC-I and `rankWildType` and `rankMutated` for MHC-II                        | Generator rate                                          |
| anchor_mutated                  | flag indicating if a mutation lies in an anchor position (i.e. position 2 or 9)                                                    | anchor/non-anchor (only available for MHC-I)            |
| hex_alignment_score             | the alignment score by HEX for `mutatedSequence`                                                                                   | HEX                                                     |
| number_of_mismatches            | number of amino acids that do no match between `mutatedSequence` and `wildTypeSequence`                                            | Priority score (only available for MHC-I)               |
| pathogen_similarity             | score representing the similarity of `mutatedSequence` to pathogen sequences in IEDB database                                      | Recognition potential                                   |
| recognition_potential           | product of `amplitude` and `pathogenSimilarity`                                                                                    | Recognition potential (only available for MHC-I)        |
| position_mutation               | indicates position of the mutation in `mutatedSequence`                                                                            | MHC I binding with netMHCpan (only available for MHC-I) |


## JSON format

Only when using the command `neofox-epitope` an output file with the suffix "*_neoepitope_candidates_annotated.json*" is created.  
This file contains neoepitope candidates information in JSON format.  
The names within the models are described in **TABLE 2**.

This is a dummy example of a "*_neoantigen_candidates.json*" file.  
This file contains a list of neoantigen candidate models (for further information, please see [here](05_models.md).  
To simplify, only one full neoantigen candidate model is shown:
```json
[{
    "patient_identifier": "Ptx",
    "gene": "BRCA2",
    "mutated_peptide": "AAAALAAAA",
    "wild_type_peptide": "AAAAFAAAA",
    "allele_mhc_i": "HLA-A*01:01",
    "rna_expression": 0.519,
    "imputed_gene_expression": 0.5,
    "dna_variant_allele_frequency": 0.294,
    "rna_variant_allele_frequency": 0.857,
    "affinity_mutated": 2.567,
    "rank_mutated": 0.898,
    "affinity_wild_type": 1.023,
    "rank_wild_type": 2.398,
    "neofox_annotations": [...],
    "external_annotations": [...]
}, {
    "patient_identifier": "Ptx",
    "gene": "BRCA2",
    "mutated_peptide": "AAAAAAAAAAAAARAAAAAAAAAAAAA",
    "wild_type_peptide": "AAAAAAAAAAAAAMAAAAAAAAAAAAA",
    "isoform_mhc_i_i": "DRB1*01:01",
    "rna_expression": 0.715,
    "imputed_gene_expression": 0.5,
    "dna_variant_allele_frequency": 0.173,
    "rna_variant_allele_frequency": 0.556,
    "affinity_mutated": 2.567,
    "rank_mutated": 0.898,
    "affinity_wild_type": 1.023,
    "rank_wild_type": 2.398,
    "neofox_annotations": [ ... ],
    "external_annotations": [ ... ]
}]
```

Notice that for simplicity purposes both fields `neofox_annotations` and `external_annotations` are not shown above.
For further information, please see [here](05_models.md).  
For an example of the NeoFox annotations section, see the previous section.
