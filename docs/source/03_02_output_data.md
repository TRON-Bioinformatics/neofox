# Output data

## General information

NeoFox returns the neoantigen candidates and their annotated features as output. Several output formats are supported: [short-wide](#short-wide-format ) or [tall-skinny](#tall-skinny-format) tabular formats or [json](#json-format) format. The user can choose one preferred format or get the neoantigen annotations in all formats. Despite different structures, all three formats provide the same content.

The following table describes each of the annotations in the output:  
  
**TABLE 1** 

| Column   Name                             | Description                                                                                                                                                                                                                                      | Feature group/ Paper            |
|-------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------|
| identifier                                | unique neoantigen id given by   NeoFox                                                                                                                                                                                                           | -                               |
| dnaVariantAlleleFrequency                 | the variant allele frequency   calculated from the DNA                                                                                                                                                                                           | -                               |
| mutation.mutatedXmer                      | the mutated amino acid sequence   defined by left flanking region, mutated amino acid and right flanking region                                                                                                                                  | -                               |
| mutation.wildTypeXmer                     | the non-mutated amino acid   sequence defined by left flanking region, wildtype amino acid and right   flanking region                                                                                                                           | -                               |
| patientIdentifier                         | the patient identifier                                                                                                                                                                                                                           | -                               |
| rnaExpression                             | the RNA expression. If   expression was imputed, the value is the same as `imputedGeneExpression`                                                                                                                                                     | expression                      |
| imputedGeneExpression                     | median gene expression in the TCGA cohort of the entity provided in the patient file.                                                                                                                                                      | expression                      |
| rnaVariantAlleleFrequency                 | the variant allele frequency   calculated from the RNA                                                                                                                                                                                           | -                               |
| gene                                      | the HGNC gene symbol                                                                                                                                                                                                                             | -                               |
| Expression_mutated_transcript             | transcript expression normalized   by the variant allele frequency of the mutation                                                                                                                                                               | expression                      |
| mutation_not_found_in_proteome            | mutated amino acid sequence not   found in the WT proteome by exact search                                                                                                                                                                       | Priority score                  |
| Best_rank_MHCI_score                      | minimal MHC I binding rank score   over all neoepitope candidate (8-11mers) and MHC I alleles                                                                                                                                                    | MHC I binding with netMHCpan    |
| Best_rank_MHCI_score_epitope              | neoepitope candidate sequence   with minimal MHC I binding rank score                                                                                                                                                                            | MHC I binding with netMHCpan    |
| Best_rank_MHCI_score_allele               | the allele with minimal MHC I   binding rank score                                                                                                                                                                                               | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_score                  | minimal MHC I binding   affinity  over all neoepitope candidate   (8-11mers) and MHC I alleles                                                                                                                                                   | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_epitope                | neoepitope candidate sequence   with minimal MHC I binding affinity                                                                                                                                                                              | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_allele                 | the allele with minimal MHC I   binding affinity                                                                                                                                                                                                 | MHC I binding with netMHCpan    |
| Best_rank_MHCI_9mer_score                 | minimal MHC I binding rank score   over all neoepitope candidate (9mers only) and MHC I alleles                                                                                                                                                  | MHC I binding with netMHCpan    |
| Best_rank_MHCI_9mer_epitope               | neoepitope candidate sequence   (9mer) with minimal MHC I binding rank score                                                                                                                                                                     | MHC I binding with netMHCpan    |
| Best_rank_MHCI_9mer_allele                | the allele with minimal MHC I   binding rank score of 9mer neoepitope candidate                                                                                                                                                                  | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_9mer_score             | minimal MHC I binding affinity   over all neoepitope candidate(9mers) and MHC I alleles                                                                                                                                                          | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_9mer_allele            | neoepitope candidate sequence   (9mer) with minimal MHC I binding affinity                                                                                                                                                                       | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_9mer_epitope           | the allele with minimal MHC I   binding affinity of 9mer neoepitope candidate                                                                                                                                                                    | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_score_WT               | minimal MHC I binding affinity   of WT epitope sequence that corresponds to the general (8-11mers) best   neoepitope candidate based on MHC I binding affinity                                                                                   | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_epitope_WT             | WT epitope with minimal MHC I   binding affinity that corresponds to the general (8-11mers) best neoepitope   candidate based on MHC I binding affinity                                                                                          | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_allele_WT              | the allele with minimal WT MHC I   binding affinity                                                                                                                                                                                              | MHC I binding with netMHCpan    |
| Best_rank_MHCI_score_WT                   | minimal MHC I binding rank score   of WT epitope sequence that corresponds to the general (8-11mers) best   neoepitope candidate based on MHC I binding rank score                                                                               | MHC I binding with netMHCpan    |
| Best_rank_MHCI_score_epitope_WT           | WT epitope with minimal MHC I   binding rank score that corresponds to the general (8-11mers) best neoepitope   candidate based on MHC I binding rank score                                                                                      | MHC I binding with netMHCpan    |
| Best_rank_MHCI_score_allele_WT            | the allele with minimal WT MHC I   binding rank score                                                                                                                                                                                            | MHC I binding with netMHCpan    |
| Best_rank_MHCI_9mer_score_WT              | minimal MHC I binding rank score   of 9mer WT epitope sequence that corresponds to the best 9mer neoepitope   candidate based on MHC I binding rank score. The WT epitope should bind to   the same MHC allele as the neoepitope candidate.      | MHC I binding with netMHCpan    |
| Best_rank_MHCI_9mer_epitope_WT            | 9mer WT epitope with minimal MHC   I binding rank score that corresponds to the best 9mer neoepitope candidate   based on MHC I binding rank score. The WT epitope should bind to the same MHC   allele as the neoepitope candidate.             | MHC I binding with netMHCpan    |
| Best_rank_MHCI_9mer_allele_WT             | the allele with minimal  MHC I binding rank score related to 9mer WT   epitope. The WT epitope should bind to the same MHC allele as the neoepitope   candidate.The WT epitope should bind to the same MHC allele as the neoepitope   candidate. | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_9mer_score_WT          | minimal MHC I binding affinity   of 9mer WT epitope sequence that corresponds to the best 9mer neoepitope   candidate based on MHC I binding affinity                                                                                            | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_9mer_allele_WT         | 9mer WT epitope with minimal MHC   I binding affinity that corresponds to the best 9mer neoepitope candidate   based on MHC I binding affinity                                                                                                   | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_9mer_epitope_WT        | the allele with minimal 9mer WT   epitope  MHC I binding affinity. The WT   epitope should bind to the same MHC allele as the neoepitope candidate.                                                                                              | MHC I binding with netMHCpan    |
| Generator_rate                            | number of neoepitope candidates   with MHC I binding affinity < 50 nM per neoantigen canidate                                                                                                                                                    | Generator rate                  |
| PHBR-I                                    | harmonic mean of minimal MHC I   binding rank scores of each MHC I alleles                                                                                                                                                                       | PHBR-I                          |
| Best_affinity_MHCI_9mer_position_mutation | position of the mutation in the   9mer neoepitope candidate with the minimal MHC I binding affinity                                                                                                                                              | MHC I binding with netMHCpan    |
| Best_affinity_MHCI_9mer_anchor_mutated    | mutation in   `Best_rank_MHCI_9mer_epitope` in an anchor position (i.e. position 2 or 9)                                                                                                                                                         | anchor/non-anchor               |
| Best_rank_MHCII_score                     | minimal MHC II binding rank   score over all neoepitope candidates (15mers) and all MHC II alleles                                                                                                                                               | MHC II binding with netMHCIIpan |
| Best_rank_MHCII_score_epitope             | neoepitope candidate sequence   (15mer) with minimal MHC II binding rank score                                                                                                                                                                   | MHC II binding with netMHCIIpan |
| Best_rank_MHCII_score_allele              | the allele with minimal MHC II   binding rank score                                                                                                                                                                                              | MHC II binding with netMHCIIpan |
| Best_affinity_MHCII_score                 | minimal MHC II binding   affinity  over all neoepitope   candidates (15mers) and all MHC II alleles                                                                                                                                              | MHC II binding with netMHCIIpan |
| Best_affinity_MHCII_epitope               | neoepitope candidate sequence   with best MHC II binding affinity                                                                                                                                                                                | MHC II binding with netMHCIIpan |
| Best_affinity_MHCII_allele                | the allele with minimal MHC II   binding affinity                                                                                                                                                                                                | MHC II binding with netMHCIIpan |
| Best_rank_MHCII_score_WT                  | minimal MHC II binding rank of   WT epitope that corresponds to the best neoepitope candidate based on MHC II   binding rank score                                                                                                               | MHC II binding with netMHCIIpan |
| Best_rank_MHCII_score_epitope_WT          | WT epitope sequence (15mer) with   minimal MHC II binding rank score that corresponds to the best neoepitope   candidate based on MHC II binding rank score                                                                                      | MHC II binding with netMHCIIpan |
| Best_rank_MHCII_score_allele_WT           | the allele with minimal WT MHC   II binding rank score                                                                                                                                                                                           | MHC II binding with netMHCIIpan |
| Best_affinity_MHCII_score_WT              | minimal MHC II binding   affinity  of the WT epitope that   corresponds to the best neoepitope candidate based on MHC II binding rank   score                                                                                                    | MHC II binding with netMHCIIpan |
| Best_affinity_MHCII_epitope_WT            | WT epitope sequence (15mer) with   minimal MHC II binding affinity that corresponds to the best neoepitope   candidate based on MHC II binding affinity                                                                                          | MHC II binding with netMHCIIpan |
| Best_affinity_MHCII_allele_WT             | the allele with minimal WT MHC   II binding affinity                                                                                                                                                                                             | MHC II binding with netMHCIIpan |
| PHBR-II                                   | harmonic mean of minimal MHC II   binding rank scores of each MHC II alleles                                                                                                                                                                     | PHBR-II                         |
| Amplitude_MHCI_affinity_9mer              | ratio of  `Best_affinity_MHCI_9mer_score_WT` and   `Best_affinity_MHCI_9mer_score`                                                                                                                                                               | Recognition Potential           |
| Amplitude_MHCI_affinity                   | ratio of   `Best_affinity_MHCI_score_WT` and `Best_affinity_MHCI_score`                                                                                                                                                                          | Generator rate                  |
| Amplitude_MHCII_rank                      | ratio of   `Best_rank_MHCII_score_WT` and `Best_rank_MHCII_score` and                                                                                                                                                                            | Generator rate                  |
| Pathogensimiliarity_MHCI_affinity_9mer    | score representing the   similarity of    `Best_affinity_MHCI_9mer_epitope` to pathongen sequences in IEDB   database                                                                                                                            | Recognition Potential           |
| Recognition_Potential_MHCI_affinity_9mer  | product of   `Amplitude_MHCI_affinity_9mer` and `Pathogensimiliarity_MHCI_affinity_9mer`                                                                                                                                                         | Recognition Potential           |
| DAI_MHCI_affinity_cutoff500nM             | difference of   `Best_affinity_MHCI_score_WT` and `Best_affinity_MHCI_score`                                                                                                                                                                     | DAI                             |
| CDN_MHCI                                  | `Best_affinity_MHCI_score` <   50 nM                                                                                                                                                                                                             | Generator rate                  |
| ADN_MHCI                                  | `Best_affinity_MHCI_score` <   5000 nM and `Amplitude_MHCI_affinity` > 10                                                                                                                                                                        | Generator rate                  |
| CDN_MHCII                                 | `Best_rank_MHCII_score` < 1                                                                                                                                                                                                                      | Generator rate                  |
| ADN_MHCII                                 | `Best_rank_MHCII_score` < 4   and `Amplitude_MHCII_rank` < 2                                                                                                                                                                                     | Generator rate                  |
| Tcell_predictor_score_cutoff500nM         | output score of T cell predictor   model                                                                                                                                                                                                         | Tcell predictor                 |
| Improved_Binder_MHCI                      | ratio of   `Best_rank_MHCI_score_WT` and `Best_rank_MHCI_score` > 1.2                                                                                                                                                                           | self-similarity                 |
| Selfsimilarity_MHCI_conserved_binder      | score representing the   similarity between `Best_rank_MHCI_score_epitope` and   `Best_affinity_MHCI_epitope_WT`                                                                                                                                 | self-similarity                 |
| Number_of_mismatches_MCHI                 | number of amino acids that do no   match between `Best_rank_MHCI_score_epitope` and   `Best_rank_MHCI_score_epitope_WT`                                                                                                                          | Priority score                  |
| Priority_score                            | combinatorial Score of several   features such as MHC binding, expression and VAF                                                                                                                                                                | Priority score                  |
| Neoag_immunogenicity                      | output score of neoag model                                                                                                                                                                                                                      | neoag                           |
| IEDB_Immunogenicity_MHCI_cutoff500nM      | IEDB Immunogenicity score                                                                                                                                                                                                                        | IEDB Immunogenicity             |
| MixMHCpred_best_peptide                   | MHC class I neoepitope candidate   sequence with maximum MixMHCpred score over all neoepitope canidates   (8-11mers) and MHC I alleles                                                                                                           | MixMHCpred                      |
| MixMHCpred_best_score                     | maximum MixMHCpred score over   all neoepitope canidates (8-11mers) and MHC I alleles                                                                                                                                                            | MixMHCpred                      |
| MixMHCpred_best_rank                      | rank that corresponds to   `MixMHCpred_best_score`                                                                                                                                                                                               | MixMHCpred                      |
| MixMHCpred_best_allele                    | the allele with maximum   MixMHCpred score                                                                                                                                                                                                       | MixMHCpred                      |
| MixMHC2pred_best_peptide                  | MHC class II neoepitope   candidate sequence with minimal MixMHC2pred score over all neoepitope   canidates (13-18mers) and MHC II alleles                                                                                                       | MixMHC2pred                     |
| MixMHC2pred_best_rank                     | minimal MixMHC2pred score over   all neoepitope canidates (13-18mers) and MHC II alleles                                                                                                                                                         | MixMHC2pred                     |
| MixMHC2pred_best_allele                   | the allele with minimum   MixMHCpred score                                                                                                                                                                                                       | MixMHC2pred                     |
| Dissimilarity_MHCI_cutoff500nM            | score reflecting the   dissimilarity of `Best_affinity_MHCI_epitope` to the self-proteome                                                                                                                                                        | dissimilarity                   |
| vaxrank_binding_score                     | total binding score of vaxrank                                                                                                                                                                                                                   | vaxrank                         |
| vaxrank_total_score                       | product of total binding score   and expression score. Originally, the root of the number of reads   supporting the mutation are used in the original implementation. To simplify,   the expression normalised to VAF is used.            | vaxrank                         |

## Short-wide format

If the `--with-short-wide-table` flag is enabled an output table with the suffix "*_neoantigen_candidates_annotated.tsv*" is created. This table contains the neoantigen candidates information, the neoantigen annotations and if some user-specific additional columns were provided in the input table, these external annotations.  

This is a dummy example:  

| identifier               | dnaVariantAlleleFrequency | gene  |imputedGeneExpression	 | mutation.mutatedXmer        | mutation.position | mutation.wildTypeXmer       | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency | +-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) | ADN_MHCI | ADN_MHCII | Amplitude_MHCII_rank | Amplitude_MHCI_affinity | Amplitude_MHCI_affinity_9mer | Best_affinity_MHCII_allele | Best_affinity_MHCII_allele_WT | Best_affinity_MHCII_epitope | Best_affinity_MHCII_epitope_WT | Best_affinity_MHCII_score | Best_affinity_MHCII_score_WT | Best_affinity_MHCI_9mer_allele | Best_affinity_MHCI_9mer_allele_WT | Best_affinity_MHCI_9mer_anchor_mutated | Best_affinity_MHCI_9mer_epitope | Best_affinity_MHCI_9mer_epitope_WT | Best_affinity_MHCI_9mer_position_mutation | Best_affinity_MHCI_9mer_score | Best_affinity_MHCI_9mer_score_WT | Best_affinity_MHCI_allele | Best_affinity_MHCI_allele_WT | Best_affinity_MHCI_epitope | Best_affinity_MHCI_epitope_WT | Best_affinity_MHCI_score | Best_affinity_MHCI_score_WT | Best_rank_MHCII_score | Best_rank_MHCII_score_WT | Best_rank_MHCII_score_allele | Best_rank_MHCII_score_allele_WT | Best_rank_MHCII_score_epitope | Best_rank_MHCII_score_epitope_WT | Best_rank_MHCI_9mer_allele | Best_rank_MHCI_9mer_allele_WT | Best_rank_MHCI_9mer_epitope | Best_rank_MHCI_9mer_epitope_WT | Best_rank_MHCI_9mer_score | Best_rank_MHCI_9mer_score_WT | Best_rank_MHCI_score | Best_rank_MHCI_score_WT | Best_rank_MHCI_score_allele | Best_rank_MHCI_score_allele_WT | Best_rank_MHCI_score_epitope | Best_rank_MHCI_score_epitope_WT | CDN_MHCI | CDN_MHCII | DAI_MHCI_affinity_cutoff500nM | Dissimilarity_MHCI_cutoff500nM | Expression_mutated_transcript | Generator_rate | IEDB_Immunogenicity_MHCI_cutoff500nM | Improved_Binder_MHCI | MixMHC2pred_best_allele | MixMHC2pred_best_peptide | MixMHC2pred_best_rank | MixMHCpred_best_allele | MixMHCpred_best_peptide | MixMHCpred_best_rank | MixMHCpred_best_score | Neoag_immunogenicity | Number_of_mismatches_MCHI | PHBR-I  | PHBR-II | Pathogensimiliarity_MHCI_affinity_9mer | Priority_score | Recognition_Potential_MHCI_affinity_9mer | Selfsimilarity_MHCI_conserved_binder | Tcell_predictor_score_cutoff500nM | VAF_in_RNA | VAF_in_tumor | [WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL) | mutation_not_found_in_proteome | patient | substitution | transcript_expression | vaxrank_binding_score | vaxrank_total_score |
|--------------------------|---------------------------|-------|-----------------------------|-----------------------------|-------------------|-----------------------------|-------------------|---------------|---------------------------|----------------------------------------|----------|-----------|----------------------|-------------------------|------------------------------|----------------------------|-------------------------------|-----------------------------|--------------------------------|---------------------------|------------------------------|--------------------------------|-----------------------------------|----------------------------------------|---------------------------------|------------------------------------|-------------------------------------------|-------------------------------|----------------------------------|---------------------------|------------------------------|----------------------------|-------------------------------|--------------------------|-----------------------------|-----------------------|--------------------------|------------------------------|---------------------------------|-------------------------------|----------------------------------|----------------------------|-------------------------------|-----------------------------|--------------------------------|---------------------------|------------------------------|----------------------|-------------------------|-----------------------------|--------------------------------|------------------------------|---------------------------------|----------|-----------|-------------------------------|--------------------------------|-------------------------------|----------------|--------------------------------------|----------------------|-------------------------|--------------------------|-----------------------|------------------------|-------------------------|----------------------|-----------------------|----------------------|---------------------------|---------|---------|----------------------------------------|----------------|------------------------------------------|--------------------------------------|-----------------------------------|------------|--------------|---------------------------------------------|--------------------------------|---------|--------------|-----------------------|-----------------------|---------------------|
| ou11p7RD+tZvjY88DA55Mw== | 0.294                     | BRCA2 | 0.5| AAAAAAAAAAAAAFAAAAAAAAAAAAA| AAAAAAAAAAAAAFAAAAAAAAAAAAA | 14                | AAAAAAAAAAAAALAAAAAAAAAAAAA | Ptx               | 0.51950689    | 0.857                     | AAAAAAAAAAAAAFAAAAAAAAAAAAA            | 0        | 1         | 28                   | 0.88723                 | 0.88723                      | HLA-DQA10401-DQB10402      | HLA-DQA10401-DQB10402         | AAAAFAAAAAAAAAA             | AAAALAAAAAAAAAA                | 251.77                    | 513.02                       | HLA-C*16:01                    | HLA-C*16:01                       | 1                                      | AAAAAAAAF                       | AAAAAAAAL                          | 9                                         | 24.3                          | 21.7                             | HLA-C*16:01               | HLA-C*16:01                  | AAAAAAAAF                  | AAAAAAAAL                     | 24.3                     | 21.7                        | 0.05                  | 1.4                      | HLA-DQA10301-DQB10402        | HLA-DQA10301-DQB10402           | AAAAFAAAAAAAAAA               | AAAALAAAAAAAAAA                  | HLA-C*16:01                | HLA-C*16:01                   | AAAAAAAAF                   | AAAAAAAAL                      | 0.0592                    | 0.0493                       | 0.0592               | 0.0493                  | HLA-C*16:01                 | HLA-C*16:01                    | AAAAAAAAF                    | AAAAAAAAL                       | 1        | 1         | -2.6                          | 1                              | 0.44522                       | 1              | 0.18288                              | 0                    | DPA1_01_03__DPB1_04_01  | AAAAFAAAAAAAAAAA         | 0.997                 | B0702                  | AAAAAAAAF               | 0.1                  | 0.50487               | 13.16998             | 1                         | 0.31193 | 0.21892 | 0                                      | 0.07017        | 0                                        | 0.99178271                           | 0.40327581                        | 0.857      | 0.294        | AAAAAAAAAAAAALAAAAAAAAAAAAA                 | 1                              | Ptx     | I547T        | 0.51950689            | 3.7689                | 1.678               |
| rzXB3nQlZ5misn6VN8EA2A== | 0.173                     | BRCA2 | 0.5| AAAAAAAAAAAAAMAAAAAAAAAAAAA | 14                | AAAAAAAAAAAAARAAAAAAAAAAAAA | Ptx               | 0.71575659    | 0.556                     | AAAAAAAAAAAAAMAAAAAAAAAAAAA            | 1        | 1         | 10                   | 90.685                  | 90.685                       | HLA-DQA10401-DQB10402      | HLA-DQA10401-DQB10402         | AAAAAAAAAMAAAAA             | AAAAAAAAARAAAAA                | 421.53                    | 554.92                       | HLA-C*16:01                    | HLA-C*16:01                       | 1                                      | AAAAAAAAM                       | AAAAAAAAR                          | 9                                         | 24.1                          | 6346.9                           | HLA-C*16:01               | HLA-C*16:01                  | AAAAAAAAM                  | AAAAAAAAR                     | 24.1                     | 6346.9                      | 0.25                  | 2.5                      | HLA-DQA10401-DQB10302        | HLA-DQA10401-DQB10302           | AAAAAAAAAAMAAAA               | AAAAAAAAAARAAAA                  | HLA-C*16:01                | HLA-C*16:01                   | AAAAAAAAM                   | AAAAAAAAR                      | 0.0587                    | 8.9317                       | 0.0587               | 8.9317                  | HLA-C*16:01                 | HLA-C*16:01                    | AAAAAAAAM                    | AAAAAAAAR                       | 1        | 1         | 6322.8                        | 1                              | 0.39796                       | 1              | 0.18288                              | 1                    | DPA1_01_03__DPB1_04_01  | AAAAMAAAAAAAAAAA         | 2.44                  | B0702                  | AAAAAAAAM               | 0.07                 | 0.5444                | 39.51379             | 1                         | 0.29303 | 1.5594  | 0                                      | 0.10626        | 0                                        | NA                                   | 0.46452844                        | 0.556      | 0.173        | AAAAAAAAAAAAARAAAAAAAAAAAAA                 | 1                              | Ptx     | E135S        | 0.71575659            | 3.8741                | 1.5417              |

## Tall-skinny format

If the `--with-tall-skinny-table` flag is enabled an output table with the suffix "*_neoantigen_candidates.tsv*" is created. This file contains neoantigen candidates information in short-wide format. Furthermore, a second file with the suffix "*_neoantigen_features.tsv*" is created. This file contains the annotated neoantigen features in tall-skinny format. 

This is a dummy example of a "*_neoantigen_candidates.tsv"*" file with headers following the descriptions in **TABLE 1**:  

| dnaVariantAlleleFrequency | gene  | identifier               | imputedGeneExpression| mutation.mutatedXmer        | mutation.position | mutation.wildTypeXmer       | patientIdentifier | rnaExpression | rnaVariantAlleleFrequency |
|---------------------------|-------|--------------------------|-----------------------------|-----------------------------|-------------------|-----------------------------|-------------------|---------------|---------------------------|
| 0.294                     | BRCA2 | ou11p7RD+tZvjY88DA55Mw== | 0.5 | AAAAAAAAAAAAAFAAAAAAAAAAAAA | 14                | AAAAAAAAAAAAALAAAAAAAAAAAAA | Ptx               | 0.51950689    | 0.857                     |
| 0.173                     | BRCA2 | rzXB3nQlZ5misn6VN8EA2A== | 0.5 | AAAAAAAAAAAAAMAAAAAAAAAAAAA | 14                | AAAAAAAAAAAAARAAAAAAAAAAAAA | Ptx               | 0.71575659    | 0.556                     |

This is the head of a dummy *"_neoantigen_features.tsv"* file:  

| name                            | neoantigen_identifier    | value       |
|---------------------------------|--------------------------|-------------|
| Best_rank_MHCI_score            | ou11p7RD+tZvjY88DA55Mw== | 0.0592      |
| Best_rank_MHCI_score_epitope    | ou11p7RD+tZvjY88DA55Mw== | AAAAAAAAF   |
| Best_rank_MHCI_score_allele     | ou11p7RD+tZvjY88DA55Mw== | HLA-C*16:01 |
| Best_affinity_MHCI_score        | ou11p7RD+tZvjY88DA55Mw== | 24.3        |
| Best_affinity_MHCI_epitope      | ou11p7RD+tZvjY88DA55Mw== | AAAAAAAAF   |
| Best_affinity_MHCI_allele       | ou11p7RD+tZvjY88DA55Mw== | HLA-C*16:01 |
| Best_rank_MHCI_9mer_score       | ou11p7RD+tZvjY88DA55Mw== | 0.0592      |
| Best_rank_MHCI_9mer_epitope     | ou11p7RD+tZvjY88DA55Mw== | AAAAAAAAF   |
| Best_rank_MHCI_9mer_allele      | ou11p7RD+tZvjY88DA55Mw== | HLA-C*16:01 |
| Best_affinity_MHCI_9mer_score   | ou11p7RD+tZvjY88DA55Mw== | 24.3        |
| Best_affinity_MHCI_9mer_allele  | ou11p7RD+tZvjY88DA55Mw== | HLA-C*16:01 |
| Best_affinity_MHCI_9mer_epitope | ou11p7RD+tZvjY88DA55Mw== | AAAAAAAAF   |
| Best_affinity_MHCI_score_WT     | ou11p7RD+tZvjY88DA55Mw== | 21.7        |

\
where:
- name: name of the neoantigen feature, which follow the descriptions given in **TABLE 1**.
- value: value of the noenantigen feaature
- neoantigen_identifier: unique neoantigen id given by NeoFox and the same as `identifier` in the  "*_neoantigen_candidates.tsv*" file.


## JSON format

If the `--with-json` flag is enabled an output file with the suffix "*_neoantigen_candidates.json*" is created. This file contains neoantigen candidates information in JSON format. Furthermore, a second file with the suffix *"_neoantigen_features.json"* is created. This file contains the annotated neoantigen features in JSON format. The names within the models are described in **TABLE 1**.   
\
This is a dummy example of a "*_neoantigen_candidates.json*" file. This file contains a list of neoantigen candidate models (for further information, please see [here](05_models.md). To simplify, only one full neoantigen candidate model is shown:
```json
[{
    "identifier": "ou11p7RD+tZvjY88DA55Mw==",
    "patient_identifier": "Ptx",
    "gene": "BRCA2",
    "mutation": {
        "position": [14],
        "wild_type_xmer": "AAAAAAAAAAAAALAAAAAAAAAAAAA",
        "mutated_xmer": "AAAAAAAAAAAAAFAAAAAAAAAAAAA"
    },
    "rna_expression": 0.5195068939999999,
    "imputed_gene_expression": 0.5,
    "dna_variant_allele_frequency": 0.294,
    "rna_variant_allele_frequency": 0.857
}, {
    "identifier": "rzXB3nQlZ5misn6VN8EA2A==",
    "patient_identifier": "Ptx",
    "gene": "BRCA2",
    "mutation": {
        "position": [14],
        "wild_type_xmer": "AAAAAAAAAAAAARAAAAAAAAAAAAA",
        "mutated_xmer": "AAAAAAAAAAAAAMAAAAAAAAAAAAA"
    },
    "rna_expression": 0.715756594,
    "imputed_gene_expression": 0.5,
    "dna_variant_allele_frequency": 0.17300000000000001,
    "rna_variant_allele_frequency": 0.556
}]
```  

This is a dummy example of a "*_neoantigen_candidates.json*" file. This file contains a list of neoantigen candidate models (for further information, please see [here](05_models.md). To simplify, only one full annotated neoantigen candidate model is shown:

````json
[{
    "neoantigen_identifier": "ou11p7RD+tZvjY88DA55Mw==",
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
        "value": 0.5195068939999999
    }, {
        "name": "+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)",
        "value": "AAAAAAAAAAAAAFAAAAAAAAAAAAA"
    }, {
        "name": "[WT]_+-13_AA_(SNV)_/_-15_AA_to_STOP_(INDEL)",
        "value": "AAAAAAAAAAAAALAAAAAAAAAAAAA"
    }, {
        "name": "VAF_in_tumor",
        "value": 0.294
    }, {
        "name": "VAF_in_RNA",
        "value": 0.857
    }, {
        "name": "Unnamed: 8",
        "value": null
    }],
    "annotator": "Neofox",
    "annotator_version": "0.4.0",
    "timestamp": "20201211115212061465"
}]
````