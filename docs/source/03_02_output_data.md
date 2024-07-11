# Output data

## Neoantigens

NeoFox returns the neoantigen candidates and their annotated features as output. 
Two output formats are supported: tabular format or JSON format.
The user gets the neoantigen annotations in all formats.
Despite different structures, both formats provide the same content with the exception of the metadata on 
the annotations which is only present in the JSON format.

> For each neoantigen candidate (`mutatedXmer`), the best predicted epitope based on MHC prediction is annotated. Columns that refer to a predicted epitope consist of 3 components, each separated by "_". The first part contains the tool name, the second part contains the selection method of the best epitope and the third part contains the name of the value. Example: "ToolName_SelectionMethod_ValueName" -> NetMHCpan_bestRank_peptide (literally: the peptide with the best predicted NetMHCpan rank).  
If no unique selection of the best predicted epitope can be made due to ties, the lexicographically first epitope is returned. If there are nevertheless two or more alleles that have equally good predicted values for the peptide, the lexicographically first one is reported. These cases can be analysed using the --with-all-epitopes options.

The following table describes each of the annotations in the output:  
  
**TABLE 1** 

| Column    name                              | Description                                                                                                                                                                                                                             | Feature group/ Paper              |
|---------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------|
| dnaVariantAlleleFrequency                   | the variant allele frequency calculated from the DNA                                                                                                                                                                                    | -                                 |
| mutatedXmer                                 | the long mutated amino acid sequence                                                                                                                                                                                                    | -                                 |
| wildTypeXmer                                | the long non-mutated amino acid sequence. This field shall be empty for alternative neoantigen classes                                                                                                                                  | -                                 |
| patientIdentifier                           | the patient identifier                                                                                                                                                                                                                  | -                                 |
| rnaExpression                               | the RNA expression. If expression was imputed, this will will be `imputedGeneExpression`                                                                                                                                                | expression                        |
| imputedGeneExpression                       | median gene expression in the TCGA cohort of the tumor entity provided in the patient file.                                                                                                                                             | expression                        |
| rnaVariantAlleleFrequency                   | the variant allele frequency calculated from the RNA                                                                                                                                                                                    | -                                 |
| gene                                        | the HGNC gene symbol                                                                                                                                                                                                                    | -                                 |
| Mutated_rnaExpression_fromRNA               | transcript expression normalized by the variant allele frequency in RNA of the mutation                                                                                                                                                        | expression                        |
| Mutated_rnaExpression_fromDNA               | transcript expression normalized by the variant allele frequency in DNA of the mutation                                                                                                                                                      |                                   |
| Mutated_imputedGeneExpression_fromRNA       | imputeted gene expression normalized by the variant allele frequency in RNA of the mutation                                                                                                                                                        | expression                        |
| Mutated_imputedGeneExpression_fromDNA       | imputed gene expression normalized by the variant allele frequency in RNA of the mutation                                                                                                                                                      |                                   |
| mutation_not_found_in_proteome              | indicates if mutated amino acid sequence was not found in the WT proteome   by exact search                                                                                                                                             | Priority score                    |
| NetMHCpan_bestRank_rank                     | minimal MHC I binding rank score over all neoepitope candidates (8-11mers) and MHC I alleles                                                                                                                                            | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank_peptide                  | neoepitope candidate sequence with minimal MHC I binding rank score                                                                                                                                                                     | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank_core                     | The minimal 9 amino acid binding core directly in contact with the MHC.                                                                                                                                                                 | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank_Icore                    | Interaction core. This is the sequence of the binding core including eventual insertions of deletions.                                                                                                                                  | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank_Of                       | Starting position offset of the core in the peptide (0 based)                                                                                                                                                                           | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank_Gp                       | Position of the deletion (0 based), if any, in the Icore compared to the core.                                                                                                                                                          | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank_Gl                       | Length of the deletion, if any, in the Icore compared to the core.                                                                                                                                                                      | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank_allele                   | the MHC I allele related to `NetMHCpan_bestRank_peptide`                                                                                                                                                                                | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_affinity             | minimal MHC I binding affinity  over all neoepitope   candidates (8-11mers) and MHC I alleles                                                                                                                                           | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_peptide              | neoepitope candidate sequence with minimal MHC I binding affinity                                                                                                                                                                       | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_core                 | The minimal 9 amino acid binding core directly in contact with the MHC.                                                                                                                                                                 | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_Icore                | Interaction core. This is the sequence of the binding core including eventual insertions of deletions.                                                                                                                                  | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_Of                   | Starting position offset of the core in the peptide (0 based)                                                                                                                                                                           | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_Gp                   | Position of the deletion (0 based), if any, in the Icore compared to the core.                                                                                                                                                          | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_Gl                   | Length of the deletion, if any, in the Icore compared to the core.                                                                                                                                                                      | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_allele               | the MHC I allele related to `NetMHCpan_bestAffinity_peptide`                                                                                                                                                                            | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank9mer_rank                 | minimal MHC I binding rank score over all neoepitope candidates (9mers only) and MHC I alleles                                                                                                                                          | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank9mer_peptide              | neoepitope candidate sequence (9mer) with minimal MHC I binding rank score                                                                                                                                                              | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank9mer_allele               | the MHC I allele related to `NetMHCpan_bestRank9mer_peptide`                                                                                                                                                                            | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity9mer_affinity         | minimal MHC I binding affinity over all neoepitope candidates (9mers) and MHC I alleles                                                                                                                                                 | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity9mer_allele           | the MHC I allele related to `NetMHCpan_bestAffinity9mer_peptide             `                                                                                                                                                           | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity9mer_peptide          | neoepitope candidate sequence (9mer) with minimal MHC I binding affinity                                                                                                                                                                | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_affinityWT           | MHC I binding affinity of `NetMHCpan_bestAffinity_peptideWT`                                                                                                                                                                            | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity_peptideWT            | WT epitope that corresponds to `NetMHCpan_bestAffinity_peptide`                                                                                                                                                                         | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank_rankWT                   | MHC I binding rank score of `NetMHCpan_bestRank_peptideWT`                                                                                                                                                                              | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank_peptideWT                | WT epitope that corresponds to `NetMHCpan_bestRank_peptide`                                                                                                                                                                             | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank9mer_rankWT               | MHC I binding rank score of `NetMHCpan_bestRank9mer_peptideWT `                                                                                                                                                                         | MHC I binding with netMHCpan      |
| NetMHCpan_bestRank9mer_peptideWT            | WT epitope that corresponds to `NetMHCpan_bestRank9mer_peptide`                                                                                                                                                                         | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity9mer_affinityWT       | MHC I binding affinity of ` NetMHCpan_bestAffinity9mer_peptideWT `                                                                                                                                                                      | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity9mer_peptideWT        | WT epitope that corresponds to `NetMHCpan_bestAffinity9mer_peptide`                                                                                                                                                                     |                                   |
| PHBR_I                                      | harmonic mean of minimal MHC I binding rank scores of all MHC I alleles of a patient                                                                                                                                                    | PHBR-I                            |
| NetMHCpan_bestAffinity9mer_positionMutation | indicates position of the mutation in ` NetMHCpan_bestRank9mer_peptide`                                                                                                                                                                 | MHC I binding with netMHCpan      |
| NetMHCpan_bestAffinity9mer_anchorMutated    | mutation in `NetMHCpan_bestRank9mer_peptide` in an anchor position (i.e. position 2 or   9)                                                                                                                                             | anchor/non-anchor                 |
| NetMHCIIpan_bestRank_rank                   | minimal MHC II binding rank score over all neoepitope candidates (15mers) and all MHC II alleles                                                                                                                                        | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestRank_peptide                | neoepitope candidate sequence with minimal MHC II binding rank score                                                                                                                                                                    | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestRank_core                   | Binding core register (9mer)                                                                                                                                                                                                            | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestRank_Of                     | Starting position offset of the optimal binding core (starting from 0)                                                                                                                                                                  | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestRank_coreRel                | Reliability of the binding core, expressed as the fraction of networks in the ensemble selecting the optimal core                                                                                                                       | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestRank_allele                 | the MHC II isoform related to `NetMHCIIpan_bestRank_peptide`                                                                                                                                                                            | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestAffinity_affinity           | minimal MHC II binding affinity over all neoepitope candidates (15mers) and all MHC II alleles                                                                                                                                          | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestAffinity_peptide            | neoepitope candidate sequence with minimal MHC II binding affinity                                                                                                                                                                      | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestAffinity_core               | Binding core register (9mer)                                                                                                                                                                                                            | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestAffinity_Of                 | Starting position offset of the optimal binding core (starting from 0)                                                                                                                                                                  | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestAffinity_coreRel            | Reliability of the binding core, expressed as the fraction of networks in the ensemble selecting the optimal core                                                                                                                       | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestAffinity_allele             | the MHC II isoform related to `NetMHCIIpan_bestAffinity_peptide`                                                                                                                                                                        | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestRank_rankWT                 | minimal MHC II binding rank  of `NetMHCIIpan_bestRank_peptideWT`                                                                                                                                                                        | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestRank_peptideWT              | WT epitope sequence (15mer) that corresponds to ` NetMHCIIpan_bestRank_peptide`                                                                                                                                                         | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestAffinity_affinityWT         | minimal MHC II binding rank of `NetMHCIIpan_bestAffinity_peptideWT`                                                                                                                                                                     | MHC II binding with   netMHCIIpan |
| NetMHCIIpan_bestAffinity_peptideWT          | WT epitope sequence (15mer) that corresponds to `NetMHCIIpan_bestAffinity_peptide`                                                                                                                                                      | MHC II binding with   netMHCIIpan |
| PHBR_II                                     | harmonic mean of minimal MHC II binding rank scores of all MHC II alleles of a patient                                                                                                                                                  | PHBR-II                           |
| Amplitude_MHCI_bestAffinity9mer             | ratio of `NetMHCpan_bestAffinity9mer_affinity` and `NetMHCpan_bestAffinity9mer_affinityWT`                                                                                                                                              | Recognition Potential             |
| Amplitude_MHCI_bestAffinity                 | ratio of `NetMHCpan_bestAffinity_affinityWT` and `NetMHCpan_bestAffinity_affinity`                                                                                                                                                      | Generator rate                    |
| Amplitude_MHCII_bestRank                    | ratio of `NetMHCIIpan_bestRank_rank` and `NetMHCIIpan_bestRank_rankWT`                                                                                                                                                                  | Generator rate                    |
| Pathogensimiliarity_MHCI_bestAffinity9mer   | score representing the similarity of `NetMHCpan_bestAffinity9mer_peptide` to pathogen sequences in IEDB database                                                                                                                        | Recognition Potential             |
| Pathogensimiliarity_MHCII_bestAffinity      | score representing the similarity of `NetMHCIIpan_bestRank_peptide` to pathogen sequences in IEDB database                                                                                                                              | Recognition Potential             |
| RecognitionPotential_MHCI_bestAffinity9mer  | product of `Amplitude_MHCI_affinity_9mer` and `Pathogensimiliarity_MHCI_affinity_9mer`                                                                                                                                                  | Recognition Potential             |
| DAI_MHCI_bestAffinity                       | difference of `NetMHCpan_bestAffinity_affinityWT` and `NetMHCpan_bestAffinity_affinity`                                                                                                                                                 | DAI                               |
| Classically_defined_neopeptide_MHCI         | `NetMHCpan_bestAffinity_peptide`< 50 nM                                                                                                                                                                                                 | Generator rate                    |
| Alternatively_defined_neopeptide_MHCI       | `NetMHCpan_bestAffinity_peptide` < 5000 nM and `Amplitude_MHCI_bestAffinity` > 10                                                                                                                                                       | Generator rate                    |
| Classically_defined_neopeptide_MHCII        | `NetMHCIIpan_bestRank_rank` < 1                                                                                                                                                                                                         | Generator rate                    |
| Alternatively_defined_neopeptide_MHCII      | `Best_rank_MHCII_score` < 4 and `Amplitude_MHCII_bestRank` < 2                                                                                                                                                                          | Generator rate                    |
| GeneratorRate_CDN_MHCI                      | number of neoepitope candidates with MHC I binding   affinity < 50 nM per neoantigen canidate                                                                                                                                           | Generator rate                    |
| GeneratorRate_ADN_MHCI                      | number of neoepitope candidates with MHC I binding affinity   < 5000 nM per neoantigen canidate 10x better affinity in comparison to corresponding WT peptide                                                                           | Generator rate                    |
| GeneratorRate_MHCI                          | sum of `GeneratorRate_CDN_MHCI` and `GeneratorRate_ADN_MHCI`                                                                                                                                                                            | Generator rate                    |
| GeneratorRate_CDN_MHCII                     | number of neoepitope candidates with MHC II binding rank score < 1 per neoantigen canidate                                                                                                                                              | Generator rate                    |
| GeneratorRate_ADN_MHCII                     | number of neoepitope candidates with MHC II binding rank score < 4 per neoantigen candidate 4x better rank in comparison to corresponding WT peptide                                                                                    | Generator rate                    |
| GeneratorRate_MHCII                         | sum of `GeneratorRate_CDN_MHCII` and `GeneratorRate_ADN_MHCII`                                                                                                                                                                          | Generator rate                    |
| ImprovedBinder_MHCI                         | ratio of `NetMHCpan_MHCI_rank_bestRankWT` and `NetMHCpan_MHCI_rank_bestRank` > 1.2                                                                                                                                                      | self-similarity                   |
| Selfsimilarity_MHCI_conserved_binder        | score representing the similarity between `NetMHCpan_bestRank_peptide` and `NetMHCpan_bestRank_peptideWT` For conservered binder only                                                                                                   | self-similarity                   |
| Selfsimilarity_MHCI                         | score representing the similarity between `NetMHCpan_bestRank_peptide` and `NetMHCpan_bestRank_peptide`                                                                                                                                 | self-similarity                   |
| Selfsimilarity_MHCII                        | score representing the similarity between `NetMHCIIpan_bestAffinity_peptide` and `NetMHCIIpan_bestAffinity_peptide`                                                                                                                     | self-similarity                   |
| Number_of_mismatches_MCHI                   | number of amino acids that do no match between `NetMHCpan_bestRank_peptide` and `NetMHCpan_bestRank_peptideWT`                                                                                                                          | Priority score                    |
| Priority_score_fromDNA                      | combinatorial score of several features such as MHC binding, transcription expression and VAF in DNA                                                                                                                                          | Priority score                    |
| Priority_score_fromRNA                      | combinatorial score of several features such as MHC binding, transcription expression and VAF in RNA                                                                                                                                          | Priority score                    |
| Priority_score_imputed_fromDNA              | combinatorial score of several features such as MHC binding, imputed gene expression and VAF  in DNA                                                                                                                                          | Priority score                    |
| Priority_score_imputed_fromRNA              | combinatorial score of several features such as MHC binding, imputed gene expression and VAF  in RNA                                                                                                                                         | Priority score                    |
| IEDB_Immunogenicity_MHCI                    | IEDB Immunogenicity score for `NetMHCpan_bestAffinity_peptide`                                                                                                                                                                          | IEDB Immunogenicity               |
| IEDB_Immunogenicity_MHCII                   | IEDB Immunogenicity score for `NetMHCIIpan_bestAffinity_peptide`                                                                                                                                                                        | IEDB Immunogenicity               |
| MixMHCpred_bestScore_peptide                | MHC class I neoepitope candidate sequence with maximum MixMHCpred score over all neoepitope canidates (8-11mers) and MHC I alleles                                                                                                      | MixMHCpred                        |
| MixMHCpred_bestScore_score                  | maximum MixMHCpred score over all neoepitope canidates (8-11mers) and MHC I alleles                                                                                                                                                     | MixMHCpred                        |
| MixMHCpred_bestScore_rank                   | rank that corresponds to `MixMHCpred_bestScore_score`                                                                                                                                                                                   | MixMHCpred                        |
| MixMHCpred_bestScore_allele                 | the allele with maximum MixMHCpred score                                                                                                                                                                                                | MixMHCpred                        |
| MixMHC2pred_bestRank_peptide                | MHC class II neoepitope candidate sequence with minimal MixMHC2pred score over all neoepitope canidates (13-18mers) and MHC II alleles                                                                                                  | MixMHC2pred                       |
| MixMHC2pred_bestRank_rank                   | minimal MixMHC2pred score over all neoepitope canidates (13-18mers) and MHC II alleles                                                                                                                                                  | MixMHC2pred                       |
| MixMHC2pred_bestRank_allele                 | the MHC II isoform with minimum MixMHC2pred rank score                                                                                                                                                                                  | MixMHC2pred                       |
| Dissimilarity_MHCI                          | score reflecting the dissimilarity of `NetMHCpan_bestAffinity_peptide` to the self-proteome                                                                                                                                             | dissimilarity                     |
| Dissimilarity_MHCII                         | score reflecting the dissimilarity of `NetMHCIIpan_bestAffinity_peptide` to the self-proteome                                                                                                                                           | dissimilarity                     |
| Vaxrank_bindingScore                        | total binding score of vaxrank                                                                                                                                                                                                          | vaxrank                           |
| Vaxrank_totalScore                          | product of total binding score and transcription expression score. Originally, the root of the number of reads supporting the mutation are used in the original implementation. To simplify, the transcript expression normalised  to VAF is used. | vaxrank                           |
| Vaxrank_totalScore_imputed                  | product of total binding score and imputed gene expression score. Originally, the root of the number of reads supporting the mutation are used in the original implementation. To simplify, the imputed gene expression normalised  to VAF is used.  | vaxrank                           |
| PRIME_bestScore_allele                      | best predicted MHC allele by PRIME model                                                                                                                                                                                                | PRIME                             |
| PRIME_bestScore_peptide                     | best predicted neoepitope candidate by PRIME model                                                                                                                                                                                      | PRIME                             |
| PRIME_bestScore_rank                        | output rank score of PRIME model                                                                                                                                                                                                        | PRIME                             |
| PRIME_bestScore_score                       | output score of PRIME model                                                                                                                                                                                                             | PRIME                             |
| HexAlignmentScore_MHCI                      | the alignment score by HEX for `NetMHCpan_bestAffinity_peptide`                                                                                                                                                                         | HEX                               |
| HexAlignmentScore_MHCII                     | the alignment score by HEX for `NetMHCIIpan_bestAffinity_peptide`                                                                                                                                                                       | HEX                               |


In addition, all logging output is appended to a log file with the suffix
"*<folder>/<prefix>.log*", where the folder is set by `--output-folder` and the
prefix can be set with `--output-prefix`.



### Tabular format

An output table with the suffix "*_neoantigen_candidates_annotated.tsv*" is created.  
This table contains the neoantigen candidates information, the neoantigen annotations and external annotations if they were present provided in the input table.

This is a dummy example:  

| patientIdentifier | gene  | mutatedXmer                 | wildTypeXmer                | position | dnaVariantAlleleFrequency | rnaVariantAlleleFrequency | rnaExpression | imputedGeneExpression | Alternatively_defined_neopeptide_MHCI | Alternatively_defined_neopeptide_MHCII | Amplitude_MHCII_bestRank | Amplitude_MHCI_bestAffinity | Amplitude_MHCI_bestAffinity9mer | Classically_defined_neopeptide_MHCI | Classically_defined_neopeptide_MHCII | DAI_MHCI_bestAffinity | Dissimilarity_MHCI | Dissimilarity_MHCII | GeneratorRate_ADN_MHCI | GeneratorRate_ADN_MHCII | GeneratorRate_CDN_MHCI | GeneratorRate_CDN_MHCII | GeneratorRate_MHCI | GeneratorRate_MHCII | HexAlignmentScore_MHCI | HexAlignmentScore_MHCII | IEDB_Immunogenicity_MHCI | IEDB_Immunogenicity_MHCII | Improved_Binder_MHCI | MixMHC2pred_bestRank_allele | MixMHC2pred_bestRank_peptide | MixMHC2pred_bestRank_rank | MixMHCpred_bestScore_allele | MixMHCpred_bestScore_peptide | MixMHCpred_bestScore_rank | MixMHCpred_bestScore_score | Mutated_imputedGeneExpression_fromDNA | Mutated_imputedGeneExpression_fromRNA | Mutated_rnaExpression_fromDNA | Mutated_rnaExpression_fromRNA | NetMHCIIpan_bestAffinity_Of | NetMHCIIpan_bestAffinity_affinity | NetMHCIIpan_bestAffinity_affinityWT | NetMHCIIpan_bestAffinity_allele | NetMHCIIpan_bestAffinity_core | NetMHCIIpan_bestAffinity_coreRel | NetMHCIIpan_bestAffinity_peptide | NetMHCIIpan_bestAffinity_peptideWT | NetMHCIIpan_bestRank_Of | NetMHCIIpan_bestRank_allele | NetMHCIIpan_bestRank_core | NetMHCIIpan_bestRank_coreRel | NetMHCIIpan_bestRank_peptide | NetMHCIIpan_bestRank_peptideWT | NetMHCIIpan_bestRank_rank | NetMHCIIpan_bestRank_rankWT | NetMHCpan_bestAffinity9mer_affinity | NetMHCpan_bestAffinity9mer_affinityWT | NetMHCpan_bestAffinity9mer_allele | NetMHCpan_bestAffinity9mer_anchorMutated | NetMHCpan_bestAffinity9mer_peptide | NetMHCpan_bestAffinity9mer_peptideWT | NetMHCpan_bestAffinity9mer_positionMutation | NetMHCpan_bestAffinity_Gl | NetMHCpan_bestAffinity_Gp | NetMHCpan_bestAffinity_Icore | NetMHCpan_bestAffinity_Of | NetMHCpan_bestAffinity_affinity | NetMHCpan_bestAffinity_affinityWT | NetMHCpan_bestAffinity_allele | NetMHCpan_bestAffinity_core | NetMHCpan_bestAffinity_peptide | NetMHCpan_bestAffinity_peptideWT | NetMHCpan_bestRank9mer_allele | NetMHCpan_bestRank9mer_peptide | NetMHCpan_bestRank9mer_peptideWT | NetMHCpan_bestRank9mer_rank | NetMHCpan_bestRank9mer_rankWT | NetMHCpan_bestRank_Gl | NetMHCpan_bestRank_Gp | NetMHCpan_bestRank_Icore | NetMHCpan_bestRank_Of | NetMHCpan_bestRank_allele | NetMHCpan_bestRank_core | NetMHCpan_bestRank_peptide | NetMHCpan_bestRank_peptideWT | NetMHCpan_bestRank_rank | NetMHCpan_bestRank_rankWT | Number_of_mismatches_MCHI | PHBR_I | PHBR_II | PRIME_best_allele | PRIME_best_peptide | PRIME_best_rank | PRIME_best_score | Pathogensimiliarity_MHCII_bestAffinity | Pathogensimiliarity_MHCI_bestAffinity | Pathogensimiliarity_MHCI_bestAffinity9mer | Priority_score_fromDNA | Priority_score_fromRNA | Priority_score_imputed_fromDNA | Priority_score_imputed_fromRNA | RecognitionPotential_MHCI_bestAffinity | RecognitionPotential_MHCI_bestAffinity9mer | Selfsimilarity_MHCI | Selfsimilarity_MHCII | Selfsimilarity_MHCI_conserved_binder | Vaxrank_bindingScore | Vaxrank_totalScore | Vaxrank_totalScore_imputed | mutation_not_found_in_proteome |
|-------------------|-------|-----------------------------|-----------------------------|----------|---------------------------|---------------------------|---------------|-----------------------|---------------------------------------|----------------------------------------|--------------------------|-----------------------------|---------------------------------|-------------------------------------|--------------------------------------|-----------------------|--------------------|---------------------|------------------------|-------------------------|------------------------|-------------------------|--------------------|---------------------|------------------------|-------------------------|--------------------------|---------------------------|----------------------|-----------------------------|------------------------------|---------------------------|-----------------------------|------------------------------|---------------------------|----------------------------|---------------------------------------|---------------------------------------|-------------------------------|-------------------------------|-----------------------------|-----------------------------------|-------------------------------------|---------------------------------|-------------------------------|----------------------------------|----------------------------------|------------------------------------|-------------------------|-----------------------------|---------------------------|------------------------------|------------------------------|--------------------------------|---------------------------|-----------------------------|-------------------------------------|---------------------------------------|-----------------------------------|------------------------------------------|------------------------------------|--------------------------------------|---------------------------------------------|---------------------------|---------------------------|------------------------------|---------------------------|---------------------------------|-----------------------------------|-------------------------------|-----------------------------|--------------------------------|----------------------------------|-------------------------------|--------------------------------|----------------------------------|-----------------------------|-------------------------------|-----------------------|-----------------------|--------------------------|-----------------------|---------------------------|-------------------------|----------------------------|------------------------------|-------------------------|---------------------------|---------------------------|--------|---------|-------------------|--------------------|-----------------|------------------|----------------------------------------|---------------------------------------|-------------------------------------------|------------------------|------------------------|--------------------------------|--------------------------------|----------------------------------------|--------------------------------------------|---------------------|----------------------|--------------------------------------|----------------------|--------------------|----------------------------|--------------------------------|
| Ptx               | BRCA2 | AAAAAAAAAAAAALAAAAAAAAAAAAA | AAAAAAAAAAAAAFAAAAAAAAAAAAA | 14       | 0.857                     | 0.294                     | 0.51950689    | 0.53659963            | 0                                     | 0                                      | 0.33333                  | 2.3262                      | 0.56811                         | 0                                   | 1                                    | 7165.9                | 1                  | 3.00E-05            | 0                      | 1                       | 0                      | 18                      | 0                  | 19                  | 177                    | 389                     | 0.18182                  | 0.36258                   | 0                    | HLA-DQA1*01:02-DQB1*06:02   | AAAAAAALAAAAAAA              | 0.00438                   | HLA-B*07:02                 | AAALAAAAA                    | 0.18911                   | -0.2706                    | 0.45987                               | 0.15776                               | 0.44522                       | 0.15274                       | 3                           | 14.84                             | 21.8                                | HLA-DQA1*01:02-DQB1*06:02       | AAAAALAAA                     | 0.74                             | AAAAAAAALAAAAAA                  | AAAAAAAAFAAAAAA                    | 3                       | HLA-DQA1*03:01-DQB1*06:02   | LAAAAAAAA                 | 0.75                         | AAALAAAAAAAAAAA              | AAAFAAAAAAAAAAA                | 0.03                      | 0.01                        | 1534.8                              | 1180.8                                | HLA-C*03:04                       | 0                                        | AAALAAAAA                          | AAAFAAAAA                            | 4                                           | 1                         | 1                         | AALAAAAAAA                   | 0                         | 1018.2                          | 8184.1                            | HLA-A*02:01                   | ALAAAAAAA                   | AALAAAAAAA                     | AAFAAAAAAA                       | HLA-B*07:02                   | AAALAAAAA                      | AAAFAAAAA                        | 4.389                       | 5.068                         | 0                     | 0                     | AAALAAAAA                | 0                     | HLA-B*07:02               | AAALAAAAA               | AAALAAAAA                  | AAAFAAAAA                    | 4.389                   | 5.068                     | 1                         | 8.5743 | NA      | HLA-B*07:02       | AAALAAAAA          | 0.114           | 0.10131          | 0                                      | 0                                     | 0                                         | 0                      | 0                      | 0                              | 0                              | 0                                      | 0                                          | 0.97634652          | 0.98126377           | 0.97634652                           | 0.01301              | 0.00199            | 0.00205                    | 1                              |
| Ptx               | BRCA2 | AAAAAAAAAAAAARAAAAAAAAAAAAA | AAAAAAAAAAAAAMAAAAAAAAAAAAA | 14       | 0.556                     | 0.173                     | 0.71575659    | 0.53659963            | 0                                     | 1                                      | 390                      | 5.6409                      | 5.6409                          | 0                                   | 1                                    | 1498                  | 1                  | 0                   | 0                      | 18                      | 0                      | 27                      | 0                  | 45                  | 132                    | 409                     | 0.19477                  | 0.42378                   | 1                    | HLA-DQA1*01:02-DQB1*06:02   | AAAARAAAAAAAAAA              | 0.00438                   | HLA-B*07:02                 | AAAAAARAA                    | 0.08828                   | 0.04913                    | 0.29835                               | 0.09283                               | 0.39796                       | 0.12383                       | 3                           | 13.55                             | 15.81                               | HLA-DQA1*01:02-DQB1*06:02       | AAAARAAAA                     | 0.97                             | AAAAAAARAAAAAAA                  | AAAAAAAMAAAAAAA                    | 3                       | HLA-DPA1*01:03-DPB1*06:01   | RAAAAAAAA                 | 1                            | AAARAAAAAAAAAAA              | AAAMAAAAAAAAAAA                | 0.01                      | 3.9                         | 199.38                              | 1697.4                                | HLA-B*07:02                       | 0                                        | AAAAARAAA                          | AAAAAMAAA                            | 6                                           | 0                         | 0                         | AAAAARAAA                    | 0                         | 199.38                          | 1697.4                            | HLA-B*07:02                   | AAAAARAAA                   | AAAAARAAA                      | AAAAAMAAA                        | HLA-B*07:02                   | AAAAARAAA                      | AAAAAMAAA                        | 1.346                       | 4.347                         | 0                     | 0                     | AAAAARAAA                | 0                     | HLA-B*07:02               | AAAAARAAA               | AAAAARAAA                  | AAAAAMAAA                    | 1.346                   | 4.347                     | 1                         | 4.6189 | NA      | HLA-B*07:02       | AAAAAARAA          | 0.129           | 0.09563          | 0                                      | 0                                     | 0                                         | 0.32903                | 0.10238                | 0.26268                        | 0.08173                        | 0                                      | 0                                          | 0.9663076           | 0.97299187           | NA                                   | 0.80336              | 0.09948            | 0.07458                    | 1                              |

### JSON format

An output file with the suffix "*_neoantigen_candidates_annotated.json*" is created.  
This file contains neoantigen candidates information in JSON format.  
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
        "name": "NetMHCpan_MHCI_rank_bestRank",
        "value": "0.0592"
    }, {
        "name": "NetMHCpan_MHCI_rank_peptide",
        "value": "AAAAAAAAF"
    }],
    "annotator": "Neofox",
    "annotator_version": "1.0.0",
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

The following table describes each of the annotations in the output. MHC-I or MHC-II specific features will be only available in the respective table:

**TABLE 2**

| Column   Name                   | Description                                                                                                                        | Feature group/ Paper                                    |
|---------------------------------|------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------|
| dnaVariantAlleleFrequency       | the variant allele frequency calculated from the DNA                                                                               | -                                                       |
| mutatedSequence                 | the mutated amino acid sequence                                                                                                    | -                                                       |
| wildTypeSequence                | the non-mutated amino acid sequence (when not provided in the input this will contain the Blastp closest sequence in the proteome) | -                                                       |
| core             |              | MHC core part of the peptide ligand that primarily interacts with the MHC binding groove                                           | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan
| alleleMhcI / isoformMhcII       | either the MHC-I allele for MHC-I neoepitopes or the MHC-II isoform for MHC-II neoepitopes                                         | -                                                       |
| patientIdentifier               | the patient identifier (optional)                                                                                                  | -                                                       |
| rnaExpression                   | the RNA expression. If expression was imputed, this will will be `imputedGeneExpression`                                           | expression                                              |
| imputedGeneExpression           | median gene expression in the TCGA cohort of the tumor entity provided in the patient file.                                        | expression                                              |
| rnaVariantAlleleFrequency       | the variant allele frequency calculated from the RNA                                                                               | -                                                       |
| gene                            | the HGNC gene symbol                                                                                                               | -                                                       |
| Mutated_rnaExpression_fromRNA               | transcript expression normalized by the variant allele frequency in RNA of the mutation                                                                                                                                                        | expression                        |
| Mutated_rnaExpression_fromDNA               | transcript expression normalized by the variant allele frequency in DNA of the mutation                                                                                                                                                      |                                   |
| Mutated_imputedGeneExpression_fromRNA       | imputeted gene expression normalized by the variant allele frequency in RNA of the mutation                                                                                                                                                        | expression                        |
| Mutated_imputedGeneExpression_fromDNA       | imputed gene expression normalized by the variant allele frequency in RNA of the mutation                                                                                                                                                      |                                   |
| affinityMutated                 | NetMHCpan / NetMHCIIpan affinity score of the mutated peptide for MHC-I or MHC-II neoepitopes respectively                         | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan     |
| affinityWildType                | NetMHCpan / NetMHCIIpan affinity score of the wild type peptide for MHC-I or MHC-II neoepitopes respectively                       | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan     |
| rankMutated                     | NetMHCpan / NetMHCIIpan rank of the mutated peptide for MHC-I or MHC-II neoepitopes respectively                                   | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan     |
| rankWildType                    | NetMHCpan / NetMHCIIpan rank of the wild type peptide for MHC-I or MHC-II neoepitopes respectively                                 | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan     |
| MixMHCpred_score                | MixMHCpred score of the mutated peptide for MHC-I neoepitopes                                                                      | MHC-I binding with mixMHCpred                           |
| MixMHCpred_rank                 | MixMHCpred rank of the mutated peptide for MHC-I neoepitopes                                                                       | MHC-I binding with mixMHCpred                           |
| MixMHCpred_WT_score             | MixMHCpred score of the wild type peptide for MHC-I neoepitopes                                                                    | MHC-I binding with mixMHCpred                           |
| MixMHCpred_WT_rank              | MixMHCpred rank of the wild type peptide for MHC-I neoepitopes                                                                     | MHC-I binding with mixMHCpred                           |
| MixMHC2pred_score               | MixMHC2pred score of the mutated peptide for MHC-II neoepitopes                                                                    | MHC-II binding with mixMHC2pred                         |
| MixMHC2pred_rank                | MixMHC2pred rank of the mutated peptide for MHC-II neoepitopes                                                                     | MHC-II binding with mixMHC2pred                         |
| MixMHC2pred_WT_score            | MixMHC2pred score of the wild type peptide for MHC-II neoepitopes                                                                  | MHC-II binding with mixMHC2pred                         |
| MixMHC2pred_WT_rank             | MixMHC2pred rank of the wild type peptide for MHC-II neoepitopes                                                                   | MHC-II binding with mixMHC2pred                         |
| PRIME_score                     | PRIME score of the mutated peptide for MHC-I neoepitopes                                                                           | MHC-I binding with PRIME                                |
| PRIME_rank                      | PRIME rank of the mutated peptide for MHC-I neoepitopes                                                                            | MHC-I binding with PRIME                                |
| PRIME_WT_score                  | PRIME score of the wild type peptide for MHC-I neoepitopes                                                                         | MHC-I binding with PRIME                                |
| PRIME_WT_rank                   | PRIME rank of the wild type peptide for MHC-I neoepitopes                                                                          | MHC-I binding with PRIME                                |
| DAI                             | difference of `affinityWildType` and `affinityMutated`                                                                             | DAI (only availble for MHC-I)                           |
| Gl                              | Length of the deletion (in the core), if any.                                                                                      | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan
| Gp                              | Position of the deletion (in the core), if any.                                                                                    | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan
| Icore                           | Interaction core. This is the sequence of the binding core including eventual insertions of deletions.                             | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan
| Of                              | The starting position of the Core within the predicted peptide                                                                     | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan
| Core_Rel                        | Reliability of the (MHCII) binding core, expressed as the fraction of networks in the ensemble selecting the optimal core          | MHC-I / MHC-II binding with netMHCpan / netMHCIIpan 
| IEDB_Immunogenicity             | IEDB Immunogenicity score for `affinityMutated`                                                                                    | IEDB immunogenicity                                     |
| Improved_Binder_MHCI            | ratio of `affinityWildType` and `affinityMutated` > 1.2                                                                            | self-similarity (only available for MHC-I)              |
| Priority_score_fromDNA                      | combinatorial score of several features such as MHC binding, transcription expression and VAF in DNA                                                                                                                                          | Priority score                    |
| Priority_score_fromRNA                      | combinatorial score of several features such as MHC binding, transcription expression and VAF in RNA                                                                                                                                          | Priority score                    |
| Priority_score_imputed_fromDNA              | combinatorial score of several features such as MHC binding, imputed gene expression and VAF  in DNA                                                                                                                                          | Priority score                    |
| Priority_score_imputed_fromRNA              | combinatorial score of several features such as MHC binding, imputed gene expression and VAF  in RNA                                                                                                                                         | Priority score                    |
| mutation_not_found_in_proteome  | indicates if mutated amino acid sequence was not found in the WT proteome by exact search                                          | Priority score                                          |
| Selfsimilarity                  | score representing the similarity between `rankMutated` and `rankWildType`                                                         | self-similarity                                         |
| Selfsimilarity_conserved_binder | score representing the similarity between `rankMutated` and `rankWildType` for conserved binder only                               | self-similarity (only available for MHC-I)              |
| dissimilarity_score             | score reflecting the dissimilarity of `affinityMutated` to the self-proteome                                                       | dissimilarity                                           |
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
