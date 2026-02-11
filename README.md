
Last Updated on 2026-02-11

-   [Supplementary Materials](#supplementary-materials)
    -   [Supplementary Data S1: NSCLC Patient Information](#supplementary-data-s1-nsclc-patient-information)
    -   [Supplementary Data S2: NSCLC Protein-Affecting Variants](#supplementary-data-s2-nsclc-protein-affecting-variants)
    -   [Supplementary Data S3: NSCLC Missense Variants](#supplementary-data-s3-nsclc-missense-variants)
    -   [Supplementary Data S4: NSCLC HLA Loss of Heterozygosity](#supplementary-data-s4-nsclc-hla-loss-of-heterozygosity)
    -   [Supplementary Data S5: NSCLC Shared Protein Lists](#supplementary-data-s5-nsclc-shared-protein-lists)
    -   [Supplementary Data S6: NSCLC pVACseq Class I Predictions](#supplementary-data-s6-nsclc-pvacseq-class-i-predictions)
    -   [Supplementary Data S7: NSCLC pVACseq Class II Predictions](#supplementary-data-s7-nsclc-pvacseq-class-ii-predictions)
    -   [Supplementary Data S8: NSCLC pVACseq Peptidome Combined Predictions](#supplementary-data-s8-nsclc-pvacseq-peptidome-combined-predictions)
    -   [Supplementary Data S9: NSCLC Tested Neoantigens](#supplementary-data-s9-nsclc-tested-neoantigens)
    -   [References](#references)

# Supplementary Materials

This repository contains the supplementary data files for:
**Immunopeptidomics-guided identification of functional neoantigens in non-small cell lung cancer** \[[1](#ref-nicholas2024)\]. This repository
is associated with
[![DOI](https://zenodo.org/badge/833647778.svg)](https://zenodo.org/doi/10.5281/zenodo.12820423)

**For Supplementary Figures S1-S11 and Supplementary Tables S1-S2, please refer to [supplement-2026-02-11.pdf](supplement-2026-02-11.pdf).**

Supplementary Data S1 to S9 are `csv` files. The column names and contents of these files are described below.

## Supplementary Data S1: NSCLC Patient Information

Supplementary Data S1, a `csv` file containing patient information with 24 rows and 21 column variables. Each row in Supplementary Data S1 represents observations for a single patient.

| Column Variable | Description |
|-----------------------------------------|-------------------------------|
| accel_id | CRUK Accelerator patient identifier |
| target_lung_id | Targeted Lung Health Check patient identifier |
| tissue | NSCLC type: Adenocarcinoma or Squamous cell carcinoma |
| n_somatic_variants | Total number of somatic variants identified by whole exome sequencing |
| mut_burden_per_mb | Mutational burden: mutations per million bases of DNA. Exome target size was 35.7 Mb |
| obs_class_I | Number of observed HLA I peptides by mass spec. immunopeptidomics |
| obs_class_II | Number of observed HLA II peptides by mass spec. immunopeptidomics |
| HLA | Class I and II HLA allotypes identified by genomic sequencing |
| wet_weight | Wet weight of tumour tissue |
| tumour_purity | Tumour purity as calculated from WES by ASCAT |
| tumour_ploidy | Tumour ploidy as calculated from WES by ASCAT |
| til_status | Tumour infiltrating T-cell status by immunohistochemistry: Low, Moderate, High or NA |
| weeks_post_surgery | Number of weeks since surgery |
| status_as_of_2021_01_19 | Status since 2021-01-19: Alive, Deceased or NA |
| date_of_diagnosis | Date of diagnosis |
| smoking_status | Smoking status |
| notes_2 | Notes about smoking history |

## Supplementary Data S2: NSCLC Protein-Affecting Variants

| Column Variable | Description |
|-----------------------------------------|-------------------------------|
| accel_id | CRUK Accelerator patient identifier |
| vid | Unique variant identifier |
| chrom | Chromosome |
| pos | Genomic coordinate |
| ref | Reference base |
| alt | Variant base |
| type | Variant type: `snv`, `ins`, `del` or `complex`. Single nucleotide variant, insertion, deletion and complex variant respectively |
| gene_name | HGNC gene name |
| ensg | Ensembl gene identifier |
| ensp | Ensembl protein identifier |
| consequence | Variant consequence annotation |
| biotype | Gene biotype filter: protein_coding |
| canonical | Canonical transcript flag |
| vaf | Variant allele frequency |
| depth | Sequencing depth |
| filter | Variant filter status |
| info | Information field from VCF file |
| format | Format of VCF variable columns |
| sample_1 | Reference sample VCF variable values corresponding with format |
| sample_2 | Tumour sample VCF variable values corresponding with format |
| tissue | Lung tumour tissue type: `Squamous` or `Adenocarcinoma` |

## Supplementary Data S3: NSCLC Missense Variants

| Column Variable | Description |
|-----------------------------------------|-------------------------------|
| accel_id | CRUK Accelerator patient identifier |
| vid | Unique variant identifier |
| chrom | Chromosome |
| pos | Genomic coordinate |
| ref | Reference base |
| alt | Variant base |
| type | Variant type: `snv`, `ins`, `del` or `complex`. Single nucleotide variant, insertion, deletion and complex variant respectively |
| gene_name | HGNC gene name |
| ensg | Ensembl gene identifier |
| ensp | Ensembl protein identifier |
| vaf | Variant allele frequency |
| depth | Sequencing depth |
| filter | Variant filter status |
| info | Information field from VCF file |
| format | Format of VCF variable columns |
| sample_1 | Reference sample VCF variable values corresponding with format |
| sample_2 | Tumour sample VCF variable values corresponding with format |
| tissue | Lung tumour tissue type: `Squamous` or `Adenocarcinoma` |

## Supplementary Data S4: NSCLC HLA Loss of Heterozygosity

| Column Variable | Description                         |
|-----------------|-------------------------------------|
| accel_id        | CRUK Accelerator patient identifier |
| message         | HLA LOH detection message           |
| HLA_A_type1     | HLA allele type 1: A, B or C        |
| HLA_A_type2     | HLA allele type 2: A, B or C        |
| LossAllele      | HLA allele that was lost            |
| KeptAllele      | HLA allele that was retained        |

## Supplementary Data S5: NSCLC Shared Protein Lists

| Column Variable | Description |
|----|----|
| hla_class | HLA class (I or II) |
| intersection | Set intersection identifier |
| n_proteins | Number of patients observations of protein in set intersection |
| uniprot_id | UniProt protein identifier |
| gene_name | HGNC gene name |

## Supplementary Data S6: NSCLC pVACseq Class I Predictions

| Column Variable | Description |
|-----------------------------------------|-------------------------------|
| sample | CRUK Accelerator patient identifier |
| chromosome | The chromosome of this variant |
| start | The start position of this variant in the zero-based, half-open coordinate system |
| stop | The stop position of this variant in the zero-based, half-open coordinate system |
| reference | The reference allele |
| variant | The alt allele |
| transcript | The Ensembl ID of the affected transcript |
| transcript_support_level | The transcript support level (TSL) of the affected transcript. `NA` if the VCF entry doesn't contain TSL information |
| ensembl_gene_id | The Ensembl ID of the affected gene |
| variant_type | The type of variant. `missense` for missense mutations, `inframe_ins` for inframe insertions, `inframe_del` for inframe deletions, and `FS` for frameshift variants |
| mutation | The amino acid change of this mutation |
| protein_position | The protein position of the mutation |
| gene_name | The Ensembl gene name of the affected gene |
| hgv_sc | The HGVS coding sequence variant name |
| hgv_sp | The HGVS protein sequence variant name |
| hla_allele | The HLA allele for this prediction |
| peptide_length | The peptide length of the epitope |
| sub_peptide_position | The one-based position of the epitope within the protein sequence used to make the prediction |
| mutation_position | The one-based position of the start of the mutation within the epitope sequence. 0 if the start of the mutation is before the epitope |
| mt_epitope_seq | The mutant epitope sequence |
| wt_epitope_seq | The wildtype (reference) epitope sequence at the same position in the full protein sequence. `NA` if there is no wildtype sequence at this position or if more than half of the amino acids of the mutant epitope are mutated |
| best_mt_score_method | Prediction algorithm with the lowest mutant ic50 binding affinity for this epitope |
| best_mt_score | Lowest ic50 binding affinity of all prediction algorithms used |
| corresponding_wt_score | ic50 binding affinity of the wildtype epitope. `NA` if there is no WT Epitope Seq |
| corresponding_fold_change | Corresponding WT Score / Best MT Score. `NA` if there is no WT Epitope Seq |
| tumor_dna_depth | Tumor DNA depth at this position. `NA` if VCF entry does not contain tumor DNA readcount annotation |
| tumor_dna_vaf | Tumor DNA variant allele frequency (VAF) at this position. `NA` if VCF entry does not contain tumor DNA readcount annotation |
| tumor_rna_depth | Tumor RNA depth at this position. `NA` if VCF entry does not contain tumor RNA readcount annotation |
| tumor_rna_vaf | Tumor RNA variant allele frequency (VAF) at this position. `NA` if VCF entry does not contain tumor RNA readcount annotation |
| normal_depth | Normal DNA depth at this position. `NA` if VCF entry does not contain normal DNA readcount annotation |
| normal_vaf | Normal DNA variant allele frequency (VAF) at this position. `NA` if VCF entry does not contain normal DNA readcount annotation |
| gene_expression | Gene expression value for the annotated gene containing the variant. `NA` if VCF entry does not contain gene expression annotation |
| transcript_expression | Transcript expression value for the annotated transcript containing the variant. `NA` if VCF entry does not contain transcript expression annotation |
| median_mt_score | Median ic50 binding affinity of the mutant epitope across all prediction algorithms used |
| median_wt_score | Median ic50 binding affinity of the wildtype epitope across all prediction algorithms used. `NA` if there is no WT Epitope Seq |
| median_fold_change | Median WT Score / Median MT Score. `NA` if there is no WT Epitope Seq |
| mh_cflurry_wt_score | MHCflurry ic50 binding affinity for the wildtype epitope |
| mh_cflurry_mt_score | MHCflurry ic50 binding affinity for the mutant epitope |
| mh_cnuggets_i_wt_score | MHCnuggets-I ic50 binding affinity for the wildtype epitope |
| mh_cnuggets_i_mt_score | MHCnuggets-I ic50 binding affinity for the mutant epitope |
| net_mhc_wt_score | NetMHC ic50 binding affinity for the wildtype epitope |
| net_mhc_mt_score | NetMHC ic50 binding affinity for the mutant epitope |
| pick_pocket_wt_score | PickPocket ic50 binding affinity for the wildtype epitope |
| pick_pocket_mt_score | PickPocket ic50 binding affinity for the mutant epitope |
| cterm_7mer_gravy_score | Mean hydropathy of last 7 residues on the C-terminus of the peptide |
| max_7mer_gravy_score | Max GRAVY score of any kmer in the amino acid sequence. Used to determine if there are any extremely hydrophobic regions within a longer amino acid sequence |
| difficult_n_terminal_residue | Is N-terminal amino acid a Glutamine, Glutamic acid, or Cysteine? (T/F) |
| c_terminal_cysteine | Is the C-terminal amino acid a Cysteine? (T/F) |
| c_terminal_proline | Is the C-terminal amino acid a Proline? (T/F) |
| cysteine_count | Number of Cysteines in the amino acid sequence. Problematic because they can form disulfide bonds across distant parts of the peptide |
| n_terminal_asparagine | Is the N-terminal amino acid an Asparagine? (T/F) |
| asparagine_proline_bond_count | Number of Asparagine-Proline bonds. Problematic because they can spontaneously cleave the peptide |
| b_rank | Rank of binding score: 1/median neoantigen binding affinity. Lower is better |
| f_rank | Rank of fold change: the difference in median binding affinity between neoantigen and wildtype peptide (agretopicity). Higher is better |
| m_rank | Ranks of mutant allele expression: the product of gene_expression and tumor_rna_vaf. Higher is better |
| d_rank | Rank of the tumor_dna_vaf. Higher is better |
| score | A score is calculated from the above ranks with the following formula: b_rank + f_rank + (m_rank * 2) + (d_rank/2). Higher is better |
| rank_score | The score converted to a rank, with the best being 1, splitting ties by first. Lower is better |
| rank_percent | The percentage rank score. Lower is better |

## Supplementary Data S7: NSCLC pVACseq Class II Predictions

| Column Variable | Description |
|-----------------------------------------|-------------------------------|
| sample | CRUK Accelerator patient identifier |
| chromosome | The chromosome of this variant |
| start | The start position of this variant in the zero-based, half-open coordinate system |
| stop | The stop position of this variant in the zero-based, half-open coordinate system |
| reference | The reference allele |
| variant | The alt allele |
| transcript | The Ensembl ID of the affected transcript |
| transcript_support_level | The transcript support level (TSL) of the affected transcript. `NA` if the VCF entry doesn't contain TSL information |
| ensembl_gene_id | The Ensembl ID of the affected gene |
| variant_type | The type of variant. `missense` for missense mutations, `inframe_ins` for inframe insertions, `inframe_del` for inframe deletions, and `FS` for frameshift variants |
| mutation | The amino acid change of this mutation |
| protein_position | The protein position of the mutation |
| gene_name | The Ensembl gene name of the affected gene |
| hgv_sc | The HGVS coding sequence variant name |
| hgv_sp | The HGVS protein sequence variant name |
| hla_allele | The HLA allele for this prediction |
| peptide_length | The peptide length of the epitope |
| sub_peptide_position | The one-based position of the epitope within the protein sequence used to make the prediction |
| mutation_position | The one-based position of the start of the mutation within the epitope sequence. 0 if the start of the mutation is before the epitope |
| mt_epitope_seq | The mutant epitope sequence |
| wt_epitope_seq | The wildtype (reference) epitope sequence at the same position in the full protein sequence. `NA` if there is no wildtype sequence at this position or if more than half of the amino acids of the mutant epitope are mutated |
| best_mt_score_method | Prediction algorithm with the lowest mutant ic50 binding affinity for this epitope |
| best_mt_score | Lowest ic50 binding affinity of all prediction algorithms used |
| corresponding_wt_score | ic50 binding affinity of the wildtype epitope. `NA` if there is no WT Epitope Seq |
| corresponding_fold_change | Corresponding WT Score / Best MT Score. `NA` if there is no WT Epitope Seq |
| tumor_dna_depth | Tumor DNA depth at this position. `NA` if VCF entry does not contain tumor DNA readcount annotation |
| tumor_dna_vaf | Tumor DNA variant allele frequency (VAF) at this position. `NA` if VCF entry does not contain tumor DNA readcount annotation |
| tumor_rna_depth | Tumor RNA depth at this position. `NA` if VCF entry does not contain tumor RNA readcount annotation |
| tumor_rna_vaf | Tumor RNA variant allele frequency (VAF) at this position. `NA` if VCF entry does not contain tumor RNA readcount annotation |
| normal_depth | Normal DNA depth at this position. `NA` if VCF entry does not contain normal DNA readcount annotation |
| normal_vaf | Normal DNA variant allele frequency (VAF) at this position. `NA` if VCF entry does not contain normal DNA readcount annotation |
| gene_expression | Gene expression value for the annotated gene containing the variant. `NA` if VCF entry does not contain gene expression annotation |
| transcript_expression | Transcript expression value for the annotated transcript containing the variant. `NA` if VCF entry does not contain transcript expression annotation |
| median_mt_score | Median ic50 binding affinity of the mutant epitope across all prediction algorithms used |
| median_wt_score | Median ic50 binding affinity of the wildtype epitope across all prediction algorithms used. `NA` if there is no WT Epitope Seq |
| median_fold_change | Median WT Score / Median MT Score. `NA` if there is no WT Epitope Seq |
| mh_cnuggets_ii_wt_score | MHCnuggets-II ic50 binding affinity for the wildtype epitope |
| mh_cnuggets_ii_mt_score | MHCnuggets-II ic50 binding affinity for the mutant epitope |
| cterm_7mer_gravy_score | Mean hydropathy of last 7 residues on the C-terminus of the peptide |
| max_7mer_gravy_score | Max GRAVY score of any kmer in the amino acid sequence. Used to determine if there are any extremely hydrophobic regions within a longer amino acid sequence |
| difficult_n_terminal_residue | Is N-terminal amino acid a Glutamine, Glutamic acid, or Cysteine? (T/F) |
| c_terminal_cysteine | Is the C-terminal amino acid a Cysteine? (T/F) |
| c_terminal_proline | Is the C-terminal amino acid a Proline? (T/F) |
| cysteine_count | Number of Cysteines in the amino acid sequence. Problematic because they can form disulfide bonds across distant parts of the peptide |
| n_terminal_asparagine | Is the N-terminal amino acid an Asparagine? (T/F) |
| asparagine_proline_bond_count | Number of Asparagine-Proline bonds. Problematic because they can spontaneously cleave the peptide |
| net_mhci_ipan_wt_score | NetMHCIIpan ic50 binding affinity for the wildtype epitope |
| net_mhci_ipan_mt_score | NetMHCIIpan ic50 binding affinity for the mutant epitope |
| n_nalign_wt_score | NNalign ic50 binding affinity for the wildtype epitope |
| n_nalign_mt_score | NNalign ic50 binding affinity for the mutant epitope |
| sm_malign_wt_score | SMMalign ic50 binding affinity for the wildtype epitope |
| sm_malign_mt_score | SMMalign ic50 binding affinity for the mutant epitope |
| b_rank | Rank of binding score: 1/median neoantigen binding affinity. Lower is better |
| f_rank | Rank of fold change: the difference in median binding affinity between neoantigen and wildtype peptide (agretopicity). Higher is better |
| m_rank | Ranks of mutant allele expression: the product of gene_expression and tumor_rna_vaf. Higher is better |
| d_rank | Rank of the tumor_dna_vaf. Higher is better |
| score | A score is calculated from the above ranks with the following formula: b_rank + f_rank + (m_rank * 2) + (d_rank/2). Higher is better |
| rank_score | The score converted to a rank, with the best being 1, splitting ties by first. Lower is better |
| rank_percent | The percentage rank score. Lower is better |

## Supplementary Data S8: NSCLC pVACseq Peptidome Combined Predictions

| Column Variable | Description |
|-----------------------------------------|-------------------------------|
| sample_id | CRUK Accelerator patient identifier |
| hla_class | HLA class (I or II) |
| table_name | Identifier in the form `accel_id`/`predicted_hla_allotype`/`peptide_length` |
| length | Peptide length |
| hla_length_pref1 | HLA allotype preferred peptide length 1 |
| hla_length_pref2 | HLA allotype preferred peptide length 2 |
| gene_name | HGNC gene name |
| mt_epitope_seq | The mutant epitope sequence |
| median_mt_score | Median ic50 binding affinity of the mutant epitope across all prediction algorithms used |
| corresponding_fold_change | `Corresponding WT Score` / `Best MT Score`. `NA` if there is no `WT Epitope Seq` |
| gene_expression | Gene expression value for the annotated gene containing the variant |
| tumor_rna_vaf | Tumor RNA variant allele frequency (VAF) at this position |
| tumor_dna_vaf | Tumor DNA variant allele frequency (VAF) at this position |
| rank_percent | The percentage rank score. Lower is better |
| rank_score | The `score` converted to a rank, with the best being 1, splitting ties by first. Lower is better |
| Obs_I | The number of peptides from the source protein observed by mass spectrometry in HLA-I immunopeptidome |
| Obs_II | The number of peptides from the source protein observed by mass spectrometry in HLA-II immunopeptidome |

## Supplementary Data S9: NSCLC Tested Neoantigens

| Column Variable | Description |
|-----------------------------------------|-------------------------------|
| accel_id | CRUK Accelerator patient identifier |
| gene_name | Gene |
| mt_epitope_seq | Mutated (neoantigen) peptide sequence |
| wt_epitope_seq | Wildtype peptide sequence |
| peptide_length | Peptide length |
| table_name | Identifier in the form `accel_id`/`predicted_hla_allotype`/`peptide_length` e.g. `A119/DRB1*04:04/15` |
| mutation | The mutation `From/To` |
| protein_position | Location of the mutation in the source protein, UNIPROT sequence number |
| Obs_I | The number of peptides from the source protein observed by mass spectrometry in HLA-I immunopeptidome |
| Obs_II | The number of peptides from the source protein observed by mass spectrometry in HLA-II immunopeptidome |
| median_mt_score | The median pVACseq predicted binding affinity of the neoantigen peptide |
| median_wt_score | The median pVACseq predicted binding affinity of the wildtype peptide |
| median_fold_change | The ratio between the median neoantigen affinity and wildtype peptide affinity |
| rank_percent | The overall rank percentage for the neoantigen from pVACseq for the peptide of that length and HLA allotype |
| mean_sfc_mt | Mean IFN-γ ELISPOT spot forming cells per million cells for the neoantigen peptide |
| mean_sfc_wt | Mean IFN-γ ELISPOT spot forming cells per million cells for the wildtype peptide |
| elispot_response | ELISPOT response category: Strong, Weak or None |

## References

<span class="csl-left-margin">1.
</span><span class="csl-right-inline">Nicholas B, Bailey A, McCann KJ,
Wood O, Currall E, Johnson P, et al. Proteogenomics guided
identification of functional neoantigens in non-small cell lung cancer.
2024. Available: <http://dx.doi.org/10.1101/2024.05.30.596609></span>

<span class="csl-left-margin">2.
</span><span class="csl-right-inline">Hundal J, Carreno BM, Petti AA,
Linette GP, Griffith OL, Mardis ER, et al. pVAC-seq: A genome-guided in
silico approach to identifying tumor neoantigens. Genome Medicine.
2016;8: 11.
doi:[10.1186/s13073-016-0264-5](https://doi.org/10.1186/s13073-016-0264-5)</span>
