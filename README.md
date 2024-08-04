
Last Updated on 2024-08-04

# Supplementary Materials

This repository contains the supplementary tables as `csv` files for:
**Proteogenomics guided identification of functional neoantigens in
non-small cell lung cancer** \[[1](#ref-nicholas2024)\]. This repository
is associated with
[![DOI](https://zenodo.org/badge/833647778.svg)](https://zenodo.org/doi/10.5281/zenodo.12820423)

The column names and contents of the `csv` files in the `tables` folder
are described below.

## Supplementary Material 1: Patient information

Supplementary Material 1 is Table S1, a `csv` file containing patient
information with 24 rows and 21 column variables. Each row in Table S1
represents observations for a single patient.

<a href="#tbl-supp-01" class="quarto-xref">Table 1</a> provides
descriptions of the values contained in each column of Table S1.

<table>
<colgroup>
<col style="width: 24%" />
<col style="width: 75%" />
</colgroup>
<thead>
<tr class="header">
<th>Column name</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>accel_id</code></td>
<td>CRUK Accelerator patient identifier</td>
</tr>
<tr class="even">
<td><code>target_lung_id</code></td>
<td>Targeted Lung Health Check patient identifier</td>
</tr>
<tr class="odd">
<td><code>tissue</code></td>
<td>NSCLC type: Adenocarcinoma or Squamous cell carcinoma</td>
</tr>
<tr class="even">
<td><code>n_somatic_variants</code></td>
<td>Total number of somatic variants identified by whole exome
sequencing</td>
</tr>
<tr class="odd">
<td><code>mut_burden_per_mb</code></td>
<td>Mutational burden: mutations per million bases of DNA.<br />
Exome target size was 35.7 Mb</td>
</tr>
<tr class="even">
<td><code>obs_class_I</code></td>
<td>Number of observed HLA I peptides by mass spec.
immunopeptidomics</td>
</tr>
<tr class="odd">
<td><code>obs_class_II</code></td>
<td>Number of observed HLA II peptides by mass spec.
immunopeptidomics</td>
</tr>
<tr class="even">
<td><code>HLA</code></td>
<td>Class I and II HLA allotypes identified by genomic sequencing</td>
</tr>
<tr class="odd">
<td><code>wet_weight</code></td>
<td>Wet weight of tumour tissue</td>
</tr>
<tr class="even">
<td><code>tumour_purity</code></td>
<td>Tumour purity as calculated from WES by ASCAT</td>
</tr>
<tr class="odd">
<td><code>tumour_ploidy</code></td>
<td>Tumour ploidy as calculated from WES by ASCAT</td>
</tr>
<tr class="even">
<td><code>til_status</code></td>
<td>Tumour infiltrating T-cell status by immunohistochemistry: Low,
Moderate, High or NA</td>
</tr>
<tr class="odd">
<td><code>weeks_post_surgery</code></td>
<td>Number of weeks since surgery</td>
</tr>
<tr class="even">
<td><code>status_as_of_2021_01_19</code></td>
<td>Status since 2021-01-19: Alive, Deceased or NA</td>
</tr>
<tr class="odd">
<td><code>sex</code></td>
<td>Patient sex</td>
</tr>
<tr class="even">
<td><code>date_of_diagnosis</code></td>
<td>Date of diagnosis</td>
</tr>
<tr class="odd">
<td><code>age_at_diagnosis</code></td>
<td>Age at diagnosis</td>
</tr>
<tr class="even">
<td><code>smoking_status</code></td>
<td>Smoking status</td>
</tr>
<tr class="odd">
<td><code>notes_2</code></td>
<td>Notes about smoking history</td>
</tr>
</tbody>
</table>

## Supplementary Material 2: NSCLC mutations

Supplementary Material 2 is Table S2, a compressed `csv` file containing
all the mutations (variant calls) from the WES comparing tumour to
normal adjacent tissue. It has 106,285 rows with 16 columns comprising
the variants from 24 donors. Variant types are single nucleotide
variant, insertion, deletion and complex variant.
<a href="#tbl-supp-02" class="quarto-xref">Table 2</a> contains the
description of the column variables.

| Column variable | Description |
|-----------|-------------------------------------------------------------|
| `accel_id` | CRUK Accelerator patient identifier |
| `vid` | Unique variant identifier |
| `chrom` | Chromosome |
| `pos` | Genomic coordinate |
| `ref` | Reference base |
| `alt` | Variant base |
| `info` | Information field from VCF file |
| `format` | Format of VCF variable columns |
| `sample_1` | Reference sample VCF variable values corresponding with format |
| `sample_2` | Tumour sample VCF variable values corresponding with format |
| `type` | Variant type: `snv`, `ins` , `del` or `complex` . Single nucleotide variant, insertion, deletion and complex variant respectively |
| `ensembl` | Ensembl gene identifier |
| `gene_name` | HGNC gene name |
| `vaf` | Variant allele frequency |
| `tissue` | Lung tumour tissue type: `Squamous` or `Adenocarcinoma` |
| `cell_compartment` | Cell compartment of protein product of gene, |

## Supplementary Material 4 and 4: pVACseq predicted neoantigens

Supplementary Material 3 and 4 are Tables S3 and S4. These are `csv`
files containing all the pVACseq \[[2](#ref-hundal2016)\] predicted
neoantigen peptides and their wildtype equivalents,
<a href="#tbl-supp-03-04" class="quarto-xref">Table 3</a> contains
descriptions of the values contained in each column. Each row in Tables
S3 and S4 represents one set of predictions i.e. one mutation and
predicted neoantigen peptide per row.

Table S3 has 27,446 rows and 59 columns. Table S4 has 127,015 rows and
59 columns.

<table>
<colgroup>
<col style="width: 22%" />
<col style="width: 77%" />
</colgroup>
<thead>
<tr class="header">
<th>Column Name</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>sample</code></td>
<td>CRUK Accelerator patient identifier</td>
</tr>
<tr class="even">
<td><code>Chromosome</code></td>
<td>The chromosome of this variant</td>
</tr>
<tr class="odd">
<td><code>Start</code></td>
<td>The start position of this variant in the zero-based, half-open
coordinate system</td>
</tr>
<tr class="even">
<td><code>Stop</code></td>
<td>The stop position of this variant in the zero-based, half-open
coordinate system</td>
</tr>
<tr class="odd">
<td><code>Reference</code></td>
<td>The reference allele</td>
</tr>
<tr class="even">
<td><code>Variant</code></td>
<td>The alt allele</td>
</tr>
<tr class="odd">
<td><code>Transcript</code></td>
<td>The Ensembl ID of the affected transcript</td>
</tr>
<tr class="even">
<td><code>Transcript Support Level</code></td>
<td>The <a
href="https://useast.ensembl.org/info/genome/genebuild/transcript_quality_tags.html#tsl">transcript
support level (TSL)</a> of the affected transcript. <code>NA</code> if
the VCF entry doesn’t contain TSL information.</td>
</tr>
<tr class="odd">
<td><code>Ensembl Gene ID</code></td>
<td>The Ensembl ID of the affected gene</td>
</tr>
<tr class="even">
<td><code>Variant Type</code></td>
<td>The type of variant. <code>missense</code> for missense mutations,
<code>inframe_ins</code> for inframe insertions,
<code>inframe_del</code> for inframe deletions, and <code>FS</code> for
frameshift variants</td>
</tr>
<tr class="odd">
<td><code>Mutation</code></td>
<td>The amnio acid change of this mutation</td>
</tr>
<tr class="even">
<td><code>Protein Position</code></td>
<td>The protein position of the mutation</td>
</tr>
<tr class="odd">
<td><code>Gene Name</code></td>
<td>The Ensembl gene name of the affected gene</td>
</tr>
<tr class="even">
<td><code>HGVSc</code></td>
<td>The HGVS coding sequence variant name</td>
</tr>
<tr class="odd">
<td><code>HGVSp</code></td>
<td>The HGVS protein sequence variant name</td>
</tr>
<tr class="even">
<td><code>HLA Allele</code></td>
<td>The HLA allele for this prediction</td>
</tr>
<tr class="odd">
<td><code>Peptide Length</code></td>
<td>The peptide length of the epitope</td>
</tr>
<tr class="even">
<td><code>Sub-peptide Position</code></td>
<td>The one-based position of the epitope within the protein sequence
used to make the prediction</td>
</tr>
<tr class="odd">
<td><code>Mutation Position</code></td>
<td>The one-based position of the start of the mutation within the
epitope sequence. <code>0</code> if the start of the mutation is before
the epitope</td>
</tr>
<tr class="even">
<td><code>MT Epitope Seq</code></td>
<td>The mutant epitope sequence</td>
</tr>
<tr class="odd">
<td><code>WT Epitope Seq</code></td>
<td>The wildtype (reference) epitope sequence at the same position in
the full protein sequence. <code>NA</code> if there is no wildtype
sequence at this position or if more than half of the amino acids of the
mutant epitope are mutated</td>
</tr>
<tr class="even">
<td><code>Best MT Score Method</code></td>
<td>Prediction algorithm with the lowest mutant ic50 binding affinity
for this epitope</td>
</tr>
<tr class="odd">
<td><code>Best MT Score</code></td>
<td>Lowest ic50 binding affinity of all prediction algorithms used</td>
</tr>
<tr class="even">
<td><code>Corresponding WT Score</code></td>
<td>ic50 binding affinity of the wildtype epitope. <code>NA</code> if
there is no <code>WT Epitope Seq</code>.</td>
</tr>
<tr class="odd">
<td><code>Corresponding Fold Change</code></td>
<td><code>Corresponding WT Score</code> / <code>Best MT Score</code>.
<code>NA</code> if there is no <code>WT Epitope Seq</code>.</td>
</tr>
<tr class="even">
<td><code>Best MT Percentile Method</code></td>
<td>Prediction algorithm with the lowest binding affinity percentile
rank for this epitope</td>
</tr>
<tr class="odd">
<td><code>Best MT Percentile</code></td>
<td>Lowest percentile rank of this epitope’s ic50 binding affinity of
all prediction algorithms used (those that provide percentile
output)</td>
</tr>
<tr class="even">
<td><code>Corresponding WT Percentile</code></td>
<td>binding affinity percentile rank of the wildtype epitope.
<code>NA</code> if there is no <code>WT Epitope Seq</code>.</td>
</tr>
<tr class="odd">
<td><code>Tumor DNA Depth</code></td>
<td>Tumor DNA depth at this position. <code>NA</code> if VCF entry does
not contain tumor DNA readcount annotation.</td>
</tr>
<tr class="even">
<td><code>Tumor DNA VAF</code></td>
<td>Tumor DNA variant allele frequency (VAF) at this position.
<code>NA</code> if VCF entry does not contain tumor DNA readcount
annotation.</td>
</tr>
<tr class="odd">
<td><code>Tumor RNA Depth</code></td>
<td>Tumor RNA depth at this position. <code>NA</code> if VCF entry does
not contain tumor RNA readcount annotation.</td>
</tr>
<tr class="even">
<td><code>Tumor RNA VAF</code></td>
<td>Tumor RNA variant allele frequency (VAF) at this position.
<code>NA</code> if VCF entry does not contain tumor RNA readcount
annotation.</td>
</tr>
<tr class="odd">
<td><code>Normal Depth</code></td>
<td>Normal DNA depth at this position. <code>NA</code> if VCF entry does
not contain normal DNA readcount annotation.</td>
</tr>
<tr class="even">
<td><code>Normal VAF</code></td>
<td>Normal DNA variant allele frequency (VAF) at this position.
<code>NA</code> if VCF entry does not contain normal DNA readcount
annotation.</td>
</tr>
<tr class="odd">
<td><code>Gene Expression</code></td>
<td>Gene expression value for the annotated gene containing the variant.
<code>NA</code> if VCF entry does not contain gene expression
annotation.</td>
</tr>
<tr class="even">
<td><code>Transcript Expression</code></td>
<td>Transcript expression value for the annotated transcript containing
the variant. <code>NA</code> if VCF entry does not contain transcript
expression annotation.</td>
</tr>
<tr class="odd">
<td><code>Median MT Score</code></td>
<td>Median ic50 binding affinity of the mutant epitope across all
prediction algorithms used</td>
</tr>
<tr class="even">
<td><code>Median WT Score</code></td>
<td>Median ic50 binding affinity of the wildtype epitope across all
prediction algorithms used. <code>NA</code> if there is no
<code>WT Epitope Seq</code>.</td>
</tr>
<tr class="odd">
<td><code>Median Fold Change</code></td>
<td><code>Median WT Score</code> / <code>Median MT Score</code>.
<code>NA</code> if there is no <code>WT Epitope Seq</code>.</td>
</tr>
<tr class="even">
<td><code>Individual Prediction Algorithm WT and MT Scores</code>
(multiple)</td>
<td><p>ic50 binding affintity for the <code>MT Epitope Seq</code> and
<code>WT Eptiope Seq</code> for the individual prediction algorithms
used.</p>
<p>Four binding algorithms were used for class I predictions (MHCflurry,
MHCnuggetsI, NNalign, NetMHC, PickPocket) and four for class II
predictions (MHCnuggetsII, NetMHCIIpan, NNalign, SMMalign).</p></td>
</tr>
<tr class="odd">
<td><code>cterm_7mer_gravy_score</code></td>
<td>Mean hydropathy of last 7 residues on the C-terminus of the
peptide</td>
</tr>
<tr class="even">
<td><code>max_7mer_gravy_score</code></td>
<td>Max GRAVY score of any kmer in the amino acid sequence. Used to
determine if there are any extremely hydrophobic regions within a longer
amino acid sequence.</td>
</tr>
<tr class="odd">
<td><code>difficult_n_terminal_residue</code> (T/F)</td>
<td>Is N-terminal amino acid a Glutamine, Glutamic acid, or
Cysteine?</td>
</tr>
<tr class="even">
<td><code>c_terminal_cysteine</code> (T/F)</td>
<td>Is the C-terminal amino acid a Cysteine?</td>
</tr>
<tr class="odd">
<td><code>c_terminal_proline</code> (T/F)</td>
<td>Is the C-terminal amino acid a Proline?</td>
</tr>
<tr class="even">
<td><code>cysteine_count</code></td>
<td>Number of Cysteines in the amino acid sequence. Problematic because
they can form disulfide bonds across distant parts of the peptide</td>
</tr>
<tr class="odd">
<td><code>n_terminal_asparagine</code> (T/F)</td>
<td>Is the N-terminal amino acid an Asparagine?</td>
</tr>
<tr class="even">
<td><code>asparagine_proline_bond_count</code></td>
<td>Number of Asparagine-Proline bonds. Problematic because they can
spontaneously cleave the peptide</td>
</tr>
<tr class="odd">
<td><code>b_rank</code></td>
<td>Rank of binding score: 1/median neoantigen binding affinity . Lower
is better</td>
</tr>
<tr class="even">
<td><code>f_rank</code></td>
<td>Rank of fold change: the difference in median binding affinity
between neoantigen and wildtype peptide (agretopicity). Higher is
better.</td>
</tr>
<tr class="odd">
<td><code>m_rank</code></td>
<td>Ranks of mutant allele expression: the product of
<code>gene_expression</code> and <code>tumor_rna_vaf</code> . Higher is
better.</td>
</tr>
<tr class="even">
<td><code>d_rank</code></td>
<td>Rank of the <code>tumor_dna_vaf</code> . Higher is better.</td>
</tr>
<tr class="odd">
<td><code>score</code></td>
<td>A <code>score</code> is calculated from the above ranks with the
following formula:
<code>b_rank + f_rank + (m_rank * 2) + (d_rank/2)</code> . Higher is
better</td>
</tr>
<tr class="even">
<td><code>rank_score</code></td>
<td>The <code>score</code> converted to a rank, with the best being 1,
splitting ties by first. Lower is better</td>
</tr>
<tr class="odd">
<td><code>rank_percent</code></td>
<td>The percentage rank score. Lower is better.</td>
</tr>
</tbody>
</table>

## Supplementary Material 5: Tested neoantigens

Supplementary Material 5 is Table S5, a `csv` file with 70 rows and 17
column variables for the neoantigen peptide predictions tested by IFN-𝛄
ELISPOT using autologous PBMCs. Each row in Table S4 represents one
neoantigen peptide and its wildtype equivalent and
<a href="#tbl-supp-05" class="quarto-xref">Table 4</a> contains
descriptions of the values contained in each column of Table S5.

|  |  |
|-------------|-----------------------------------------------------------|
| **Column name** | **Description** |
| `accel_id` | CRUK Accelerator patient identifier |
| `gene_name` | Gene |
| `mt_epitope_seq` | Mutated (neoantigen) peptide sequeunce |
| `wt_epitope_seq` | Wildtype peptide sequence |
| `peptide_length` | Peptide length |
| `table_name` | Identifier in the form `accel_id` / `predicted_hla_allotype` / `peptide_length` e.g. `A119/DRB1*04:04/15` |
| `mutation` | The mutation `From/To` |
| `protein_position` | Location of the mutation in the source protein, UNIPROT sequence number. |
| `Obs_I` | The number of peptides from the source protein observed by mass spectrometry observed in HLA-I immunopeptidome |
| `Obs_II` | The number of peptides from the source protein observed by mass spectrometry observed in HLA-II immunopeptidome |
| `median_mt_score` | The median pVACseq predicted binding affinity of the neoantigen peptide |
| `median_wt_score` | The median pVACseq predicted binding affinity of the wildtype peptide |
| `median_fold_change` | The ratio between the median neoantigen affinity and wildtype peptide affinity |
| `rank_percent` | The overall rank percentage for the neoantigen from pVACseq for the peptide of that length and HLA allotype. |
| `mean_sfc_mt` | Mean IFN-𝛄 ELISPOT spot forming cells per million cells for the neoantigen peptide |
| `mean_sfc_wt` | Mean IFN-𝛄 ELISPOT spot forming cells per million cells for the wildtype peptide |
| `elispot_response` | ELISPOT response category: Strong, Weak or None |

### **Table S6** List of patient samples selected for single-cell RNA and TCR sequencing and TotalSeq C antibodies (Biolegend).

|                              |                           |                     |
|---------------------------|-------------------------|--------------------|
| **Patient ID and condition** | **TotalSeq C Hashtag ID** | **Hashtag barcode** |
| A119_PTPRT-12\_**MUT**       | C0255                     | AAGTATCGTTTCGCA     |
| A119_PTPRT-12\_**WT**        | C0256                     | GGTTGCCAGATGTCA     |

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
