# Fine-Mapping for Samples of Related Individuals

Most fine-mapping methods (e.g., FINEMAP and SuSiE) are designed for samples of unrelated individuals, which can be problematic when dealing with related individuals, such as in farm animal populations. To address this, we previously developed BFMAP, a method utilizing individual-level data. BFMAP demonstrates higher power for detecting true causal mutations and lower false positive rates compared to existing methods when applied to related individuals.

Building on this work, we introduce two new methods, FINEMAP-adj and SuSiE-adj, which apply FINEMAP and SuSiE, respectively, by incorporating a relatedness-adjusted LD matrix for fine-mapping in samples of related individuals. Both methods utilize adjusted summary statistics and can achieve performance comparable to BFMAP-SSS.

## Methods
- **FINEMAP-adj**: Applies FINEMAP with relatedness-adjusted inputs
- **SuSiE-adj**: Applies SuSiE with relatedness-adjusted inputs  
- **BFMAP**: Individual-level fine-mapping method ([GitHub](https://github.com/jiang18/bfmap/))

## FINEMAP-adj and SuSiE-adj
The core of FINEMAP-adj and SuSiE-adj is the use of a relatedness-adjusted LD matrix, an adjusted sample size (effective sample size), and mixed-model association statistics.

While mixed models are the standard method for GWAS, special attention is required for the relatedness-adjusted LD matrix and effective sample size, which can be computed using our [LD Adjuster](ld_adjuster/) tool. 

These serve as standard inputs to [FINEMAP](http://www.christianbenner.com/) and [SuSiE](https://stephenslab.github.io/susieR/index.html) to enable FINEMAP-adj and SuSiE-adj, respectively.

### [LD Adjuster](ld_adjuster/): relatedness-adjusted LD matrix

### FINEMAP-adj: [FINEMAP](http://www.christianbenner.com/) with relatedness-adjusted inputs
`finemap --sss` is the major routine for FINEMAP-adj.

### SuSiE-adj: [SuSiE](https://stephenslab.github.io/susieR/index.html) with relatedness-adjusted inputs
[`susie_rss`](https://stephenslab.github.io/susieR/reference/susie_rss.html) is the major function for SuSiE-adj. Set `estimate_residual_variance = TRUE` as recommended.

## Gene PIP
The Posterior Inclusion Probability (PIP) for a gene can be calculated by summing the posterior probabilities of all models (or configs) that contain any variants within that gene. This calculation can be performed using results from BFMAP-SSS or FINEMAP/FINEMAP-adj.

#### Usage
We provide an R function to calculate gene PIPs from either BFMAP-SSS or FINEMAP/FINEMAP-adj results.
```r
# Load required packages
library(data.table)

# Source the function
source("calc_gene_pip.R")

# Calculate gene PIPs
genepip <- calc_gene_pip(gtf, pip, model)
```

#### Input
- `gtf`: A data.table or data.frame of Ensembl GTF data with columns: seqname, start, end, and attribute
- `pip`: A data.table containing SNP PIPs from BFMAP-SSS or FINEMAP/FINEMAP-adj
  - For BFMAP-SSS: read from `*.pip.csv` (requires columns: SNPname, Chr, and Pos)
  - For FINEMAP-adj: read from `*.snp` (requires columns: rsid, chromosome, and position)
- `model`: A data.table containing model configurations from BFMAP-SSS or FINEMAP/FINEMAP-adj
  - For BFMAP-SSS: read from `*.model.csv`
  - For FINEMAP-adj: read from `*.config` (requires columns: config and prob)

#### Output
Returns a data.frame with columns: Chr, Start, End, PIP, and Attribute, sorted by descending PIP value.

## [Example](https://github.com/JJWang259/FineMapping-RelatedIndividuals/tree/main/example)
