# Fine-Mapping for Samples of Related Individuals

Most fine-mapping methods have been designed for samples of unrelated individuals, which can be problematic when dealing with related individuals, such as in farm animal populations. To address this, we previously developed BFMAP, a method utilizing individual-level data. BFMAP demonstrates higher power in detecting true causal mutations and lower false positive rates compared to existing methods when used with related individuals.
Building on this work, we introduce two new methods, FINEMAP-adj and SuSiE-adj, which apply FINEMAP and SuSiE, respectively, by incorporating a relatedness-adjusted genotype correlation matrix for fine-mapping in samples of related individuals. Both methods utilize the adapted summary statistics and can achieve performance comparable to BFMAP-SSS.

## Methods

- **FINEMAP-adj**: A computational workflow for applying FINEMAP to related individuals by incorporating adjusted correlation matrices.
- **SuSiE-adj**: A computational workflow for applying SuSiE to related individuals by incorporating adjusted correlation matrices.
- **BFMAP**: https://github.com/jiang18/bfmap/

## [BFMAP](https://github.com/jiang18/bfmap/)


## FINEMAP-adj and SuSiE-adj
The computational workflows FINEMAP-adj and SuSiE-adj apply FINEMAP and SuSiE to related individuals by incorporating adapted summary statistics.

### Relatedness-adjusted LD matrix

[LD Adjuster](ld_adjuster/)

### FINEMAP-adj
This is a computational workflow for applying [FINEMAP](http://www.christianbenner.com/) with adapted summary statistics.

``` bash
./finemap --sss --in-files <data> --dataset <num>
```
- `--dataset`: Specify the row number of datesets for fine-mapping in the data file. 
- `--in-files`: A semicolon-separated text file and could look as follows:

```plaintext
z;ld;snp;config;cred;log;n_samples
dataset.z;dataset.adj.ld;dataset.snp;dataset.config;dataset.cred;dataset.log;n_eff
```

- `z`: Contains the names of Z files (input)
- `ld`: Contains the names of relatedness-adjusted genotype correlation matrix (input)
- `snp`: Contains the names of SNP files (output)
- `config`: Contains the names of CONFIG files (output)
- `cred`: Contains the names of CRED files (output)
- `n_samples`: Contains the effective sample size.
- `k`: Contains the names of optional K files (optional input)
- `log`: Contains the names of optional LOG files (optional output)

The `dataset.z` file is a space-delimited text file and contains the GWAS summary statistics one SNP per line.

```plaintext
rsid chromosome position allele1 allele2 maf beta se
```

- **rsid**: The SNP identifier.
- **chromosome**: The chromosome number where the SNP is located.
- **position**: The position of the SNP on the chromosome.
- **allele1**: The  reference allele (the coded effect allele).
- **allele2**: The other allele.
- **maf**: The minor allele frequency.
- **beta**: The effect size estimate.
- **se**: The standard error of the effect size estimate.

### SuSiE-adj
This is a computational workflow for applying [SuSiE](https://stephenslab.github.io/susieR/index.html) with adapted summary statistics.

The `susieR` package needs to be installed in advance. 

Fine-mapping with susieR using adjusted summary statistics
``` R
fitted_rss <- susie_rss(bhat = betahat, shat = sebetahat, n = n_eff, R = R_adj, var_y = var(y), L = 10,
                         estimate_residual_variance = TRUE)
```
#### Parameters

- `bhat`: Vector of effect size estimates.
- `shat`: Vector of standard errors of the effect size estimates.
- `n`: Effective sample size.
- `R`: Relatedness-adjusted genotype correlation matrix
- `var_y`: Variance of the phenotype.
- `L`: Maximum number of causal variants to consider (default is 10).
- `estimate_residual_variance`: Boolean flag to estimate residual variance (default is `FALSE`).

## Gene PIP
The Posterior Inclusion Probability (PIP) for a gene can be calculated by summing the posterior probabilities of all models that include any variants within that gene. This calculation can be performed using the results from BFMAP-SSS and FINEMAP-adj.

#### Usage
We provide a unified R function that can calculate gene PIP from either BFMAP-SSS or FINEMAP-adj results.
```r
# Load required packages
library(data.table)
library(dplyr)

# Source the function
source("genepip.R")

# Calculate gene PIP
genepip <- calc_gene_pip(genes, pip, model, method)
```

#### Input
- `genes`: A data.table or data.frame containing gene annotations, with columns: chr, start, end, attributes.
- `pip`: A data.table containing SNP PIP output from BFMAP-SSS or FINEMAP-adj.
  - For BFMAP-SSS: read from `*.pip.csv`, requires columns Pos, SNPname, Chr
  - For FINEMAP-adj: read from `*.snp`, requires columns position, rsid, chromosome
- `model`: A data.table containing model configurations output from BFMAP-SSS or FINEMAP-adj.
  - For BFMAP-SSS: read from `*.model.csv`, with one column per SNP in the model and the last column containing PIP values
  - For FINEMAP-adj: read from `*.config`, requires columns 'config' (comma-separated SNP names) and 'prob' (probability)
- `method`: Which method's output format is being used. Options: "BFMAP-SSS" or "FINEMAP-adj".

#### Output
Returns a data.frame with columns: Chr, Start, End, PIP, Attributes, sorted by descending PIP value.

## [Example](https://github.com/JJWang259/FineMapping-RelatedIndividuals/tree/main/example)
