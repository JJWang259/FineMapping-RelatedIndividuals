# Fine-Mapping for Samples of related individuals

Most fine-mapping methods have been designed for samples of unrelated individuals, which can be problematic when dealing with related individuals, such as in farm animal populations. To address this, we previously developed BFMAP, a method utilizing individual-level data. BFMAP demonstrates higher power in detecting true causal mutations and lower false positive rates compared to existing methods when used with related individuals. 
Building on this work, we introduce two new methods, FINEMAP-Adj and SuSiE-Adj, which extend FINEMAP and SuSiE, respectively, by incorporating a relatedness-adjusted genotype correlation matrix for fine-mapping in samples of related individuals. Both methods utilize summary statistics and can achieve performance comparable to BFMAP-SSS.

## Methods


- **FINEMAP-adj**: An adaptation of FINEMAP with a relatedness-adjusted genotype correlation matrix.
- **SuSiE-adj**: An adaptation of SuSiE with a relatedness-adjusted genotype correlation matrix.
- **BFMAP**: (https://github.com/jiang18/bfmap/)

## [BFMAP](https://github.com/jiang18/bfmap/)


## FINEMAP-adj and SuSiE-adj
### Relatedness-adjusted genotype correlation matrix


```
# R script to generate Relatedness-adjusted genotype correlation matrix
ld_adj <- function(raw, h2, G) {
  eta <- h2 / (1 - h2)
  diag(G) <- eta * diag(G) + 1
  L <- t(chol(G))
  X <- as.matrix(raw)
  
  # Efficient scaling and centering
  X <- scale(X, center = TRUE, scale = TRUE)
  
  # Efficient solve
  Xstar <- solve(L, X)
  Rstar <- crossprod(Xstar)  
  
  n_eff <- mean(diag(Rstar))
  print(paste0("Mean effective sample size is: ", n_eff))
  
  Ds <- sqrt(diag(Rstar))
  D_sol <- diag(1 / Ds, nrow = length(Ds), ncol = length(Ds))
  R_adj <- D_sol %*% Rstar %*% D_sol
  
  return(list(ld_matrix = R_adj, n_eff = n_eff))
}
```
#### Input
- `raw`: A matrix of raw genotype data where rows represent individuals and columns represent SNPs.
- `h2`: A numeric value representing the heritability estimate.
- `G`: A genetic relationship matrix (GRM).

### FINEMAP-Adj
This is an extension of [FINEMAP](http://www.christianbenner.com/) with adjusted summary statistics.

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

### SuSiE-Adj
This is an extension of [SuSiE](https://stephenslab.github.io/susieR/index.html) with adjusted summary statistics.

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
The Posterior Inclusion Probability (PIP) for a gene can be calculated by summing the posterior probabilities of all models that include any variants within that gene. This calculation can be performed using the results from BFMAP-SSS.

#### Usage

To calculate the Gene PIP, run the following command:
```bash
Rscript gene_ppc.R --chr <num> --gtf <gft file> --m <model file> --p <pip file> --o <prefix for output>
```

#### Input

- `--chr <num>`: Chromosome number.
- `--gtf <gtf file>`: Path to the GTF file containing gene annotations.
- `--m <model file>`: Path to the model file generated by BFMAP-SSS.
- `--p <pip file>`: Path to the PIP file containing posterior inclusion probabilities.
- `--o <prefix for output>`: Prefix for the output files.

## [Example](https://github.com/JJWang259/FineMapping-RelatedIndividuals/tree/main/Example)
