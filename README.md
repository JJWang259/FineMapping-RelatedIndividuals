# Fine-Mapping for Farm Animal Populations

Most fine-mapping methods have been designed for samples of unrelated individuals, which can be problematic when dealing with related individuals, such as in farm animal populations. To address this, we previously developed BFMAP, a method utilizing individual-level data. BFMAP demonstrates higher power in detecting true causal mutations and lower false positive rates compared to existing methods when used with related individuals. 
Building on this work, we introduce two new methods, FINEMAP-Adj and SuSiE-Adj, which extend FINEMAP and SuSiE, respectively, by incorporating a relatedness-adjusted genotype correlation matrix for fine-mapping in samples of related individuals. Both methods utilize summary statistics and can achieve performance comparable to BFMAP-SSS.

## Key feature

- **BFMAP**: (https://github.com/jiang18/bfmap/)
- **FINEMAP-Adj**: An extension of FINEMAP with a relatedness-adjusted genotype correlation matrix.
- **SuSiE-Adj**: An extension of SuSiE with a relatedness-adjusted genotype correlation matrix.

## [BFMAP](https://github.com/jiang18/bfmap/)

## Relatedness-adjusted genotype correlation matrix


## FINEMAP-Adj

## SuSiE-Adj

## Gene PIP
The PIP for a gene can be calculated by summing the posterior probabilities of all models that include any variants within that gene. We will use the results from BFMAP-SSS.

```
./bfmap --compute_grm 1 --binary_genotype_file geno --snp_info_file all.snp_info.csv --output grm1 --num_threads 10
./bfmap --varcomp --phenotype phen.csv --trait milk --binary_grm_file grm1 --output milk --num_threads 10
./bfmap --sss --phenotype phen.csv --trait milk --snp_info_file topQTL.snp_info.csv --snp_weight weight --binary_genotype_file geno --binary_grm grm1 --heritability 0.3 --output milk.topQTL --num_threads 10
Rscript gene_ppc.R --chr 1 --gtf Sus_scrofa.Sscrofa11.1.111.gtf --m example.topQTL.model.csv --p example.topQTL.pip.csv --o example
```

