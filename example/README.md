# Example: American Duroc data

Below is a complete workflow example using American Duroc pig data. While fine-mapping is typically performed with sequence data, we use SNP chip data in this example for demonstration purposes.

Refer to https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0218263 for the details.

Download at https://figshare.com/articles/dataset/Porcine_50K_SNP_genotypes_and_phenotypes_of_American_and_Canadian_Duroc_pig_populations/8019551

## GWAS

Constrct the GRM with [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#GREML).

```bash
gcta64 --make-grm  --bfile American_Duroc_pigs_genotypes_qc  --thread-num 10  --out gcta_grm 
gcta64 --mlma --bfile candidate_region --grm gcta_grm --pheno pheno.sim.txt --thread-num 10  --out out.gwa
````

## Heritability estimation and relatedness-adjusted genotype correlation matrix construction
Estimate heritability using [BFMAP](https://github.com/jiang18/bfmap/tree/master).
The heritability of this simulated phenotype is estimated to be 0.542256.

```bash
mph --make_grm --binary_genotype data/American_Duroc_pigs_genotypes_qc --snp_info snp_info.csv --num_threads 10 --out mph_grm

````



## FINEMAP-adj

```R
library(data.table)
gwa_result <- fread("out.gwa,mlma", head =T)
z <- gwa_result[, .(SNP, Chr, bp, A1, A2, Freq, b, se)]
z[Freq > 0.5, Freq := 1 - Freq]

colnames(z) = c("rsid", "chromosome","position","allele1","allele2","maf", "beta","se")
fwrite(z, "data.finemap.z", sep = " ")
```

```bash
finemap --sss --prior-std 0.1 --in-files data --dataset 1
```

## SuSiE-adj

```R
library(susieR)
gwa_result <- read.table("out.gwa.mlma",head=T)
y <- read.csv("pheno.1.csv")
n_eff=960
betahat <- gwa_result[,'b']
sebetahat <- gwa_result[,'se']
fitted_rss1 <- susie_rss(bhat = betahat, shat = sebetahat, n = n_eff, R = R_adj, var_y = var(y[,2]), L = 5,
                         estimate_residual_variance = TRUE)
print(fitted_rss1$converged)
out <- data.frame(SNP = gwa_result[,2], pip = fitted_rss1$pip)
write.table(out,"out.susieadj.pip",quote=F,row.names=F,sep=",")
```

## BFMAP
```bash
bfmap --compute_grm 2 --binary_genotype_file American_Duroc_pigs_genotypes_qc --snp_info_file snp_info.csv --output bfmap_grm --num_threads 10
bfmap --sss --phenotype pheno.sim.csv --trait pheno --snp_info_file snp_info.csv --binary_genotype_file candidate_region --binary_grm bfmap_grm --heritability 0.524424 --output sss --num_threads 10
```

## Gene PIP calculation
Download gene annotation gtf file from ensembl dataset (https://ftp.ensembl.org/pub/release-112/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.112.gtf.gz).

```R
source("finemap_functs.R")
library(data.table)
model <- fread("sss.model.csv",head=T)
genes <- fread("Sus_scrofa.Sscrofa11.1.113.gtf")
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )
genes <- genes[type == "gene"]
snp_pip <- fread("sss.pip.csv",head=T)
result <- calc_gene_pip(genes, snp_pip, model, "BFMAP-SSS", chr = 1)
```
