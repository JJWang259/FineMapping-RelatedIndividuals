# Example: American Duroc data
Refer to https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0218263 for the details.

Download at https://figshare.com/articles/dataset/Porcine_50K_SNP_genotypes_and_phenotypes_of_American_and_Canadian_Duroc_pig_populations/8019551

## GRM construction

Constrct the GRM with [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#GREML).
```bash
gcta64 --make-grm  --bfile American_Duroc_pigs_genotypes  --thread-num 20  --out us_duroc 
````

## Heritability estimation

Estimate heritability using [BFMAP](https://github.com/jiang18/bfmap/tree/master).

```bash
bfmap --compute_grm 1 --binary_genotype_file ./American_Duroc_pigs_genotypes --snp_info_file ./all.snp_info.csv --output grm1 --num_threads 10
bfmap --varcomp --phenotype ./pheno.1.csv --trait pheno --binary_grm_file ./grm1 --output pheno1
```
The heritability of this simulated phenotype is estimated to be 0.508427.

## Relatedness-adjusted genotype correlation matrix

```R
ld_adj <- function(raw, h2, G) {
  eta <- h2 / (1 - h2)
  diag(G) <- eta * diag(G) + 1
  L <- t(chol(G))
  X <- as.matrix(raw)
  
  # Efficient scaling and centering
  X <- scale(X, center = TRUE, scale = TRUE)
  
  # Efficient solve
  Xstar <- solve(L, X)
  Rstar <- crossprod(Xstar)  # More efficient way to compute t(Xstar) %*% Xstar
  
  n_eff <- mean(diag(Rstar))
  print(paste0("Mean effective sample size is: ", n_eff))
  
  Ds <- sqrt(diag(Rstar))
  D_sol <- diag(1 / Ds, nrow = length(Ds), ncol = length(Ds))
  R_adj <- D_sol %*% Rstar %*% D_sol
  
  return(list(ld_matrix = R_adj, n_eff = n_eff))
}

# read grm
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  closeAllConnections()
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}


raw = read.table("add.raw",head=T)
raw = raw[,-(1:6)]
bin = ReadGRMBin( "us_duroc")
np = length(bin$diag)
G = matrix(0, nrow=np, ncol=np)
G[upper.tri(G)] = bin$off
G = G + t(G)
diag(G) = bin$diag
h2 = 0.508427
result <- ld_adj(raw, h2, G)
R_adj <- result$ld_matrix
n_eff <- result$n_eff
write.table(R_adj, "adj.ld", quote = FALSE, row.names = FALSE, col.names = FALSE)
```
The mean effective sample size is 960.

## SuSiE-Adj

```R
library(susieR)
gwa_result <- read.table("out.gwa.1.txt.mlma",head=T)
y <- read.csv("pheno.1.csv")
n_eff=960
betahat <- gwa_result[,'b']
sebetahat <- gwa_result[,'se']
fitted_rss1 <- susie_rss(bhat = betahat, shat = sebetahat, n = n_eff, R = R_adj, var_y = var(y[,2]), L = 5,
                         estimate_residual_variance = TRUE)
print(fitted_rss1$converged)
out <- cbind(gwa_result[,2],fitted_rss1$pip)
colnames(out) = c("SNP","pip")
write.table(out,"pheno1.susieadj.pip",quote=F,row.names=F,sep=",")
```

## FINEMAP-Adj

```bash
finemap --sss --prior-std 0.1 --in-files data --dataset 1
```

## Gene PIP calculation
Download gene annotation gtf file from ensembl dataset (https://ftp.ensembl.org/pub/release-112/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.112.gtf.gz).

```bash
bfmap --sss --phenotype pheno.1.csv --trait pheno --snp_info_file snp_info.csv --binary_genotype_file geno_region --binary_grm grm1 --heritability 0.508427 --output pheno1 --num_threads 10
Rscript gene_ppc.R --chr 6 --gtf Sus_scrofa.Sscrofa11.1.112.gtf --m pheno1.topQTL.model.csv --p pheno1.topQTL.pip.csv --o pheno1
```
