# Example: American Duroc pig data
Below is a workflow example using American Duroc pig data. While fine-mapping is typically performed with sequence data, we use SNP chip data in this example for demonstration.

## Data
The example data is provided as [`data.zip`](./data.zip) in the current directory, with original genotype data downloaded from [Zhuang et al. (2019)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0218263).

**Simulated phenotype**: For this demonstration, we simulated a phenotype with:
- Heritability (h²) = 0.5
- Two causal variants: `WU_10.2_1_29501954` and `ALGA0001958`
- Total proportion of variance explained by causal variants = 0.04
- Causal gene: *ARFGEF3*

The following files are included in the zip file:
- American_Duroc_pigs_genotypes_qc.bed (.bim/.fam)
- pheno.sim.txt

## Tools
- **GWAS**: [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#MLMA)
- **LD matrix adjustment**: [MPH](https://jiang18.github.io/mph/) and [LD Adjuster](../ld_adjuster/)
- **Genotype file conversion**: [PLINK 1.9](https://www.cog-genomics.org/plink/)
- **Fine-mapping**:
  - **BFMAP** (https://github.com/jiang18/bfmap)
  - **FINEMAP** (http://www.christianbenner.com/)
  - **susieR** (install in R with `install.packages("susieR")`)


## GWAS
For this demonstration, we perform associations of SNPs in a pre-selected region (Chr1:20,000,000-30,000,000) rather than a genome-wide analysis.

```bash
# Set the data folder as the working directory.

# Extract candidate region
plink --bfile American_Duroc_pigs_genotypes_qc --chr 1 --from-mb 20 --to-mb 30 --make-bed --out candidate_region
# Construct GRM using full dataset
gcta64 --make-grm  --bfile American_Duroc_pigs_genotypes_qc  --thread-num 10  --out gcta_grm
# Run GWAS
gcta64 --mlma --bfile candidate_region --grm gcta_grm --pheno pheno.sim.txt --thread-num 10  --out out.gwa
````

## Relatedness-adjusted LD matrix

[GRM construction](https://jiang18.github.io/mph/options/#making-a-grm-from-snps) and [heritability estimation](https://jiang18.github.io/mph/options/#remlminque) using MPH.
```bash
echo "SNP" > snp_info.csv && awk '{print $2}' American_Duroc_pigs_genotypes_qc.bim >> snp_info.csv
mph --make_grm --binary_genotype American_Duroc_pigs_genotypes_qc --snp_info snp_info.csv --num_threads 10 --out mph_grm
echo "mph_grm 1" > grm_list.txt
awk 'NR==1 {print "ID,pheno"} NR>1 {print $2","$3}' pheno.sim.txt > pheno.sim.csv
mph --reml --grm_list grm_list.txt --phenotype pheno.sim.csv --trait pheno --num_threads 10 --out mph_h2
````

The heritability estimate is found in `mph_h2.mq.vc.csv`:
```
trait_x,trait_y,vc_name,m,var,seV,pve,seP,enrichment,seE,mph_grm,mph_grm,err
pheno,pheno,mph_grm,38646,0.531154,0.0364273,0.525258,0.0214742,1,NA,NA,0.00132695,-0.000221669
pheno,pheno,err,NA,0.480072,0.0149826,0.474742,0.0214742,NA,NA,NA,-0.000221669,0.000224478
```
The estimated heritability for this simulated phenotype is **h² = 0.525258**


### Relatedness-adjusted LD matrix
```bash
plink --bfile American_Duroc_pigs_genotypes_qc --extract candidate_list.csv --recode A --out candidate_region 
ld_adjuster --raw candidate_region.raw --grm mph_grm --h2 0.525258 --out adj --threads 10
````
LD Adjuster generates three outputs that are essential for fine-mapping:
- `adj.summary` → Contains **effective sample size** (804 in this example) for fine-mapping
- `adj.ld` → **LD correlation matrix** input for summary-statistics based fine-mapping
- `adj.snpids` → **SNP identifiers** with counted alleles (e.g., `rs1234567_T`)


## FINEMAP-adj
Prepare summary statistics for FINEMAP.
```R
library(data.table)
gwa_result <- fread("out.gwa.mlma", head =T)
z <- gwa_result[, .(SNP, Chr, bp, A1, A2, Freq, b, se)]
snpids <- fread("adj.snpids", header = FALSE)
snpids[, `:=`( SNP = sub("_[A-Z]$", "", V1), counted_allele = sub(".*_", "", V1))]
z <- snpids[z, on = "SNP", nomatch = 0]
z[A1 != counted_allele, `:=`(A1 = A2, A2 = A1, b = -b)]
z[Freq > 0.5, Freq := 1 - Freq]
z <- z[, .(SNP, Chr, bp, A1, A2, Freq, b, se)]
colnames(z) = c("rsid", "chromosome","position","allele1","allele2","maf", "beta","se")
fwrite(z, "data.finemap.z", sep = " ")

# Create FINEMAP master file
finemap_master <- data.table(
  z = "data.finemap.z",
  ld = "adj.ld",
  snp = "out.finemap.snp",
  config = "out.finemap.config", 
  cred = "out.finemap.cred",
  log = "out.finemap.log",
  n_samples = 804  # Use the effective sample size from adj.summary
)

# Write master file
fwrite(finemap_master, "data", sep = ";")
```
Run FINEMAP:
```bash
finemap --sss --prior-std 0.1 --in-files data --dataset 1
```

## SuSiE-adj

```R
library(susieR)
library(data.table)
gwa_result <- fread("out.gwa.mlma", head =T)
z <- gwa_result[, .(SNP, Chr, bp, A1, A2, Freq, b, se)]
snpids <- fread("adj.snpids", header = FALSE)
snpids[, `:=`( SNP = sub("_[A-Z]$", "", V1), counted_allele = sub(".*_", "", V1))]
z <- snpids[z, on = "SNP", nomatch = 0]
z[A1 != counted_allele, `:=`(A1 = A2, A2 = A1, b = -b)]
y <- fread("pheno.sim.csv")
R_adj = fread("adj.ld")
n_eff=804
betahat <- z[,b]
sebetahat <- z[,se]
fitted_rss1 <- susie_rss(bhat = betahat, shat = sebetahat, n = n_eff, R = R_adj, var_y = var(y[,2]), L = 5,
                         estimate_residual_variance = TRUE)
print(fitted_rss1$converged)
out <- data.frame(SNP = gwa_result[,2], pip = fitted_rss1$pip)
out <- out[order(-out$pip),]
write.table(out,"out.susieadj.pip",quote=F,row.names=F,sep=",")
```

## BFMAP
```bash
echo "SNP" > snp_info.csv && awk '{print $2}' American_Duroc_pigs_genotypes_qc.bim >> snp_info.csv
bfmap --compute_grm 2 --binary_genotype_file American_Duroc_pigs_genotypes_qc --snp_info_file snp_info.csv --output bfmap_grm --num_threads 10
bfmap --sss --phenotype pheno.sim.csv --trait pheno --snp_info_file snp_info.csv --binary_genotype_file candidate_region --binary_grm bfmap_grm --heritability 0.525258 --output sss --num_threads 10
```

## Gene PIP calculation

### Prerequisites
- `calc_gene_pip.R` containing the `calc_gene_pip` function (provided in this repository)
- Completed BFMAP-SSS and FINEMAP-adj analyses from previous steps
Download gene annotation gtf file from ensembl dataset (https://ftp.ensembl.org/pub/release-113/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.113.gtf.gz).


### Download gene annotation
```bash
wget https://ftp.ensembl.org/pub/release-113/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.113.gtf.gz
gunzip Sus_scrofa.Sscrofa11.1.113.gtf.gz
````
### Calculate gene PIPs
```R
source("calc_gene_pip.R")
library(data.table)
genes <- fread("Sus_scrofa.Sscrofa11.1.113.gtf")
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )
genes <- genes[type == "gene"]

# BFMAP-SSS gene PIP calculation
sss_pip <- fread("sss.pip.csv", head = TRUE)
sss_model <- fread("sss.model.csv", head = FALSE)
sss_genepip <- calc_gene_pip(genes, sss_pip, sss_model)


# FINEMAP gene PIP calculation  
finemap_pip <- fread("out.finemap.snp")
finemap_model <- fread("out.finemap.config", head = TRUE)
finemap_genepip <- calc_gene_pip(genes, finemap_pip, finemap_model)
```
### Output format
The function returns a data frame with columns:
- `Chr`: Chromosome
- `Start`: Gene start position
- `End`: Gene end position  
- `PIP`: Gene-level posterior inclusion probability
- `Attributes`: Gene annotation details from GTF
