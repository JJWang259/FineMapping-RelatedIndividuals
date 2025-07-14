# Example: American Duroc pig data
Below is a workflow example using American Duroc pig data. While fine-mapping is typically performed with sequence data, we use SNP chip data in this example for demonstration.

## Data
The example dataset is provided as [`data.zip`](./data.zip) in the current directory, with original genotype data downloaded from [Zhuang et al. (2019)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0218263).

The following files are included in the zip file:
- American_Duroc_pigs_genotypes_qc.bed (.bim/.fam)
- simulated_pheno.csv

The phenotypic values in `simulated_pheno.csv` are simulated:
- Heritability (*h*²) = 0.5
- Two causal variants: `WU_10.2_1_29501954` (chr1:26247570) and `ALGA0001958` (chr1:27032086)
- Total proportion of variance explained by the two causal variants = 0.04
- Causal gene: *ARFGEF3*

## Tools
- **Genotype data manipulation**: [PLINK 1.9](https://www.cog-genomics.org/plink/)
- **GWAS**: [SLEMM](https://github.com/jiang18/slemm/) or [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#MLMA)
- **Fine-mapping**:
  - **Individual-level data**: [BFMAP](https://github.com/jiang18/bfmap)
  - **Summary statistics**: FINEMAP-adj and SuSiE-adj
    1. **LD matrix adjustment**: [MPH](https://jiang18.github.io/mph/) and [LD Adjuster](../ld_adjuster/)
    2. [FINEMAP](http://www.christianbenner.com/) and/or [susieR](https://stephenslab.github.io/susieR/)

## GWAS

```bash
# Download and unzip data.zip.
# Set the data folder as the working directory.

# Create a file named snp_info.csv listing the SNPs to be included in the GRM.
echo "SNP" > snp_info.csv && awk '{print $2}' American_Duroc_pigs_genotypes_qc.bim >> snp_info.csv

# Run GWAS
slemm --lmm --phenotype_file simulated_pheno.csv --trait trait1 --bfile American_Duroc_pigs_genotypes_qc --snp_info_file snp_info.csv --out trait1 --num_threads 10 --num_qf 100
OMP_NUM_THREADS=10 slemm_gwa --pfile American_Duroc_pigs_genotypes_qc --slemm trait1 --out trait1.chr1.txt --chr 1

# Note that the command above performs association tests for only SNPs on chromosome 1.
````

## Fine-mapping
### LD matrix adjustment
Computing relatedness-adjusted LD matrix requires a GRM and a heritability estimate obtained using MPH:
- [GRM construction](https://jiang18.github.io/mph/options/#making-a-grm-from-snps)
- [Heritability estimation](https://jiang18.github.io/mph/options/#remlminque)

```bash
# Construct GRM, reusing snp_info.csv used by SLEMM
mph --make_grm --binary_genotype American_Duroc_pigs_genotypes_qc --snp_info snp_info.csv --num_threads 10 --out mph_chip

# Estimate variance components
echo "mph_chip 1" > grm_list.txt
mph --reml --grm_list grm_list.txt --phenotype simulated_pheno.csv --trait trait1 --num_threads 10 --out trait1
````

The heritability estimate is shown in `trait1.mq.vc.csv`:
```csv
trait_x,trait_y,vc_name,m,var,seV,pve,seP,enrichment,seE,mph_grm,mph_grm,err
pheno,pheno,mph_grm,38646,0.531154,0.0364273,0.525258,0.0214742,1,NA,NA,0.00132695,-0.000221669
pheno,pheno,err,NA,0.480072,0.0149826,0.474742,0.0214742,NA,NA,NA,-0.000221669,0.000224478
```
The estimated *h*² for this simulated trait is `0.525258`.

Based on the association analysis, we defined chr1:26,000,000-30,000,000 as the candidate region for fine-mapping.

```bash
plink --bfile American_Duroc_pigs_genotypes_qc --chr 1 --from-mb 26 --to-mb 30 --recode A --out candidate_region 
ld_adjuster --raw candidate_region.raw --grm chip --h2 0.525258 --out ld_adjusted --threads 10
````
The `ld_adjuster` command above generates three output files for fine-mapping using summary statistics:
- `ld_adjusted.summary` → effective sample size (`750.727` in this example)
- `ld_adjusted.ld` → LD matrix 
- `ld_adjusted.snpids` → SNP identifiers with counted alleles (e.g., `rs1234567_T`)

### FINEMAP-adj
Prepare summary statistics for FINEMAP:

```sh
plink --bfile American_Duroc_pigs_genotypes_qc --freq --out pigs
```

```R
library(data.table)

ld_adjusted_prefix = "ld_adjusted"
out_prefix = "finemap_adj"
gwa_file = "trait1.chr1.txt"
frq_file = "pigs.frq"
n_eff = 751        # Must be an integer (required by FINEMAP)

gwa_result <- fread(gwa_file, head =T)
maf <- fread(frq_file, head=T)
gwa_result <- gwa_result[maf[, .(SNP, MAF)], on = "SNP", nomatch = 0]
z <- gwa_result[, .(SNP, CHR, BP, A1, A2, MAF, BETA, SE)]

snpids <- fread(paste0(ld_adjusted_prefix,".snpids"), header = FALSE)
snpids[, `:=`( SNP = sub("_[^_]*$", "", V1), counted_allele = sub(".*_", "", V1))]

z <- snpids[z, on = "SNP", nomatch = 0]
z[A2 != counted_allele, `:=`(BETA = -BETA)]
z <- z[, .(SNP, CHR, BP, A1, A2, MAF, BETA, SE)]
colnames(z) = c("rsid", "chromosome","position","allele1","allele2","maf", "beta","se")
fwrite(z, paste0(out_prefix, ".z"), sep = " ")

# Create FINEMAP master file
finemap_master <- data.table(
  z = paste0(out_prefix, ".z"),
  ld = paste0(ld_adjusted_prefix, ".ld"),
  snp = paste0(out_prefix, ".snp"),
  config = paste0(out_prefix, ".config"), 
  cred = paste0(out_prefix, ".cred"),
  log = paste0(out_prefix, ".log"),
  n_samples = n_eff  # Use the effective sample size from ld_adjusted.summary
)

# Write master file
fwrite(finemap_master, paste0(out_prefix, ".data"), sep = ";")
```
Run FINEMAP:
```bash
finemap --sss --prior-std 0.1 --in-files finemap_adj.data --dataset 1
```

### SuSiE-adj

```R
library(susieR)
library(data.table)

ld_adjusted_prefix = "ld_adjusted"
out_prefix = "susie_adj"
gwa_file = "trait1.chr1.txt"
var_y = 1          # Phenotypic variance, needed by susieR
n_eff = 751        

snpids <- fread(paste0(ld_adjusted_prefix,".snpids"), header = FALSE)
snpids[, `:=`( SNP = sub("_[^_]*$", "", V1), counted_allele = sub(".*_", "", V1))]

gwa_result <- fread(gwa_file, head =T)
gwa_result <- snpids[gwa_result, on = "SNP", nomatch = 0]
gwa_result[A2 != counted_allele, `:=`(BETA = -BETA)]

R_adj = fread(paste0(ld_adjusted_prefix, ".ld"))
fit <- susie_rss(bhat = gwa_result[,BETA], shat = gwa_result[,SE], n = n_eff, R = R_adj, var_y = var_y, L = 5, estimate_residual_variance = TRUE)
print(fit$converged)
out <- data.frame(SNP = gwa_result[,SNP], pip = fit$pip)
out <- out[order(-out$pip),]
write.csv(out, paste0(out_prefix, ".pip.csv"),quote=F,row.names=F)
```

### BFMAP
```bash
# Construct GRM with BFMAP, reusing snp_info.csv used by SLEMM and MPH
# Note that BFMAP is currently not compatible with MPH's GRM format.
bfmap --compute_grm 2 --binary_genotype_file American_Duroc_pigs_genotypes_qc --snp_info_file snp_info.csv --output bfmap --num_threads 10

# Extract SNPs in the candidate region
plink --bfile American_Duroc_pigs_genotypes_qc --chr 1 --from-mb 26 --to-mb 30 --make-bed --out candidate_region
echo "SNP" > candidate_snp_info.csv && awk '{print $2}' candidate_region.bim >> candidate_snp_info.csv

# Perform shotgun stochastic search with BFMAP
bfmap --sss --phenotype simulated_pheno.csv --trait trait1 --snp_info_file candidate_snp_info.csv --binary_genotype_file candidate_region --binary_grm bfmap --heritability 0.525258 --output sss --num_threads 10
```

## Gene PIP

### Prerequisites
- [`calc_gene_pip.R`](../calc_gene_pip.R)
- Output files from BFMAP-SSS or FINEMAP-adj
- [Ensembl GTF file](https://www.ensembl.org/info/website/upload/gff.html?redirect=no) for gene annotations

### Download gene annotation
```bash
wget https://ftp.ensembl.org/pub/release-113/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.113.gtf.gz
gunzip Sus_scrofa.Sscrofa11.1.113.gtf.gz
````

### Calculate gene PIPs
```R
source("calc_gene_pip.R")
library(data.table)
gtf <- fread("Sus_scrofa.Sscrofa11.1.113.gtf", sep="\t", head = FALSE)
setnames(gtf, names(gtf), c("seqname","source","feature","start","end","score","strand","frame","attribute") )
gtf <- gtf[feature == "gene"]

# BFMAP-SSS gene PIP calculation
sss_pip <- fread("sss.pip.csv")
sss_model <- fread("sss.model.csv")
sss_gene_pip <- calc_gene_pip(gtf, sss_pip, sss_model)

# FINEMAP gene PIP calculation  
finemap_pip <- fread("finemap_adj.snp")
finemap_model <- fread("finemap_adj.config")
finemap_gene_pip <- calc_gene_pip(gtf, finemap_pip, finemap_model)
```
### Output format
The function returns a data frame with the following columns:
- `Chr`: Chromosome
- `Start`: Gene start position
- `End`: Gene end position  
- `PIP`: Gene-level posterior inclusion probability
- `Attribute`: Attribute from GTF
