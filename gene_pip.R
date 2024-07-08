library(argparse)
parser <- ArgumentParser(description = "PPC calculation")
parser$add_argument("--chr", help = "specify chromosome", required=TRUE)
parser$add_argument("--gtf", help = "input gtf file", required=TRUE)
parser$add_argument("--m", help = "input model file", required=TRUE)
parser$add_argument("--p", help = "input pip file", required=TRUE)
parser$add_argument("--o", help = "specify the prefix of output file", required=TRUE)

args <- parser$parse_args()

chr = args$chr
output_file <- args$o
gtf_file <- args$gtf
pip_file <- args$p
model_file <- args$m
library(data.table)
library(dplyr)
library(purrr)
genes <- fread(gtf_file, skip = 5)

setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes") )



pip <- fread(pip_file)
model <- fread(model_file)

#define range
st <- min(pip$Pos)
ed <- max(pip$Pos)
chr_ref = pip$Chr[1]
gene_sub <- subset(genes, type == "gene" & chr == chr_ref & ((start >= st & start <= ed)|(end >= st & end <= ed)))
output <- data.frame(matrix(ncol=6,nrow=(nrow(gene_sub))))
colnames(output) <- c("Chr", "Start", "End", "Strand","PPC","Attributes")



for (i in 1:nrow(gene_sub)){
  snp_list <- subset(pip, Pos >= gene_sub$start[i] & Pos <= gene_sub$end[i])
  in_snp_list <- data.frame(matrix(nrow =nrow(model),ncol =1))
  for (j in 1:(ncol(model)-2)){
    in_snp_list[,j] <- model[, j, with =FALSE][[1]] %in% snp_list$SNPname
  }
  in_snp_list$result <- rowSums(in_snp_list) > 0

  model_sub <- model[in_snp_list$result, ]
  n = ncol(model_sub)

  output[i,c(1,2,3,4,6)] = gene_sub[i,c("chr", "start", "end", "strand","attributes")]
  output[i,5] = sum(model_sub[,..n])
}
output <- output[order(-output$PPC),]

output_file <- paste0(output_file,".geneppc.csv")
write.csv(output, output_file, quote =F, row.names =F)
