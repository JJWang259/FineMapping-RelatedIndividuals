###########################################
#                                         #
#  Utility functions for finemap-related  #
#                                         #
###########################################

library(data.table)
library(dplyr)

#---------------------------------------
# Relatedness-adjusted LD matrix generation
#---------------------------------------

#' Generate Relatedness-adjusted genotype correlation matrix
#' 
#' @param raw Raw genotype matrix
#' @param h2 Heritability estimate
#' @param G Genomic relationship matrix
#' @return List containing adjusted LD matrix and effective sample size
#' @export
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
  
  n_eff <- round(mean(diag(Rstar)))
  message(paste0("Mean effective sample size is: ", n_eff))
  
  Ds <- sqrt(diag(Rstar))
  D_sol <- diag(1 / Ds, nrow = length(Ds), ncol = length(Ds))
  R_adj <- D_sol %*% Rstar %*% D_sol
  
  return(list(ld_matrix = R_adj, n_eff = n_eff))
}



#---------------------------------------
# Gene PIP Calculation
#---------------------------------------

#' Calculate Gene-level Posterior Inclusion Probability
#' 
#' @param genes A data.table or data.frame containing gene annotations
#' @param pip A data.table containing SNP PIP output
#' @param model A data.table containing model configurations
#' @param method Which method's output format is being used ("BFMAP-SSS" or "FINEMAP-adj")
#' @param chr Chromosome value (optional)
#' @return Data frame with gene-level PIP values
#' @export
calc_gene_pip <- function(genes, pip, model, method = c("BFMAP-SSS", "FINEMAP-adj"), chr = NULL) {
  # Validate method parameter
  method <- match.arg(method)
  genes$chr  = as.character(genes$chr)
  target_chr <- as.character(chr)
  # Convert to data.table if not already
  if (!is.data.table(genes)) genes <- as.data.table(genes)
  if (!is.data.table(pip)) pip <- as.data.table(pip)
  if (!is.data.table(model)) model <- as.data.table(model)
  # Set up method-specific variables
  if (method == "BFMAP-SSS") {
    # Get chromosome if not provided
    if (is.null(chr)) chr <- as.character(pip$Chr[1])
    
    # Define SNP range
    st <- min(pip$Pos)
    ed <- max(pip$Pos)
    
    # Column mappings
    pos_col <- "Pos"
    snp_col <- "SNPname"
  } else { # FINEMAP-adj
    # Get chromosome if not provided
    if (is.null(chr)) chr <- pip$chromosome[1]
    
    # Define SNP range
    st <- min(pip$position)
    ed <- max(pip$position)
    
    # Column mappings
    pos_col <- "position"
    snp_col <- "rsid"
  }
  
  # Subset genes
  gene_sub <- subset(genes,
                     chr == target_chr &
                     ((start >= st & start <= ed) | (end >= st & end <= ed)))
  
  # Create an output data frame
  output <- data.frame(matrix(ncol = 5, nrow = nrow(gene_sub)))
  colnames(output) <- c("Chr", "Start", "End", "PIP", "Attributes")
  
  # Loop over each gene in gene_sub
  for (i in seq_len(nrow(gene_sub))) {
    # Select SNPs in the pip file that fall within the gene's boundaries
    if (method == "BFMAP-SSS") {
      snp_list <- subset(pip, get(pos_col) >= gene_sub$start[i] & get(pos_col) <= gene_sub$end[i])
      gene_snps <- snp_list[[snp_col]]
      
      # Calculate gene PIP
      nSNPcols <- ncol(model) - 2
      if (nSNPcols > 0 && length(gene_snps) > 0) {
        in_snp_mat <- sapply(1:nSNPcols, function(j) model[[j]] %in% gene_snps)
        in_snp_result <- rowSums(in_snp_mat) > 0
        model_sub <- model[in_snp_result, ]
        gene_pip <- if (nrow(model_sub) > 0) sum(model_sub[[ncol(model_sub)]]) else 0
      } else {
        gene_pip <- 0
      }
    } else { # FINEMAP-adj
      snp_list <- subset(pip, get(pos_col) >= gene_sub$start[i] & get(pos_col) <= gene_sub$end[i])
      gene_snps <- unique(snp_list[[snp_col]])
      
      # Calculate gene PIP
      if (length(gene_snps) == 0) {
        gene_pip <- 0
      } else {
        include_flag <- sapply(model$config, function(x) {
          conf_snps <- unlist(strsplit(x, split = ","))
          any(conf_snps %in% gene_snps)
        })
        gene_pip <- sum(model$prob[include_flag])
      }
    }
    
    # Assign gene information and calculated gene PIP to the output
    output[i, c("Chr", "Start", "End", "Attributes")] <- gene_sub[i, c("chr", "start", "end", "attributes")]
    output[i, "PIP"] <- gene_pip
  }
  
  # Order the output by descending gene PIP
  output <- subset(output, PIP >0)
  output <- output[order(-output$PIP), ]
  
  # Return the output
  return(output)
}