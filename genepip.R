library(data.table)
library(dplyr)

# Define a unified function for calculating gene PIP from either BFMAP-SSS or FINEMAP-adj results
calcGenePIP <- function(genes, pip, model, method = c("BFMAP-SSS", "FINEMAP-adj"), chr = NULL) {
  # Parameters:
  # genes: a data.table or data.frame containing gene annotations,
  #        with columns: chr, source, type, start, end, score, phase, attributes.
  # pip: a data.table or data.frame containing SNP PIP data
  #      - For BFMAP-SSS: requires columns Pos, SNPname, Chr
  #      - For FINEMAP-adj: requires columns position, rsid, chromosome
  # model: a data.table or data.frame containing the model data
  #      - For BFMAP-SSS: first (ncol-2) columns are SNP names, last column contains PIP values
  #      - For FINEMAP-adj: requires columns config (comma-separated SNP names) and prob (probability)
  # method: specify which method's output format is being used ("BFMAP-SSS" or "FINEMAP-adj")
  # chr: (optional) chromosome value. If not provided, inferred from pip data
  
  # Validate method parameter
  method <- match.arg(method)
  
  # Convert to data.table if not already
  if (!is.data.table(genes)) genes <- as.data.table(genes)
  if (!is.data.table(pip)) pip <- as.data.table(pip)
  if (!is.data.table(model)) model <- as.data.table(model)
  
  # Set up method-specific variables
  if (method == "BFMAP-SSS") {
    # Get chromosome if not provided
    if (is.null(chr)) chr <- pip$Chr[1]
    
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
                     chr == chr &
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
  output <- output[order(-output$PIP), ]
  
  # Return the output
  return(output)
}