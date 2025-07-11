calc_gene_pip <- function(genes, pip, model, method = NULL, ext = 3000) {
  # Convert to data.table if not already
  if (!is.data.table(genes)) genes <- as.data.table(genes)
  if (!is.data.table(pip)) pip <- as.data.table(pip)
  if (!is.data.table(model)) model <- as.data.table(model)
  
  # Auto-detect method if not provided
  if (is.null(method)) {
    pip_cols <- names(pip)
    
    # Check for method-specific columns
    if ("SNPname" %in% pip_cols && "Chr" %in% pip_cols && "Pos" %in% pip_cols) {
      method <- "BFMAP-SSS"
      message("Auto-detected BFMAP-SSS format")
    } else if ("rsid" %in% pip_cols && "chromosome" %in% pip_cols && "position" %in% pip_cols) {
      method <- "FINEMAP-adj"
      message("Auto-detected FINEMAP-adj format")
    } else {
      stop("Could not auto-detect file format. Please specify method = 'BFMAP-SSS' or 'FINEMAP-adj'")
    }
  }
  
  # Validate required columns based on method
  if (method == "BFMAP-SSS") {
    required_cols <- c("SNPname", "Chr", "Pos")
    missing_cols <- setdiff(required_cols, names(pip))
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns for BFMAP-SSS:", paste(missing_cols, collapse = ", ")))
    }
    
    # Set column mappings
    chr_col <- "Chr"
    pos_col <- "Pos"
    snp_col <- "SNPname"
    
    # Convert Chr to character
    pip[[chr_col]] <- as.character(pip[[chr_col]])
    
  } else if (method == "FINEMAP-adj") {
    required_cols <- c("rsid", "chromosome", "position", "prob")
    missing_cols <- setdiff(required_cols, names(pip))
    if (length(missing_cols) > 0) {
      stop(paste("Missing required columns for FINEMAP-adj:", paste(missing_cols, collapse = ", ")))
    }
    
    # Validate model columns
    if (!"config" %in% names(model) || !"prob" %in% names(model)) {
      stop("FINEMAP-adj model file must contain 'config' and 'prob' columns")
    }
    
    # Set column mappings
    chr_col <- "chromosome"
    pos_col <- "position"
    snp_col <- "rsid"
    
    # Convert chromosome to character
    pip[[chr_col]] <- as.character(pip[[chr_col]])
    
  } else {
    stop("Method must be 'BFMAP-SSS' or 'FINEMAP-adj'")
  }
  
  # Ensure genes chr is character
  genes$chr <- as.character(genes$chr)
  
  # Get chromosome list and SNP range
  chrlist <- unique(pip[[chr_col]])
  st <- min(pip[[pos_col]])
  ed <- max(pip[[pos_col]])
  
  # Subset genes
  gene_sub <- subset(genes,
                     chr %in% chrlist &
                     ((start >= st & start <= ed) | (end >= st & end <= ed)))
  
  # Create an output data frame
  output <- data.frame(matrix(ncol = 5, nrow = nrow(gene_sub)))
  colnames(output) <- c("Chr", "Start", "End", "PIP", "Attributes")
  
  # Loop over each gene in gene_sub
  for (i in seq_len(nrow(gene_sub))) {
    # Select SNPs in the pip file that fall within the gene's boundaries
    snp_list <- pip[get(chr_col) == gene_sub$chr[i] & 
                    get(pos_col) >= (gene_sub$start[i] - ext) & 
                    get(pos_col) <= (gene_sub$end[i] + ext)]
    
    if (method == "BFMAP-SSS") {
      gene_snps <- snp_list[[snp_col]]
      
      # Calculate gene PIP
      nSNPcols <- ncol(model) - 2
      if (nSNPcols > 0 && length(gene_snps) > 0) {
        # Create matrix of whether each model SNP is in gene_snps
        in_snp_mat <- sapply(1:nSNPcols, function(j) model[[j]] %in% gene_snps)
        
        # Handle case where sapply returns a vector (single model)
        if (is.vector(in_snp_mat)) {
          in_snp_mat <- matrix(in_snp_mat, nrow = 1)
        }
        
        # Find models containing at least one gene SNP
        in_snp_result <- rowSums(in_snp_mat) > 0
        model_sub <- model[in_snp_result, ]
        
        # Sum probabilities (last column)
        gene_pip <- if (nrow(model_sub) > 0) sum(model_sub[[ncol(model_sub)]]) else 0
      } else {
        gene_pip <- 0
      }
      
    } else { # FINEMAP-adj
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
  output <- subset(output, PIP > 0)
  output <- output[order(-output$PIP), ]
  
  # Return the output
  return(output)
}