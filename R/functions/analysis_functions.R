# R/functions/analysis_functions.R

# Function to run differential expression analysis using limma
run_limma_analysis <- function(processed_data, metadata_file) {
  meta <- read_csv(metadata_file)
  
  # Ensure data and metadata are aligned
  mat <- as.matrix(processed_data[, meta$sample])
  
  # Create design matrix
  condition <- factor(meta$condition)
  design <- model.matrix(~0 + condition)
  colnames(design) <- levels(condition)
  
  # Create contrast matrix
  contrast.matrix <- makeContrasts(treatment - control, levels = design)
  
  # Fit the linear model
  fit <- lmFit(mat, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  # Get the results
  results <- topTable(fit2, number = Inf, sort.by = "P")
  results <- rownames_to_column(results, var = "Gene")
  
  return(results)
}
