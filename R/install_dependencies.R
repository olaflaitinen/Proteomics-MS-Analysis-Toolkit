# R/install_dependencies.R

# --- CRAN Packages ---
cran_packages <- c(
  "tidyverse", 
  "ggplot2", 
  "ggrepel", 
  "pheatmap",
  "remotes"
)

# Install CRAN packages if they are not already installed
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# --- Bioconductor Packages ---
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "limma",
  "clusterProfiler",
  "org.Hs.eg.db"  # Example organism annotation package
)

# Install Bioconductor packages
BiocManager::install(bioc_packages, update = FALSE)

cat("\nAll required R packages have been installed.\n")
