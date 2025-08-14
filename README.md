# Proteomics MS Analysis Toolkit

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
![Languages](https://img.shields.io/badge/languages-Python%20%7C%20R-orange.svg)
[![Python CI](https://github.com/olaflaitinen/Proteomics-MS-Analysis-Toolkit/actions/workflows/python-ci.yml/badge.svg)](https://github.com/olaflaitinen/Proteomics-MS-Analysis-Toolkit/actions/workflows/python-ci.yml)

A collection of Python and R scripts for the analysis of quantitative proteomics data from mass spectrometry, specifically tailored for MaxQuant `proteinGroups.txt` files.

This toolkit implements a standard workflow including data cleaning, quality control, differential abundance analysis, and functional enrichment analysis.

![Example Volcano Plot](assets/volcano_plot_example.png)

### Overview

The goal of this toolkit is to provide a reproducible and easy-to-use set of tools for researchers to analyze their proteomics data. The workflow is demonstrated using both Python (in a Jupyter Notebook) and R (in an R Markdown document).

### Features

-   **Data Preprocessing**: Filters out contaminants, reverse hits, and proteins only identified by site. Handles missing values and performs log2 transformation.
-   **Normalization**: Implements median normalization to adjust for technical variability.
-   **Differential Abundance Analysis**:
    -   **Python**: Welch's t-test using `scipy`.
    -   **R**: Linear Models for Microarray Data (`limma`) for robust statistical analysis.
-   **Functional Enrichment Analysis**:
    -   **Python**: Gene Set Enrichment Analysis (GSEA) using `gseapy`.
    -   **R**: Over-Representation Analysis (ORA) using `clusterProfiler`.
-   **Visualization**: Generates publication-quality plots, including Volcano Plots, PCA plots, and Heatmaps.

### Project Structure

```
├── R/                  # R implementation with R Markdown walkthrough
├── python/             # Python implementation with Jupyter Notebook
├── data/               # Example MaxQuant output and metadata
├── assets/             # Images for documentation
└── README.md           # This file
```

### Installation

**1. Clone the repository:**

```bash
git clone https://github.com/olaflaitinen/Proteomics-MS-Analysis-Toolkit.git
cd Proteomics-MS-Analysis-Toolkit
```

**2. Set up the Python environment:**

It is recommended to use a virtual environment.

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r python/requirements.txt
```

**3. Set up the R environment:**

Run the dependency installation script within an R session:

```R
# In your R console
source("R/install_dependencies.R")
```

### Quick Start

The easiest way to use this toolkit is to follow the provided example walkthroughs with the sample data.

-   **For Python users**:
    Open and run the Jupyter Notebook:
    ```bash
    cd python
    jupyter-lab analysis_notebook.ipynb
    ```

-   **For R users**:
    Open `R/analysis_walkthrough.Rmd` in RStudio and click "Knit" to generate a full HTML report, or run the code chunks interactively.

### Citation

The methods implemented in this toolkit are based on standard practices in the field and are described in publications such as:

> Jensen, Lars, Sofie Rasmussen and O. Yunus L. Imanov. (2025). A Scalable and Reproducible Bioinformatics Pipeline for Differential Analysis of Mass Spectrometry-based Proteomics Data. *Journal of Proteome Research 24(2): 112-125*.

### License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
