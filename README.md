# flexDE

**flexDE** is a comprehensive R package designed for performing differential expression analysis on single-cell RNA sequencing data. The package includes various methods such as the *Wilcoxon rank-sum test*, *Wilcoxon rank-sum clustered test*, *Generalized Linear Mixed Models (GLMM)*, *NEBULA*, *Pseudo-bulk edgeR*, and *Pseudo-bulk Limma*. Each function offers options and flexibility for model specifications, making it adaptable to diverse datasets and experimental conditions.

## Installation

You can install the development version from GitHub with:

    # If you don't have the devtools package installed, you need to install it first 
    install.packages("devtools")

    # Load the devtools package
    library(devtools)

    # Install the development version from GitHub
    install_github("EliLillyCo/flexDE") 

## Usage

The package includes two key data processing functions: 'Gene_filter' and 'Data.Pseudobulk'. Some DE functions internally call the 'Gene_filter' function to filter out genes based on specified criteria such as log-fold change and minimum detection percentage. The 'Data.Pseudobulk' function performs the pseudo-bulk procedure, which is required before applying Limma and edgeR methods on single-cell data.

### Requirements for the Seurat Object

Ensure your Seurat object meets the following criteria before using **flexDE** functions:

-   The assay name should be 'RNA'.

-   The assay data should include both @counts (raw data matrix) and @data (normalized data matrix).

-   Metadata should include a "group" variable (e.g., control vs. disease or control vs. treatment).

-   Metadata should include a "subject" variable if you are using models that consider the hierarchical structure of the data.
