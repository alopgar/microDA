[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![Last version](https://img.shields.io/github/tag/alopgar/microDA.svg)](https://img.shields.io/github/tag/alopgar/microDA.svg)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/microDA)](https://cran.r-project.org/package=microDA)

# microDA
R package with expanded applications of classic R functions for microbiome datasets

This repository contains miscellaneous functions for microbiome data analysis using utilities from other packages.
Still in **ALPHA** state. 

# Getting started
## Installation
    install.packages("devtools")
    devtools::install_github("alopgar/microDA")

## Dependencies
- **tidyverse**: `install.packages("tidyverse")`
- **microbiome**: `BiocManager::install("microbiome")` or `devtools::install_github("microbiome/microbiome")`
- **phyloseq**: `source('http://bioconductor.org/biocLite.R'); biocLite('phyloseq')`
- **vegan**: `install.packages("vegan")`
- **nortest**: `install.packages("nortest")`
- **DESeq2**: `BiocManager::install("DESeq2")`
- **ALDEx2**: `BiocManager::install("ALDEx2")`
- **sunburstR**: `install.packages("sunburstR")`
- **ellipse**: `install.packages("ellipse")`
- **VennDiagram**: `install.packages("VennDiagram")`
- **ggrepel**: `install.packages("ggrepel")`
- **ggordiplots**: `remotes::install_github("jfq3/ggordiplots")`
- **colortools**: `devtools::install_github('gastonstat/colortools')`

# Acknowledgements
This repository uses functions from multiple R packages. Please cite R and those R packages when using it.

Some included functions have been adapted from other public repositories:
- get_max_taxonomic_rank: From [metagMisc](https://github.com/vmikk/metagMisc)
- phyloseq_to_df: From [metagMisc](https://github.com/vmikk/metagMisc)
- CIplot_biv: From [easyCODA](https://github.com/cran/easyCODA)
