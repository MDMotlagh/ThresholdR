# Overview
ThresholdR uses Gaussian Mixture Models concept to identify the noise population in each surface marker across cells in CITE-seq experiments. Once the noise distribution is identified, it calculates the upper threshold of the noise component to separate expressing and non-expressing cells for each surface marker. 
# Citation
# Installation
```r
#1. Install required packages:
required_packages <- c("mixtools", "mclust", "foreach", "doParallel", "ggplot2", "AdaptGauss")

# Check and install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
#2. Install ThresholdR:
remotes::install_guithub('MDMotlagh/ThresholdR')
```
# Instructions
## Contents
| Folder | Description | Details |
| --------------- | --------------- | --------------- |
| data   | PBMC10K and bmcite datasets   | PBMC10K immune profiling dataset acquired from 10x genomics website. bmcite is available through SeuratData R package.   |
| vignettes   | The example codes  |    |

