# Overview
ThresholdR uses Gaussian Mixture Models concept to identify the noise population in each surface marker across cells in CITE-seq experiments. Once the noise distribution is identified, it calculates the upper threshold of the noise component to separate expressing and non-expressing cells for each surface marker. 
# Citation
# Installation
```r
#1. Install required packages:
required_packages <- c("mixtools", "mclust", "foreach", "ggplot2", "AdaptGauss")

# Check and install missing packages
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
#2. Install ThresholdR:
devtools::install_github('MDMotlagh/ThresholdR')
```
# Instructions
## Contents
| Folder | Description | Details |
| --------------- | --------------- | --------------- |
| data   | bmcite dataset   | bmcite is available through SeuratData R package.   |
| vignettes   | The example codes  |  Vignette/bmcite_example.Rmd and Vignette/bmcite_example.html  |
| output   | Fitting plots and threshold values  |  Vignette/02Plots/  and Vignette/03AllThresholds.csv|

