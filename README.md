# Overview
ThresholdR uses Gaussian Mixture Models concept to identify the noise population in each surface marker across cells in CITE-seq experiments. Once the noise distribution is identified, it calculates the upper threshold of the noise component to separate expressing and non-expressing cells for each surface marker. 
# Citation
# Installation
```r
devtools::install_guithub('MDMotlagh/ThresholdR')
```
# Instructions
## Contents
| Folder | Description | Details |
| --------------- | --------------- | --------------- |
| data   | PBMC10K and bmcite datasets   | PBMC10K immune profiling dataset acquired from 10x genomics website. bmcite is available through SeuratData R package.   |
| vignettes   | The example codes  | PBMC10K immune profiling dataset acquired from 10x genomics website. bmcite is available through SeuratData R package.   |

