The PBMC10K dataset was downloaded from: https://www.10xgenomics.com/datasets/10k-human-pbmcs-stained-with-totalseq-C-human-TBNK-cocktail-GEM-X.
```r
pbmc10k <- read10x_h5('10k_Human_PBMC_TotalSeqC_5p_gemx_10k_Human_PBMC_TotalSeqC_5p_gemx_count_sample_filtered_feature_bc_matrix.h5')
```

```r
Bmcite <- SeuratData::LoadData('bmcite')
```
