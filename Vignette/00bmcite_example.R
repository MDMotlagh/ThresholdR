library(Seurat)
set.seed(42)
#so <- SeuratData::LoadData('bmcite')
so <- readRDS('Projects/Packages/ThresholdR/bmcite_so_v4.RDS')
# set the default assay to ADT
#so[["ADT4"]] <- as(object = so[["ADT"]], Class = "Assay")
#DefaultAssay(so) <- "ADT4"
# normalize the data
#DefaultAssay(so) <- "ADT"

so <- NormalizeData(so, normalization.method="CLR", margin=2)

#adt <- as.data.frame(t(as.matrix(so@assays$ADT@data)))
#adt <- as.data.frame(t(as.matrix(so@assays$ADT4@data)))
setwd('/Users/oliaei/Projects/Packages/ThresholdR/R/')
source('calculateBIC.R')
source('fit_k.R')
source('get_k_list.R')
source('get_thresholds.R')
source('publish_fittings_thr.R')

# we assume that the maximum number of underlying poulations (components) across all ADT amrkers is equal to 3 (trimodal). Then we will calculate BIC values which indicate how well a model is fitted. More positive BIC value in this setting means better fitting:
bic.values <- calculateBIC(data = adt, k = 3)
# Next, to choose the possible fits (k=2 and/or k=3) for each ADT, the user should assign a set of few markers which they believe to have bimodal or trimodal (non-unimodal) distributions. With this input, a margin is defined as the minumum value of BIC2/BIC1 devided by 3 across the set of user-provided main markers. If BIC(k+1) is greater BIC(k) by at least that calculated margin in any ADT marker, k+1 is also included in the list of possible k number of components for the corresponding marker.
k.list <- get_k_list(bic_values = bic.values, main.markers = c("CD3", "CD14", "CD19", "CD45RA"))
# If the user wishes to include k=2 and k=3 for all markers in the data matrix, they could modify k_list as the following:
# k.list <- list()
# for (i in names(adt)) {
#   k.list[[i]] <- c(2,3)
# }

#to see if a marker is not meeting the threshold for any k, and is dropped out of the k_list, we can have a look:
#names(bic.values)[!names(bic.values) %in% names(k.list)]


# with the estimated list of possible k(s) for markers, we can run the fitting algorithm:
fittings <- fit_k(data = adt, k_list = k.list, margin.den = 0.1, epsilon = 0.01, maxit = 1000, maxrestarts = 10)
# we then proceed to calculate different thresholds in bimodal and trimodal fittings
thresholds <- get_thresholds(fittings = fittings, k_list = k.list)
publish_fittings_thr(fittings = fittings, k_list = k.list, thresholds = thresholds, path = '/Users/oliaei/Projects/Packages/ThresholdR/Vignette/02Plots/', seed = 42)

markers <- c("CD3", "CD4", "CD8a", "CD14", "CD16", "CD19", "CD56")
markers.thr <- data.frame(row.names = markers,
                          value=rep(NA, length(markers)))
markers.thr[c("CD3", "CD56"),"value"] <- thresholds[c("CD3", "CD56"),"bi_thr"]
markers.thr[c("CD14", "CD19"),"value"] <- thresholds[c("CD14", "CD19"),"bi_cut"]
markers.thr[c("CD4", "CD8a"),"value"] <- thresholds[c("CD4", "CD8a"),"tri_thr1"]
markers.thr["CD16","value"] <- thresholds["CD16","tri_thr2"]
write.csv(thresholds, '/Users/oliaei/Projects/Packages/ThresholdR/Vignette/03AllThresholds.csv')
write.csv(markers.thr, '/Users/oliaei/Projects/Packages/ThresholdR/Vignette/04Deciding_Markers.csv')

thAb <- adt[,markers]
for (i in markers) {
  thAb[,i] <- thAb[,i]-markers.thr[i,"value"]
  thAb[thAb[,i]<0,i] <- 0
}

so[["thAb"]] <- CreateAssayObject(data=t(as.matrix(thAb)))
DefaultAssay(so) <- "thAb"
so$annotation <- ifelse(colnames(so) %in% WhichCells(so, expression = CD3==0 & CD4==0 & CD8a==0 & CD14==0 & CD16==0 & CD19>0 & CD56==0), "B",
                        ifelse(colnames(so) %in% WhichCells(so, expression = CD3>0 & CD4>0 & CD8a==0 & CD14==0 & CD16==0 & CD19==0 & CD56==0), "CD4T",
                               ifelse(colnames(so) %in% WhichCells(so, expression = CD3>0 & CD4==0 & CD8a>0 & CD14==0 & CD19==0 & CD56==0), "CD8T",
                                      ifelse(colnames(so) %in% WhichCells(so, expression = CD3==0 & CD14>0 & CD16==0 & CD19==0 & CD56==0), "cMo",
                                             ifelse(colnames(so) %in% WhichCells(so, expression = CD3==0 & CD14>0 & CD16>0 & CD19==0 & CD56==0), "iMo",
                                                    ifelse(colnames(so) %in% WhichCells(so, expression = CD3==0 & CD14==0 & CD16>0 & CD19==0 & CD56==0), "nMo",
                                                           ifelse(colnames(so) %in% union(WhichCells(so, expression = CD3==0 & CD14==0 & CD16==0 & CD19==0 & CD56>0),
                                                                                          WhichCells(so, expression = CD3==0 & CD14==0 & CD16>0 & CD19==0 & CD56>0)), "NK",
                                                                  "Rem")))))))
table(so$annotation)
DefaultAssay(so) <- "ADT4"
library(dplyr)
so <- so %>%
  NormalizeData(normalization.method = "CLR", margin = 2) %>%
  ScaleData(features = rownames(so)) %>%
  RunPCA(features = rownames(so))
ElbowPlot(so)
so <- so %>%
  FindNeighbors(dims = 1:20)%>%
  FindClusters(resolution = 0.3, random.seed = 42) %>%
  RunUMAP(dims = 1:20, spread = 1.5, min.dist = 0.5, seed.use = 42)
DimPlot(so)
DefaultAssay(so) <- "thAb"
FeaturePlot(so, features = rownames(so), interactive = T)

DimPlot(so, group.by = "annotation")




