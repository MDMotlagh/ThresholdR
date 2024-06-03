# Free the RAM
vars <- ls()
rm(list = vars)
rm(vars)
gc()

# Libraries
library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)
set.seed(42)


###########################################################################
# Paths and loading the object
###########################################################################
object_dir <- "/Users/ataraskin/Git/DataAnalysisProjects/CITE-Seq/Melissa_Human_CITE-Seq_2023/Modified_pipeline/Objects/MT_BPC"
object_dir_current <- "/Users/ataraskin/Git/DataAnalysisProjects/CITE-Seq/Melissa_Human_CITE-Seq_2023/Modified_pipeline/Objects/ThresholdeR"
folder_with_plots <- "/Users/ataraskin/Git/DataAnalysisProjects/CITE-Seq/Melissa_Human_CITE-Seq_2023/Modified_pipeline/Objects/ThresholdeR/plots"
object_postfix <- "_DoubletFinder_Demulti_MT_BPC"
setwd("/Users/ataraskin/Git/DataAnalysisProjects/CITE-Seq/Melissa_Human_CITE-Seq_2023/Modified_pipeline/Objects/ThresholdeR") 
source("/Users/ataraskin/Git/DataAnalysisProjects/CITE-Seq/Melissa_Human_CITE-Seq_2023/Modified_pipeline/05_ThresholdeR/functions.R") 
demult.dblfinder <- readRDS("test.RDS") 
DefaultAssay(demult.dblfinder) <- "ADT"
###########################################################################



###########################################################################
# Normalization
demult.dblfinder <- NormalizeData(demult.dblfinder, assay="ADT", normalization.method="CLR", margin=2)

adt <- as.data.frame(t(as.matrix(demult.dblfinder@assays$ADT@data)))
###########################################################################  

# main_markers could be different 
main_markers <- c("CD14.1", "CD19.1", "CD3", "CD4.1", "CD8a", "CD16", "CD56")
all_markers <- sort(rownames(demult.dblfinder@assays$ADT@counts))

###########################################################################
# Some info
###########################################################################
# markers_with_fits <- c("CD14.1", "CD19.1", "CD3", "CD4.1", "CD8a")
########################
# AB causing the problem
# main_markers <- main_markers[!(main_markers == "CD120b")]
# main_markers <- main_markers[!(main_markers == "Siglec-10")]
###########################################################################


bic_values <- c()
error_abs <- c()
warning_abs <- c()


markers <- all_markers
suppressWarnings({
  
for (num in seq_along(markers)) {
  tryCatch({
    withCallingHandlers({
      invisible(capture.output({
        adt <- as.data.frame(t(as.matrix(demult.dblfinder@assays$ADT@data)))
        marker_data <- as.data.frame(adt[, markers[[num]], drop = FALSE]) 
        bic_values[markers[[num]]] <- calculateBICvalues(marker_data, k = 4, parallel = FALSE)
      }))
    }, warning = function(w) {
      # Without "<<-" sign instead of "<-" you won't be able to access those vectors outside the loop.
      # Say thanks to R.
      warning_abs <<- c(warning_abs, markers[[num]])
      message(paste("Warning occurred for ", markers[[num]], "\nWarning:", conditionMessage(w)))
      invokeRestart("muffleWarning")  # to continue execution without stopping at the warning
    })
  }, error = function(e) {
    # Without "<<-" sign instead of "<-" you won't be able to access those vectors outside the loop.
    # Say thanks to R.
    error_abs <<- c(error_abs, markers[[num]])
    message(paste("Error occurred for ", markers[[num]], "\nError:", conditionMessage(e)))
    bic_values[markers[[num]]] <- NA  # value indicating the error
  })
}

})

print(paste("Antibodies that are causing errors:", paste(error_abs, collapse = ", ")))
print(paste("Antibodies that are causing warnings:", paste(warning_abs, collapse = ", ")))
print(paste("Total amount of BIC values are:", length(bic_values)))
print(paste("The amount of ABs for which BIC values cannot be obtained:", length(all_markers)-length(bic_values)))



for (marker in main_markers) {
  
  p <- ggplot(adt[adt[,marker]>0,], aes(x=adt[adt[,marker]>0, marker])) +
    geom_density(adjust=2)+
    labs(title = paste(marker))
  print(p)
}

#For CD3 and CD19, we use the domain knowledge that these abs have bimodal distributions

############################################################
# Margin of threshold estimation
###########################################################
all_margins <- defineMargin(main_markers, bic_values) # in percents

############################################################
# Get all possible fits
###########################################################

# k_list output usually shorter than main_markers, but they should match!!!
k_list <- list()
for (n in 1:length(bic_values)) {
  k_list <- c(k_list, GetAllPossibleFits(bic_values = bic_values[n], margin_thr = 0.10))
}


# DIFFERENT
length(bic_values)
length(k_list)
length(main_markers)
names(k_list)

#####################################
# Additional list of "bad" antibodies
#####################################
bad_abs <- setdiff(all_markers, names(k_list))
bad_abs <- unique(c(bad_abs, error_abs, warning_abs))
length(bad_abs)

# List of antibodies for further analysis
markers_with_fits <- intersect(main_markers, names(k_list))
    
fitting <- fit_k(data = adt[ , markers_with_fits], k_list = k_list[markers_with_fits], margin.den = 0.1, epsilon = 0.000001, seed = 42)

# Error
thresholds <- get_thresholds(fittings = fittings, k_list = k_list)

publish_fittings_thr(fittings = fittings, k_list =k_list ,thresholds = thresholds, path = "03 Fittings/", seed = 42)
