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
fittings_path <- "/Users/ataraskin/Git/DataAnalysisProjects/CITE-Seq/Melissa_Human_CITE-Seq_2023/Modified_pipeline/Objects/ThresholdeR/Fittings"
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

for (marker in main_markers) {
  p <- ggplot(adt[adt[, marker] > 0, ], aes(x = adt[adt[, marker] > 0, marker])) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = paste("Histogram of", marker, "Values"),
       x = "Marker Values",
       y = "Count")
  print(p)
}

#For CD3 and CD19, we use the domain knowledge that these abs have bimodal distributions

############################################################
# Margin of threshold estimation (values in percents)
###########################################################
main_margins <- defineMargin(main_markers, bic_values) # in percents
all_margins <- defineMargin(names(bic_values), bic_values) # in percents


############################################################
# Get all possible fits
###########################################################

# k_list output usually shorter than main_markers, but they should match!!!
k_list <- list()
for (n in 1:length(bic_values)) {
  k_list <- c(k_list, GetAllPossibleFits(bic_values = bic_values[n], 
                                         margin_thr = 0.01 # Based on info from main_margins table (NOT percents, part of 1)
                                         ))
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

##########################################################################
# 1. List of antibodies for further analysis. Here ONLY main markers
##########################################################################
markers_with_fits <- intersect(main_markers, names(k_list))
    
fittings <- fit_k(data = adt[ , markers_with_fits], k_list = k_list[markers_with_fits], margin.den = 0.1, epsilon = 0.000001, seed = 42)


thresholds <- get_thresholds(fittings = fittings, k_list = k_list)

##########################################################################
# Publishing
dir_name <- "main_markers_fittings"
path_to_save_fittings <- file.path(fittings_path, dir_name)
if (!dir.exists(path_to_save_fittings)) {
  # Create the directory
  dir.create(path_to_save_fittings)
  cat("Directory created:", path_to_save_fittings, "\n")
} else {
  cat("Directory already exists:", path_to_save_fittings, "\n")
}
publish_fittings_thr(fittings = fittings, k_list =k_list ,thresholds = thresholds, path = path_to_save_fittings, seed = 42)
##########################################################################



##########################################################################
# 2. List of antibodies for further analysis. ALL markers
##########################################################################

markers_with_fits <- intersect(all_markers, names(k_list))

# Error
fittings <- fit_k(data = adt[ , markers_with_fits], 
                  k_list = k_list[markers_with_fits], 
                  margin.den = 0.1, 
                  epsilon = 0.000001, 
                  maxrestarts = 5000,
                  seed = 42)
# Markers, discarded during the current step:
discarded_markers <- setdiff(markers_with_fits, names(fittings))

bad_abs <- unique(c(bad_abs, discarded_markers))
length(bad_abs)
  
thresholds_all <- get_thresholds(fittings = fittings, k_list = k_list)


##########################################################################
# Publishing
dir_name <- "all_markers_fittings"
path_to_save_fittings <- file.path(fittings_path, dir_name)
if (!dir.exists(path_to_save_fittings)) {
  # Create the directory
  dir.create(path_to_save_fittings)
  cat("Directory created:", path_to_save_fittings, "\n")
} else {
  cat("Directory already exists:", path_to_save_fittings, "\n")
}
publish_fittings_thr(fittings = fittings,
                     k_list = k_list,
                     thresholds = thresholds,
                     path = path_to_save_fittings,
                     seed = 42)

# thresholds_table <- look_fittings_thr(fittings = fittings,
#                                        k_list = k_list,
#                                        thresholds = thresholds_all,
#                                        path = path_to_save_fittings,
#                                        seed = 42)
##########################################################################



##########################################################################
# 3. Automatic processing. ALL markers
##########################################################################


# Count NA values in each row
na_counts <- rowSums(is.na(thresholds_all))

# Separate rows based on NA count
one_value <- thresholds_all[na_counts == 1, ]
two_values <- thresholds_all[na_counts == 0, ]

# No need to filter
one_value$threshold <- apply(one_value, 1, function(row) row[!is.na(row)])
set_of_thresholds_1 <- one_value["threshold"]


two_values$threshold <- ifelse(two_values$bimodal_thr>=two_values$trimodal_thr, two_values$bimodal_thr, two_values$trimodal_thr)
set_of_thresholds_2 <- two_values["threshold"]

set_of_thresholds <- rbind(set_of_thresholds_1, set_of_thresholds_2)
nrow(set_of_thresholds)


# Saving the results as .RDS
saveRDS(set_of_thresholds, file.path(object_dir_current, "set_of_thresholds.rds"))



















