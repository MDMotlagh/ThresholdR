
#' Get Thresholds
#'
#' This function calculates the threshold based on mean and standard deviation of the
#' first component (M1+3SD)
#'
#'
#' @param fittings dataframe of normalized ADT values: Cells as rows and Protein as columns
#' @param k_list a list of vectors, each vector contains values for possible components fitted best between k=1 to k=k for normalized distribution of each protein
#'
#' @return a dataframe, `thresholds` including columns for thresholds of each fitting, namely:
#'
#' * `bimodal_thr` - thresold for bimodal (k=2) fittings
#' * `trimodal_thr` - thresold for trimodal (k=3) fittings
#'
#'
#' @examples
#' data("adt.df")
#' thresholds.test <- get_thresholds(fittings = fittings.test, k_list = k_list.test)
#'
get_thresholds <- function(fittings, k_list){
  thresholds <- data.frame(row.names = names(bic_values),
                           bimodal_thr = rep(NA, length(names(bic_values))),
                           trimodal_thr = rep(NA, length(names(bic_values))))
  for (i in names(bic_values)) {
    for (j in 1:length(k_list[[i]])) {
      if (k_list[[i]][j]==2) {
        thresholds[i,"bimodal_thr"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
      }else {
        if (k_list[[i]][j]==3) {
          thresholds[i,"trimodal_thr"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
        }
      }
    }
  }
  return(thresholds)
}
