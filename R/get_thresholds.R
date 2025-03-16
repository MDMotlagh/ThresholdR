#' Calculate different thresholds for bi- and tri-modal fittings
#'
#' @param fittings A list of fitted parameters for each distribution (component) under each k across all markers.
#' @param k_list A list with set of possible k values for each marker as an element.
#' @returns A data frame with calculated thresholds: bi_thr1 defines the mean+3SD of the first component and bi_cut refers to the cutpoint on x axis where the two components meet in a bimodal (k=2) distribution. tri_thr1 and tri_thr2 indicate mean+3SD of the first and second components in a trimodal (k=3) distribution. Subsequently, tri_cut1 and tri_cut2 represent the cutpoints between first and second component, and second and third components, repectively, in a trimodal distribution. The user should choose which thereshold to use the recomended value is mean+3SD of the first component in k=2.
#' @examples
#' thresholds <- get_thresholds(fittings = fittings, k_list = k.list)



get_thresholds <- function(fittings, k_list){
  suppressPackageStartupMessages(require(AdaptGauss))
  thresholds <- data.frame(row.names = names(fittings),
                           bi_thr = rep(NA, length(names(fittings))),
                           bi_cut = rep(NA, length(names(fittings))),
                           tri_thr1 = rep(NA, length(names(fittings))),
                           tri_cut1 = rep(NA, length(names(fittings))),
                           tri_thr2 = rep(NA, length(names(fittings))),
                           tri_cut2 = rep(NA, length(names(fittings))))
  for (i in names(fittings)) {
    for (j in 1:length(k_list[[i]])) {
      if (k_list[[i]][j]==2) {
        thresholds[i,"bi_thr"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
        thresholds[i,"bi_cut"] <- AdaptGauss::Intersect2Mixtures(Mean1 = min(fittings[[i]][[j]]$mu), SD1 = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu==min(fittings[[i]][[j]]$mu))], Weight1 = fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu==min(fittings[[i]][[j]]$mu))],Mean2 = max(fittings[[i]][[j]]$mu), SD2 = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu==max(fittings[[i]][[j]]$mu))], Weight2 = fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu==max(fittings[[i]][[j]]$mu))])$CutX

      } else {
        if (k_list[[i]][j]==3) {
          thresholds[i,"tri_thr1"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
          thresholds[i,"tri_thr2"] <- fittings[[i]][[j]]$mu[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu)&fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))]+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu)&fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))]
          thresholds[i,"tri_cut1"] <- AdaptGauss::Intersect2Mixtures(Mean1 = min(fittings[[i]][[j]]$mu), SD1 = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu==min(fittings[[i]][[j]]$mu))], Weight1 = fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu==min(fittings[[i]][[j]]$mu))],Mean2 = fittings[[i]][[j]]$mu[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))], SD2 = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))], Weight2 = fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))])$CutX
          thresholds[i,"tri_cut2"] <- AdaptGauss::Intersect2Mixtures(Mean1 = fittings[[i]][[j]]$mu[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))], SD1 = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))], Weight1 = fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))],Mean2 = max(fittings[[i]][[j]]$mu), SD2 = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu==max(fittings[[i]][[j]]$mu))], Weight2 = fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu==max(fittings[[i]][[j]]$mu))])$CutX

        }
      }
    }
  }
  return(thresholds)
}
