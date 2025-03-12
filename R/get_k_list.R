#' The possible number of underlying components across all markers
#'
#' @param bic_values A data frame of CLR normalized ADT expression values. Cells are in the rows and ADT probes (surface markers) are in columns.
#' @param main.markers A user-defined set of few markers with known bimodal or trimodal distributions.
#' @returns A list with set of possible k values for each marker as an element.
#' @examples
#' k.list <- get_k_list(bic_values = bic.values, main.markers = c("CD3", "CD14", "CD19", "CD45RA"))


get_k_list <- function(bic_values, main.markers) {
  require(mixtools)
  require(mclust)
  # Initialize an empty data frame for BIC ratios
  bic_ratios <- data.frame()
  k_list <- c()

  # Set up the dimensions and names for the data frame
  bic_ratios[1:length(bic_values), 1:2] <- NA
  rownames(bic_ratios) <- names(bic_values)
  colnames(bic_ratios) <- c("BIC2overBIC1", "BIC3overBIC2")

  # Calculate the BIC ratios
  for (i in names(bic_values)) {
    print(i)
    bic_ratios[i, "BIC2overBIC1"] <- abs((bic_values[[i]][2] - bic_values[[i]][1]) / bic_values[[i]][1])
    bic_ratios[i, "BIC3overBIC2"] <- abs((bic_values[[i]][3] - bic_values[[i]][2]) / bic_values[[i]][2])

  }
  margin.thr <- min(bic_ratios[main.markers,"BIC2overBIC1"])/3
  for (i in names(bic_values)) {
    if (bic_ratios[i, "BIC2overBIC1"]>margin.thr&bic_ratios[i, "BIC3overBIC2"]>margin.thr) {
      k_list[[i]] <- c(k_list[[i]],c(2,3))
    }else {
      if (bic_ratios[i, "BIC2overBIC1"]>margin.thr&bic_ratios[i, "BIC3overBIC2"]<margin.thr) {
        k_list[[i]] <- c(k_list[[i]],c(2))
      }else {
        if (bic_ratios[i, "BIC2overBIC1"]<margin.thr&bic_ratios[i, "BIC3overBIC2"]>margin.thr) {
          k_list[[i]] <- c(k_list[[i]],c(3))
        }else {
          k_list[[i]] <- NULL
        }
      }
    }
  }
  return(k_list)
}
