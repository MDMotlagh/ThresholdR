
#' Get All Possible Fits
#'
#' This function extracts the possible fits for distribution of each protein, by evaluating in each pair of subsequent
#' fittings whether k=i+1 has greater BIC value than k=i.
#'
#'
#'
#'
#' @param bic_values a list of vectors, each contatining BIC values for k number of fittings for each protein distribution
#' @param margin_thr a minimum margin threshold for BIC of k=i+1 to be greater than BIC of k=i (BIC(k=i+1)>(1+margin_thr)*BIC(k=i))
#'
#' @return Returns a list, `k_list` of vectors:
#'
#' * each vector contains values for possible components fitted best between k=1 to k=k
#'
#'
#' @examples
#' data("adt.df")
#' GetAllPossibleFits(bic_values, margin_thr=0.1)
#' @export
#'
#'
getAllPossibleFits <- function(bic_values, margin_thr) {
  k_list <- c()
  bic_values_indices <- lapply(bic_values, function(vec) {
    indices <- NULL
    if (all(vec >= 0)) {
      sapply(1:(length(vec) - 1), function(i) {
        percent_diff <- (vec[i + 1] - vec[i]) / vec[i]
        if (percent_diff >= margin_thr) {
          indices <- c(indices, i + 1)
        }
      })
    } else {
      if (all(vec <= 0)) {
        sapply(1:(length(vec) - 1), function(i) {
          percent_diff <- (vec[i] - vec[i + 1]) / vec[i]
          if (percent_diff >= margin_thr) {
            indices <- c(indices, i + 1)
          }
        })
      }else {
        positive_values <- vec[vec >= 0]
        if (length(positive_values)>1) {
          sapply(1:(length(positive_values) - 1), function(i) {
            percent_diff <- (positive_values[i + 1] - positive_values[i]) / positive_values[i]
            if (percent_diff >= margin_thr) {
              indices <- c(indices, which(vec == positive_values[i + 1]))
            }
          })
        }else {
          indices <- c(indices, which(vec == positive_values[1]))
        }
      }

    }
  })
  for (i in names(bic_values_indices)) {
    if (class(bic_values_indices[[i]])=="list") {
      for (j in 1:length(bic_values_indices[[i]])) {
        if (!is.null(bic_values_indices[[i]][[j]])) {
          k_list[[i]][j] <- bic_values_indices[[i]][[j]]
        }
      }
    }else {
      if (class(bic_values_indices[[i]])=="numeric") {
        k_list[[i]] <- bic_values_indices[[i]]
      }
    }
  }
  return(k_list)
}
