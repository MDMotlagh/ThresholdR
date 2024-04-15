#' Calculate the BIC values
#'
#' This function calculates BIC values for each fitting from k=1 to k=k using Mclust function from
#' mclust R package. The higher BIC value points to better fit through mclust package.
#'
#'
#'
#'
#' @param data dataframe of normalized ADT values: Cells as rows and Protein as columns
#' @param k maximum number of underlying components (cell populations) to consider fitting through Gaussian Mixture Models (3 by default)
#' @param parallel to run the process in parallel (FALSE by default)
#'
#' @return Returns a list, `bic_values` of vectors each corresponding to a column in the input dataframe:
#'
#' * each vector has k number observations corresponding to k=1 to k=k (maximum number of underlying components)
#'
#'
#' @examples
#' data("adt.df")
#' calculateBICvalues(adt.df, k=3, parallel = FALSE)
#'#' @export
#'
calculateBICvalues <- function(data, k=3, parallel = FALSE) {
  require(foreach)
  require(doParallel)

  if (parallel) {
    # Initialize cluster
    cl <- makeCluster(detectCores()-2)
    registerDoParallel(cl)
  }

  bic_values <- if (parallel) {
    foreach(i = colnames(data), .packages="mclust") %dopar% {
      print(i)
      temp_bic_values <- numeric(k)
      for (j in 1:k) {
        print(j)
        model <- Mclust(data[, i], G = j)
        temp_bic_values[j] <- max(model$BIC)
      }
      return(temp_bic_values)
    }
  } else {
    foreach(i = colnames(data), .packages="mclust") %do% {
      print(i)
      temp_bic_values <- numeric(k)
      for (j in 1:k) {
        print(j)
        model <- Mclust(data[, i], G = j)
        temp_bic_values[j] <- max(model$BIC)
      }
      return(temp_bic_values)
    }
  }

  if (parallel) {
    stopCluster(cl)
  }

  names(bic_values) <- colnames(data)
  return(bic_values)
}
