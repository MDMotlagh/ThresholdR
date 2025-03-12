#' Claculate BIC values
#'
#' @param data A data frame of CLR normalized ADT expression values. Cells are in the rows and ADT probes (surface markers) are in columns.
#' @param k The maximimum expected number of underlying components across all markers. Default is set to 3.
#' @returns A list of numerical vectors of BIC values corresponding to the markers. If BIC value calculation for a marker encounters an error, the function returns NA for that marker in the list.
#' @examples
#' bic_values <- calculateBIC(data = adt, k = 3)



calculateBIC <- function(data, k = 3) {
  require(mixtools)
  require(mclust)
  require(foreach)

  # Initialize vector to store BIC values
  bic_values <- c()

  # Loop over columns (markers) of the data
  bic_values <- foreach(i = colnames(data), .packages = "mclust") %do% {
    # Try-catch to handle errors
    tryCatch({
      print(i)
      temp_bic_values <- numeric(k)

      for (j in 1:k) {
        print(j)
        # Fit the Mclust model for each value of G
        model <- Mclust(data[data[, i] > 0, i], G = j)
        # Store the BIC values, excluding NA values
        temp_bic_values[j] <- max(model$BIC[!is.na(model$BIC)])
      }

      # Return calculated BIC values for the marker
      return(temp_bic_values)

    }, warning = function(w) {
      # Print warning and continue
      message(paste("Warning occurred for", i, "\nWarning:", conditionMessage(w)))
      invokeRestart("muffleWarning")  # Continue execution without stopping at the warning
    }, error = function(e) {
      # Print error and continue
      message(paste("Error occurred for", i, "\nError:", conditionMessage(e)))
      return(NA)  # Return NA for the BIC value in case of an error
    })
  }

  # Assign column names to the bic_values list
  names(bic_values) <- colnames(data)

  # Print errors and warnings (not stored as variables anymore)
  if (length(bic_values) != length(colnames(data))) {
    missing_markers <- setdiff(colnames(data), names(bic_values))
    if (length(missing_markers) > 0) {
      message(paste("Error occurred for the following markers:", paste(missing_markers, collapse = ", ")))
    }
  }

  # Return only bic_values
  return(bic_values)
}
