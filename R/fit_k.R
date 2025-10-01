#' Fit for the possible number of underlying components (k) across all markers
#'
#' @param data A data frame of CLR normalized ADT expression values. Cells are in the rows and ADT probes (surface markers) are in columns.
#' @param k_list A list of the possible number of underlying components across all markers.
#' @param margin.den The margin for difference between the density of the mean in the first fitted component and the density value of the input data at x=mean(C1). This is needed to correct the underestimated standard deviation of the first component which happens sometimes with normalmixEM function from mixtools package.
#' @param maxit Maximum number of iterations for fitting of each k. The argument come from the 'mixtools' package.
#' @param maxrestarts Maximum number of restarts for fitting of each k. The argument come from the 'mixtools' package.
#' @param seed Random number generator. Is set to a user provided value to promote reproducibility.
#' @returns A list of parameters such as mean, standard-deviation, weight (lambda) for each component in different k values across the markers. This will later be used to define thresholds.
#' @examples
#' fittings <- fit_k(data = adt, k_list = k.list, margin.den = 0.1, epsilon = 0.01, maxit = 1000, maxrestarts = 10)



fit_k <- function(data,
                  k_list,
                  margin.den=0.1,
                  epsilon = 0.01,
                  maxit = 10000,
                  maxrestarts = 20,
                  seed=42) {
  suppressPackageStartupMessages(require(mixtools))
  suppressPackageStartupMessages(require(mclust))
  suppressPackageStartupMessages(require(ggplot2))

  set.seed(seed)
  fittings <- list()
  corrected_fittings <- list()

    # Function to optimize sigma for component 1 to match empirical density
  optimize_sigma <- function(mu, lambda, target_y, x_val, init_sd = 1) {
    objective <- function(sd1) {
      d1_y <- dnorm(x_val, mean = mu, sd = sd1) * lambda
      (d1_y - target_y)^2
    }
    res <- optim(par = init_sd, fn = objective, method = "L-BFGS-B", lower = 1e-4, upper = 10)
    return(res$par)
  }

  for (var in names(k_list)) {
    variable_data <- data[data[, var] > 0, var]
    kde <- density(variable_data, from = 0)

    fittings[[var]] <- list()
    corrected_fittings[[var]] <- list()

    for (k_idx in seq_along(k_list[[var]])) {
      k <- k_list[[var]][k_idx]
      cat(paste0(var, ": ", k, " components: "))


      # Fit GMM
      fit <- tryCatch(
        {
          normalmixEM(variable_data,
                      k = k,
                      epsilon = epsilon,
                      maxit = maxit,
                      maxrestarts = maxrestarts)
          },
        error = function(e) {
          cat(paste("Error - skipping", k, "components for", var, "\n"))
          return(NULL)
        }
        )

      if (is.null(fit)) next

      fittings[[var]][[k_idx]] <- fit
      corrected_fit <- fit  # Copy to allow correction

      # Focus only on fixing the lowest-mean component
      mu1_idx <- which.min(fit$mu)
      mu1 <- fit$mu[mu1_idx]
      lambda1 <- fit$lambda[mu1_idx]
      sd1_init <- fit$sigma[mu1_idx]

      # Get empirical KDE at mu1
      emp_y_at_mu1 <- approx(kde$x, kde$y, xout = mu1)$y

      # Model-predicted density at mu1
      model_y_at_mu1 <- dnorm(mu1, mean = mu1, sd = sd1_init) * lambda1

      if (model_y_at_mu1 > (1 + margin.den) * emp_y_at_mu1 ||
          model_y_at_mu1 < (1 - margin.den) * emp_y_at_mu1) {
        cat("Fixing sigma for C1 in", var, "- k =", k, "\n")

        sd1_new <- optimize_sigma(
          mu = mu1,
          lambda = lambda1,
          target_y = emp_y_at_mu1,
          x_val = mu1,
          init_sd = sd1_init
          )

        # Update corrected fit
        corrected_fit$sigma[mu1_idx] <- sd1_new
      }

      corrected_fittings[[var]][[k_idx]] <- corrected_fit
    }
  }

  return(corrected_fittings)
}
