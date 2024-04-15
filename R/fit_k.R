#' Fit for all Ks from k_list
#'
#' This function fits for all possible Ks extracted through GetAllPossibleFits function.
#' Fittig is done through normalmixEm function from mixtools R package
#'
#'
#' @param data dataframe of normalized ADT values: Cells as rows and Protein as columns
#' @param k_list a list of vectors, each vector contains values for possible components fitted best between k=1 to k=k for normalized distribution of each protein
#' @param margin.den a margin to evaluate the density of the first component at its mean value (x=mean(C1)) whether it is at most 0.1 times higher than the density of the whole population at x=mean(C1).
#' if density of C1 is greater than 1.1*density of the whole data at the point of x=mean(C1), a correction step is needed for the standard deviation of the first component (C1), because in some cases that left side of the distribution is truncated, SD is underestimated which leads to overshooting of the density of the first component.
#' @param epsilon inherent to normalmixEM function: the convergence criterion.
#' @param maxit The maximum number of iterations (Default=10000)
#' @param maxrestarts The maximum number of restarts (Default=5)
#'
#' @return Returns a list, `fittings.corrected` including a list for each k fittings with values:
#'
#' * `x` - The normalized data.
#' * `lambda` - mixing weights or probabilities.
#' * `mu` - The mean parameters
#' * `sigma` - The standard deviation parameters.
#' * `loglik` - log-likelihood.
#' * `posterior` - A matrix of posterior probabilities.
#' * `all.loglik` - A vector of each iteration's log-likelihood
#' * `restarts` - The number of restarts.
#' * `ft` - function name, a character vector.
#'
#' @examples
#' data("adt.df")
#' fittings.test <- fit_k(data = adt.df[,1:3], k_list = k_list.test, margin.den = 0.1, seed = 42)
#'#' @export
#'
fit_k <- function(data, k_list, margin.den=0.1,epsilon = 0.01, maxit = 10000, maxrestarts = 5, seed) {
  require(mixtools)
  require(mclust)

  fittings <- list()
  set.seed(seed)
  ggdensity <- c()
  dfs <- c()
  dfs.corrected <- c()
  fittings.corrected <- c()
  for (i in colnames(data)) {
    for (j in 1:length(k_list[[i]])) {
      cat(paste0(i, ": ",k_list[[i]][j]," modal: "))
      fittings[[i]][[j]] <- normalmixEM(data[, i], k = k_list[[i]][j], epsilon = epsilon, maxit = maxit, maxrestarts = maxrestarts)
    }
  }
  for (i in colnames(data)) {
    ggdensity[[i]] <- ggplot_build(ggplot(data, aes(x=data[,i]))+
                                     geom_density(adjust=1))
    for (j in 1:length(k_list[[i]])) {
      if (k_list[[i]][j] == 2) {
        x <- fittings[[i]][[j]]$x
        dfs[[i]][[j]] <- data.frame(x = x,
                                    density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                                    density2 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
        print(paste("Adding gg.x and gg.y for bimodal",i,  "DF-No:", which(names(dfs) == i)))
        dfs[[i]][[j]]$gg.x <- sapply(dfs[[i]][[j]]$x, find_closest_value, reference_vector = ggdensity[[i]]$data[[1]]$x)
        dfs[[i]][[j]]$gg.y <- ggdensity[[i]]$data[[1]]$y[match(dfs[[i]][[j]]$gg.x, ggdensity[[i]]$data[[1]]$x)]
        dfs.corrected[[i]][[j]] <- dfs[[i]][[j]]
        fittings.corrected[[i]][[j]] <- fittings[[i]][[j]]
        # Find the closes x value in gg.x
        x.mean <- find_closest_x(value = min(fittings[[i]][[1]]$mu), reference_vector = fittings[[i]][[1]]$x)
        C1.y <- unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$density1)
        if (C1.y > (1+margin.den)*unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$gg.y)) {
          x <- fittings[[i]][[j]]$x
          sd1 <- fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
          cat(paste("Fixing Sigma for C1 in",i,"-Bimodal"))
          while (C1.y > unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$gg.y)) {
            sd1 <- sd1 + 0.005
            fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu == min(fittings.corrected[[i]][[j]]$mu))] <- sd1
            dfs.corrected[[i]][[j]] <- data.frame(x = x,
                                                  density1 = dnorm(x, mean = min(fittings.corrected[[i]][[j]]$mu), sd = fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu == min(fittings.corrected[[i]][[j]]$mu))]) * fittings.corrected[[i]][[j]]$lambda[which(fittings.corrected[[i]][[j]]$mu == min(fittings.corrected[[i]][[j]]$mu))],
                                                  density2 = dnorm(x, mean = max(fittings.corrected[[i]][[j]]$mu), sd = fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu == max(fittings.corrected[[i]][[j]]$mu))]) * fittings.corrected[[i]][[j]]$lambda[which(fittings.corrected[[i]][[j]]$mu == max(fittings.corrected[[i]][[j]]$mu))])
            # Find the closest x value in gg.x
            C1.y <- unique(dfs.corrected[[i]][[j]][dfs.corrected[[i]][[j]]$x==x.mean,]$density1)
          }
        }
      }else {
        if (k_list[[i]][j] == 3) {
          x <- fittings[[i]][[j]]$x
          dfs[[i]][[j]] <- data.frame(x = x,
                                      density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                                      density2 = dnorm(x, mean = fittings[[i]][[j]]$mu[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                                                       sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                                      density3 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
          print(paste("Adding gg.x and gg.y for trimodal",i,  "DF-No:", which(names(dfs) == i)))
          dfs[[i]][[j]]$gg.x <- sapply(dfs[[i]][[j]]$x, find_closest_value, reference_vector = ggdensity[[i]]$data[[1]]$x)
          dfs[[i]][[j]]$gg.y <- ggdensity[[i]]$data[[1]]$y[match(dfs[[i]][[j]]$gg.x, ggdensity[[i]]$data[[1]]$x)]
          dfs.corrected[[i]][[j]] <- dfs[[i]][[j]]
          fittings.corrected[[i]][[j]] <- fittings[[i]][[j]]
          # Find the closes x value in gg.x
          x.mean <- find_closest_x(value = min(fittings[[i]][[j]]$mu), reference_vector = fittings[[i]][[j]]$x)
          C1.y <- unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$density1)
          if (C1.y > (1+margin.den)*unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$gg.y)) {
            x <- fittings[[i]][[j]]$x
            sd1 <- fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
            cat(paste("Fixing Sigma for C1 in",i,"-Trimodal"))
            while (C1.y > unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$gg.y)) {
              sd1 <- sd1+0.005
              fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu == min(fittings.corrected[[i]][[j]]$mu))] <- sd1
              dfs.corrected[[i]][[j]] <- data.frame(x = x,
                                                    density1 = dnorm(x, mean = min(fittings.corrected[[i]][[j]]$mu), sd = fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu == min(fittings.corrected[[i]][[j]]$mu))]) * fittings.corrected[[i]][[j]]$lambda[which(fittings.corrected[[i]][[j]]$mu == min(fittings.corrected[[i]][[j]]$mu))],
                                                    density2 = dnorm(x, mean = fittings.corrected[[i]][[j]]$mu[which(fittings.corrected[[i]][[j]]$mu != min(fittings.corrected[[i]][[j]]$mu) & fittings.corrected[[i]][[j]]$mu != max(fittings.corrected[[i]][[j]]$mu))],
                                                                     sd = fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu != min(fittings.corrected[[i]][[j]]$mu) & fittings.corrected[[i]][[j]]$mu != max(fittings.corrected[[i]][[j]]$mu))]) * fittings.corrected[[i]][[j]]$lambda[which(fittings.corrected[[i]][[j]]$mu != min(fittings.corrected[[i]][[j]]$mu) & fittings.corrected[[i]][[j]]$mu != max(fittings.corrected[[i]][[j]]$mu))],
                                                    density3 = dnorm(x, mean = max(fittings.corrected[[i]][[j]]$mu), sd = fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu == max(fittings.corrected[[i]][[j]]$mu))]) * fittings.corrected[[i]][[j]]$lambda[which(fittings.corrected[[i]][[j]]$mu == max(fittings.corrected[[i]][[j]]$mu))])
              # Find the closes x value in gg.x
              C1.y <- unique(dfs.corrected[[i]][[j]][dfs.corrected[[i]][[j]]$x==x.mean,]$density1)

            }
          }
        }
      }
    }
  }
  return(fittings.corrected)
}
