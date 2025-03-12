#' Fit for the possible number of unrlying components (k) across all markers
#'
#' @param data A data frame of CLR normalized ADT expression values. Cells are in the rows and ADT probes (surface markers) are in columns.
#' @param k_list The possible number of underlying components across all markers.
#' @param margin.den The margin for difference between the density of the mean in first fitted component and the density value of the input data at x=mean(C1). This is needed to correct the underestimated standard deviation of the first component which happens sometimes.
#' @param maxit Maximum number of iterations for fitting of each k. The argument come from the 'mixtools' package.
#' @param maxrestarts Maximum number of restarts for fitting of each k. The argument come from the 'mixtools' package.
#' @param seed Random number generator. Is set to a user provided value to promote reproducibility.
#' @returns A list of parameters such as mean, standard-deviation, weight (lambda) for each k across the markers. This will later be used to define thresholds.
#' @examples
#' fittings <- fit_k(data = adt, k_list = k.list, margin.den = 0.1, epsilon = 0.01, maxit = 1000, maxrestarts = 10)



fit_k <- function(data,
                  k_list,
                  margin.den=0.1,
                  epsilon = 0.01,
                  maxit = 10000,
                  maxrestarts = 20,
                  seed=42) {
  require(mixtools)
  require(mclust)
  require(ggplot2)
  seed <- 42
  find_closest_value <- function(value, reference_vector){
    closest_index <- which.min(abs(reference_vector-value))
    closest_value <- reference_vector[closest_index]
    return(closest_value)
  }

  find_closest_x <- function(value, reference_vector) {
    closest_value <- reference_vector[which.min(abs(reference_vector -value))]
    return(closest_value)
  }
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
      fittings[[i]][[j]] <- normalmixEM(data[data[, i]>0, i], k = k_list[[i]][j], epsilon = epsilon, maxit = maxit, maxrestarts = maxrestarts)
    }
  }
  for (i in colnames(data)) {
    ggdensity[[i]] <- ggplot_build(ggplot(data[data[, i]>0, ], aes(x=data[data[, i]>0, i]))+
                                     geom_density(adjust=1))
    for (j in 1:length(k_list[[i]])) {
      if (k_list[[i]][j] == 2) {
        x <- fittings[[i]][[j]]$x
        dfs[[i]][[j]] <- data.frame(x = x,
                                    density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                                    density2 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
        #print(paste("Adding gg.x and gg.y for bimodal",i,  "DF-No:", which(names(dfs) == i)))
        dfs[[i]][[j]]$gg.x <- sapply(dfs[[i]][[j]]$x, find_closest_value, reference_vector = ggdensity[[i]]$data[[1]]$x)
        dfs[[i]][[j]]$gg.y <- ggdensity[[i]]$data[[1]]$y[match(dfs[[i]][[j]]$gg.x, ggdensity[[i]]$data[[1]]$x)]
        dfs.corrected[[i]][[j]] <- dfs[[i]][[j]]
        fittings.corrected[[i]][[j]] <- fittings[[i]][[j]]
        # Find the closes x value in gg.x
        x.mean <- find_closest_x(value = min(fittings[[i]][[1]]$mu), reference_vector = fittings[[i]][[1]]$x)
        C1.y <- unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$density1)
        C2.y <- unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$density2)

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
        }else {
          if (C1.y < (1-margin.den)*unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$gg.y)) {
            x <- fittings[[i]][[j]]$x
            sd1 <- fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
            cat(paste("Fixing Sigma for C1 in",i,"-Bimodal"))
            while (C1.y < unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$gg.y)) {
              sd1 <- sd1 - 0.005
              fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu == min(fittings.corrected[[i]][[j]]$mu))] <- sd1
              dfs.corrected[[i]][[j]] <- data.frame(x = x,
                                                    density1 = dnorm(x, mean = min(fittings.corrected[[i]][[j]]$mu), sd = fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu == min(fittings.corrected[[i]][[j]]$mu))]) * fittings.corrected[[i]][[j]]$lambda[which(fittings.corrected[[i]][[j]]$mu == min(fittings.corrected[[i]][[j]]$mu))],
                                                    density2 = dnorm(x, mean = max(fittings.corrected[[i]][[j]]$mu), sd = fittings.corrected[[i]][[j]]$sigma[which(fittings.corrected[[i]][[j]]$mu == max(fittings.corrected[[i]][[j]]$mu))]) * fittings.corrected[[i]][[j]]$lambda[which(fittings.corrected[[i]][[j]]$mu == max(fittings.corrected[[i]][[j]]$mu))])
              # Find the closest x value in gg.x
              C1.y <- unique(dfs.corrected[[i]][[j]][dfs.corrected[[i]][[j]]$x==x.mean,]$density1)
            }
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
          #print(paste("Adding gg.x and gg.y for trimodal",i,  "DF-No:", which(names(dfs) == i)))
          dfs[[i]][[j]]$gg.x <- sapply(dfs[[i]][[j]]$x, find_closest_value, reference_vector = ggdensity[[i]]$data[[1]]$x)
          dfs[[i]][[j]]$gg.y <- ggdensity[[i]]$data[[1]]$y[match(dfs[[i]][[j]]$gg.x, ggdensity[[i]]$data[[1]]$x)]
          dfs.corrected[[i]][[j]] <- dfs[[i]][[j]]
          fittings.corrected[[i]][[j]] <- fittings[[i]][[j]]
          # Find the closes x value in gg.x
          x.mean <- find_closest_x(value = min(fittings[[i]][[j]]$mu), reference_vector = fittings[[i]][[j]]$x)
          C1.y <- unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$density1)
          C2.y <- unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$density2)
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
          }else {
            if (C1.y < (1-margin.den)*unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$gg.y)) {
              x <- fittings[[i]][[j]]$x
              sd1 <- fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
              cat(paste("Fixing Sigma for C1 in",i,"-Trimodal"))
              while (C1.y < unique(dfs[[i]][[j]][dfs[[i]][[j]]$x==x.mean,]$gg.y)) {
                sd1 <- sd1-0.005
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
  }

  return(fittings.corrected)
}
