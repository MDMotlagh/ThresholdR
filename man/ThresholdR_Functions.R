#Functions for ThresholdR
calculateBICvalues <- function(data, k=3, parallel = FALSE) {
  require(mixtools)
  require(mclust)
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
        model <- Mclust(data[data[,i]>0, i], G = j)
        temp_bic_values[j] <- max(model$BIC[!is.na(model$BIC)])
      }
      return(temp_bic_values)
    }
  } else {
    foreach(i = colnames(data), .packages="mclust") %do% {
      print(i)
      temp_bic_values <- numeric(k)
      for (j in 1:k) {
        print(j)
        model <- Mclust(data[data[,i]>0, i], G = j)
        temp_bic_values[j] <- max(model$BIC[!is.na(model$BIC)])
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

#F2: Get the K_list####
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


#F3: Perform fittings####
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
        print(paste("Adding gg.x and gg.y for bimodal",i,  "DF-No:", which(names(dfs) == i)))
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
          print(paste("Adding gg.x and gg.y for trimodal",i,  "DF-No:", which(names(dfs) == i)))
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
#F4: Get thresholds####
get_thresholds <- function(fittings, k_list){
  thresholds <- data.frame(row.names = names(fittings),
                           bimodal_thr = rep(NA, length(names(fittings))),
                           trimodal_thr_C1 = rep(NA, length(names(fittings))),
                           trimodal_thr_C2 = rep(NA, length(names(fittings))))
  for (i in names(fittings)) {
    for (j in 1:length(k_list[[i]])) {
      if (k_list[[i]][j]==2) {
        thresholds[i,"bimodal_thr"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
      } else {
        if (k_list[[i]][j]==3) {
          thresholds[i,"trimodal_thr_C1"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
          thresholds[i,"trimodal_thr_C2"] <- fittings[[i]][[j]]$mu[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu)&fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))]+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu!=min(fittings[[i]][[j]]$mu)&fittings[[i]][[j]]$mu!=max(fittings[[i]][[j]]$mu))]
        }
      }
    }
  }
  return(thresholds)
}
#F5: Publish Fittings####
publish_fittings_thr <- function(fittings, k_list,thresholds, path=path, seed=42) {
  require(ggplot2)
  for (i in names(fittings)) {
    if (grepl("/", i)) {
      for (j in 1:length(k_list[[i]])) {
        if (k_list[[i]][j] == 3) {
          x <- fittings[[i]][[j]]$x
          df <- data.frame(x = x,
                           density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                           density2 = dnorm(x, mean = fittings[[i]][[j]]$mu[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                                            sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                           density3 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
          pdf(paste0(path,"/", gsub("/", "_", i), "_Trimodal.pdf"), width = 7, height = 5)
          print(ggplot(df, aes(x = x)) +
                  geom_density(aes(x=x,y=after_stat(density)), adjust=1) +
                  ggtitle(paste("Trimodal Fitting for", i)) +
                  geom_line(aes(y = density1, color = "C1")) +
                  geom_line(aes(y = density2, color = "C2")) +
                  geom_line(aes(y = density3, color = "C3")) +
                  labs(x = paste(i), y = "Density", colour="Components") +
                  scale_color_manual(values = c("C1" = "blue", "C2" = "red", "C3" = "purple")) +
                  theme(legend.text = element_text(size = 16),
                        legend.title = element_text(size = 18),
                        axis.text = element_text(size = 16),
                        axis.title = element_text(size = 20),
                        plot.title = element_text(size = 24),
                        panel.background = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill=NA, size=1))+
                  geom_vline(xintercept = thresholds[i,"trimodal_thr_C1"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"trimodal_thr_C1"],y=+Inf,
                           label=paste("C1=",
                                       round(thresholds[i,"trimodal_thr_C1"], digits =2)),vjust=2,geom="label")+
                  geom_vline(xintercept = thresholds[i,"trimodal_thr_C2"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"trimodal_thr_C2"],y=max(density(fittings.wc[[i]][[j]]$x)$y)*0.8,
                           label=paste("C2=",
                                       round(thresholds[i,"trimodal_thr_C2"], digits =2)),vjust=2,geom="label")
          )
          
          dev.off()
        } else {
          if (k_list[[i]][j] == 2) {
            x <- fittings[[i]][[j]]$x
            df <- data.frame(x = x,
                             density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                             density2 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
            pdf(paste0(path,"/", gsub("/", "_", i), "_Bimodal.pdf"), width = 7, height = 5 )
            print(ggplot(df, aes(x = x)) +
                    geom_density(aes(x=x,y=after_stat(density)), adjust=1) +
                    ggtitle(paste("Bimodal Fitting for", i)) +
                    geom_line(aes(y = density1, color = "C1")) +
                    geom_line(aes(y = density2, color = "C2")) +
                    labs(x = paste(i), y = "Density", colour="Components") +
                    scale_color_manual(values = c("C1" = "blue", "C2" = "red")) +
                    theme(legend.text = element_text(size = 16),
                          legend.title = element_text(size = 18),
                          axis.text = element_text(size = 16),
                          axis.title = element_text(size = 20),
                          plot.title = element_text(size = 24),
                          panel.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.line = element_line(colour = "black"),
                          panel.border = element_rect(colour = "black", fill=NA, size=1))+
                    geom_vline(xintercept = thresholds[i,"bimodal_thr"],
                               linetype="dashed",
                               color = "black", linewidth=0.5)+
                    annotate(x= thresholds[i,"bimodal_thr"],y=+Inf,
                             label=paste("M1+3SD=",
                                         round(thresholds[i,"bimodal_thr"], digits =2)),vjust=2,geom="label")
                  
            )
            dev.off()
          }
        }
        
      }
    }else {
      for (j in 1:length(k_list[[i]])) {
        if (k_list[[i]][j] == 3) {
          x <- fittings[[i]][[j]]$x
          df <- data.frame(x = x,
                           density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                           density2 = dnorm(x, mean = fittings[[i]][[j]]$mu[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                                            sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                           density3 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
          pdf(paste0(path,"/", i, "_Trimodal.pdf"), width = 7, height = 5)
          print(ggplot(df, aes(x = x)) +
                  geom_density(aes(x=x,y=after_stat(density)), adjust=1) +
                  ggtitle(paste("Trimodal Fitting for", i)) +
                  geom_line(aes(y = density1, color = "C1")) +
                  geom_line(aes(y = density2, color = "C2")) +
                  geom_line(aes(y = density3, color = "C3")) +
                  labs(x = paste(i), y = "Density", colour="Components") +
                  scale_color_manual(values = c("C1" = "blue", "C2" = "red", "C3" = "purple")) +
                  theme(legend.text = element_text(size = 16),
                        legend.title = element_text(size = 18),
                        axis.text = element_text(size = 16),
                        axis.title = element_text(size = 20),
                        plot.title = element_text(size = 24),
                        panel.background = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        axis.line = element_line(colour = "black"),
                        panel.border = element_rect(colour = "black", fill=NA, size=1))+
                  geom_vline(xintercept = thresholds[i,"trimodal_thr_C1"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"trimodal_thr_C1"],y=+Inf,
                           label=paste("C1=",
                                       round(thresholds[i,"trimodal_thr_C1"], digits =2)),vjust=2,geom="label")+
                  geom_vline(xintercept = thresholds[i,"trimodal_thr_C2"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"trimodal_thr_C2"],y=max(density(fittings.wc[[i]][[j]]$x)$y)*0.8,
                           label=paste("C2=",
                                       round(thresholds[i,"trimodal_thr_C2"], digits =2)),vjust=2,geom="label")
          )
          
          dev.off()
        } else {
          if (k_list[[i]][j] == 2) {
            x <- fittings[[i]][[j]]$x
            df <- data.frame(x = x,
                             density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                             density2 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
            pdf(paste0(path,"/", i, "_Bimodal.pdf"), width = 7, height = 5 )
            print(ggplot(df, aes(x = x)) +
                    geom_density(aes(x=x,y=after_stat(density)), adjust=1) +
                    ggtitle(paste("Bimodal Fitting for", i)) +
                    geom_line(aes(y = density1, color = "C1")) +
                    geom_line(aes(y = density2, color = "C2")) +
                    labs(x = paste(i), y = "Density", colour="Components") +
                    scale_color_manual(values = c("C1" = "blue", "C2" = "red")) +
                    theme(legend.text = element_text(size = 16),
                          legend.title = element_text(size = 18),
                          axis.text = element_text(size = 16),
                          axis.title = element_text(size = 20),
                          plot.title = element_text(size = 24),
                          panel.background = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.line = element_line(colour = "black"),
                          panel.border = element_rect(colour = "black", fill=NA, size=1))+
                    geom_vline(xintercept = thresholds[i,"bimodal_thr"],
                               linetype="dashed",
                               color = "black", linewidth=0.5)+
                    annotate(x= thresholds[i,"bimodal_thr"],y=+Inf,
                             label=paste("M1+3SD=",
                                         round(thresholds[i,"bimodal_thr"], digits =2)),vjust=2,geom="label")
                  
            )
            dev.off()
          }
        }
        
      }
    }
    
  }
}

