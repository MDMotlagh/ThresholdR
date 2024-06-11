# Function to calculate the bic_values
calculateBICvalues <- function(data, k=4, parallel = FALSE) {
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



defineMargin <- function(markers, bic_values) {
  max_columns <- max(sapply(bic_values, function(vec) length(vec) - 1))
  
  master_output <- lapply(markers, function(name) {
    vec <- bic_values[[name]]
    column_length <- length(vec) - 1
    output <- data.frame(matrix(nrow = 1, ncol = max_columns))
    colnames(output) <- paste0("bic.", 1:max_columns, "-", 2:(max_columns + 1))
    for (i in 1:column_length) {
      if (vec[i] == 0) {
        output[1, i] <- "BIC value is zero!"
      } else {
        output[1, i] <- round(abs((vec[i + 1] - vec[i]) / vec[i])*100, digits = 2)
        output[1, i] <- ifelse(vec[i]>vec[i + 1], (-1)*output[1, i], output[1, i])
      }
      
    }
    output$max_index <- which.max(vec)
    x <- output$max_index
    output$modality_based_on_BIC_max <- ifelse(x == 1, "Unimodal",
                                        ifelse(x == 2, "Bimodal",
                                        ifelse(x == 3, "Trimodal",
                                        ifelse(x == 4, "Quadrimodal",
                                                       "Multimodal"))))
    return(output)
  })
  
  master_output <- do.call(rbind, master_output)
  rownames(master_output) <- markers
  
  # Check that in all_margins there are no markers with unimodal distribution.
  check <- ifelse(all(master_output$max_index != 1), 
         "There are no main markers with unimodal distribution!",
         paste0("Discard those main markers and run again: ", 
                paste(rownames(master_output)[master_output$max_index == 1], collapse = ", ")))
  print(check)
    return(master_output)
}

# defineMargin <- function(markers, bic_values) {
#   bic_values_indices <- lapply(names(bic_values), function(name) {
#     vec <- bic_values[[name]]
#     column_length <- (length(vec)-1)
#     max_columns <<- max(column_length)
#     if (all(vec >= 0)) {
#       sapply(1:(length(vec) - 1), function(i) {
#         percent_diff <- (vec[i + 1] - vec[i]) / vec[i]
#         # print(paste(name, percent_diff))
#       })
#     } else {
#       if (all(vec <= 0)) {
#         sapply(1:(length(vec) - 1), function(i) {
#           percent_diff <- (vec[i] - vec[i + 1]) / vec[i]
#           # print(paste(name, percent_diff))
#         })
#       } else {
#         positive_values <- vec[vec >= 0]
#         if (length(positive_values)>1) {
#           sapply(1:(length(positive_values) - 1), function(i) {
#             percent_diff <- (positive_values[i + 1] - positive_values[i]) / positive_values[i]
#             # print(paste(name, percent_diff))
#           })
#         } else {
#           # print(paste(name, "Positive value", positive_values[1]))
#           # positive_values[1]
#         }
#       }
#       
#     }
#   })
#   output <- data.frame()
#   output[1,1:max_columns] <- NA
#   for (n in 1:max_columns) {
#     colnames(output)[n] <- paste0("bic.", n, "-", n+1) 
#   }
#   
#   print(colnames(output))
#   names(bic_values_indices) <- names(bic_values)
#   return(bic_values_indices)
# }


library(mclust)
# NA values removed!
GetAllPossibleFits <- function(bic_values, margin_thr) {
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
  # This line was added.
  k_list <- lapply(k_list, function(x) x[!is.na(x)])
  return(k_list)
}

fit_k <- function(data, 
                  k_list, 
                  margin.den=0.1,
                  epsilon = 0.01, 
                  maxit = 10000, 
                  maxrestarts = 5, 
                  seed) {
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

get_thresholds <- function(fittings, k_list){
  thresholds <- data.frame(row.names = names(fittings),
                           bimodal_thr = rep(NA, length(names(fittings))),
                           trimodal_thr = rep(NA, length(names(fittings))))
  for (i in names(fittings)) {
    for (j in 1:length(k_list[[i]])) {
      if (k_list[[i]][j]==2) {
        thresholds[i,"bimodal_thr"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
      } else {
        if (k_list[[i]][j]==3) {
          thresholds[i,"trimodal_thr"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
        }
      }
    }
  }
  return(thresholds)
}

publish_fittings_thr <- function(fittings, k_list,thresholds, path=path, seed) {
  require(ggplot2)
  for (i in names(fittings)) {
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
                geom_vline(xintercept = thresholds[i,"trimodal_thr"],
                           linetype="dashed",
                           color = "black", linewidth=0.5)+
                annotate(x= thresholds[i,"trimodal_thr"],y=+Inf,
                         label=paste("M1+3SD=",
                                     round(thresholds[i,"trimodal_thr"], digits =2)),vjust=2,geom="label")
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

modality_decision <- function(marker_name, thresholds_value, modality_value, output_table) {
  response <- readline(prompt = "Press [y] or [n] to accept or decline the fitting: ")
  # marker_name <- marker_name
  # thresholds_value <- thresholds_value
  modality_value <- modality_value
  if (response == "y") {
    m <- ifelse(modality_value == 3, "trimodal",
           ifelse(modality_value == 2, "bimodal", "error"))
    df_row <- data.frame(marker = marker_name,
                         threshold = thresholds_value,
                         modality = m
                         )
    output_table <- rbind(output_table, df_row) 
    return(output_table)
  } else if (response == "n") {
    message("Skipped!")
    return(output_table)
  } else {
    message("You need to press 'y' or 'n' to accept or decline the fitting.\nOthewise press 'Ctr+C' to exit.")
    res <- modality_decision(marker_name, thresholds_value, modality_value, output_table)
    return(res)
  }
}

look_fittings_thr <- function(fittings, k_list, thresholds, path = path, seed) {
  require(ggplot2)
  require(gridExtra)
  output_table <- data.frame(marker = c(),
                             threshold = c(),
                             modality = c())
  
  for (i in names(fittings)) {
    # Loop splits into two cases, where 
    if ((length(k_list[[i]]) > 1) & all(c(2,3) %in% k_list[[i]])) {
      for (j in 1:length(k_list[[i]])) {
        if (k_list[[i]][j] == 3) {
          x <- fittings[[i]][[j]]$x
          df <- data.frame(x = x,
                           density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), 
                                            sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                           density2 = dnorm(x, mean = fittings[[i]][[j]]$mu[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                                            sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                           density3 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), 
                                            sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
          
          plot1 <- ggplot(df, aes(x = x)) +
            geom_density(aes(x = x, y = after_stat(density)), adjust = 1) +
            ggtitle(paste("Trimodal Fitting for", i)) +
            geom_line(aes(y = density1, color = "C1")) +
            geom_line(aes(y = density2, color = "C2")) +
            geom_line(aes(y = density3, color = "C3")) +
            labs(x = paste(i), y = "Density", colour = "Components") +
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
                  panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
            geom_vline(xintercept = thresholds[i, "trimodal_thr"],
                       linetype = "dashed",
                       color = "black", linewidth = 0.5) +
            annotate(x = thresholds[i, "trimodal_thr"], y = +Inf,
                     label = paste("M1+3SD=",
                                   round(thresholds[i, "trimodal_thr"], digits = 2)), vjust = 2, geom = "label")
          
          print(plot1)
          print(paste0("Marker ", i, " (trimodal) ", " (#", which(names(fittings)==i), ") out of ", length(names(fittings))))
          output_table <- modality_decision(marker_name=i, 
                                            thresholds_value=thresholds[i, "trimodal_thr"], 
                                            modality_value=k_list[[i]][j], 
                                            output_table)
          
        } else if (k_list[[i]][j] == 2) {
          x <- fittings[[i]][[j]]$x
          df <- data.frame(x = x,
                           density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), 
                                            sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                           density2 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), 
                                            sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
          
          plot2 <- ggplot(df, aes(x = x)) +
            geom_density(aes(x = x, y = after_stat(density)), adjust = 1) +
            ggtitle(paste("Bimodal Fitting for", i)) +
            geom_line(aes(y = density1, color = "C1")) +
            geom_line(aes(y = density2, color = "C2")) +
            labs(x = paste(i), y = "Density", colour = "Components") +
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
                  panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
            geom_vline(xintercept = thresholds[i, "bimodal_thr"],
                       linetype = "dashed",
                       color = "black", linewidth = 0.5) +
            annotate(x = thresholds[i, "bimodal_thr"], y = +Inf,
                     label = paste("M1+3SD=",
                                   round(thresholds[i, "bimodal_thr"], digits = 2)), vjust = 2, geom = "label")
          
          print(plot2)
          print(paste0("Marker ", i, " (bimodal) ", " (#", which(names(fittings)==i), ") out of ", length(names(fittings))))
          output_table <- modality_decision(marker_name=i, 
                                            thresholds_value=thresholds[i, "bimodal_thr"], 
                                            modality_value=k_list[[i]][j], 
                                            output_table)
          
        }
      }
    } else {
    for (j in 1:length(k_list[[i]])) {
      if (k_list[[i]][j] == 3) {
        x <- fittings[[i]][[j]]$x
        df <- data.frame(x = x,
                         density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), 
                                          sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                         density2 = dnorm(x, mean = fittings[[i]][[j]]$mu[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                                          sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu != min(fittings[[i]][[j]]$mu) & fittings[[i]][[j]]$mu != max(fittings[[i]][[j]]$mu))],
                         density3 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), 
                                          sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
        
        plot1 <- ggplot(df, aes(x = x)) +
          geom_density(aes(x = x, y = after_stat(density)), adjust = 1) +
          ggtitle(paste("Trimodal Fitting for", i)) +
          geom_line(aes(y = density1, color = "C1")) +
          geom_line(aes(y = density2, color = "C2")) +
          geom_line(aes(y = density3, color = "C3")) +
          labs(x = paste(i), y = "Density", colour = "Components") +
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
                panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
          geom_vline(xintercept = thresholds[i, "trimodal_thr"],
                     linetype = "dashed",
                     color = "black", linewidth = 0.5) +
          annotate(x = thresholds[i, "trimodal_thr"], y = +Inf,
                   label = paste("M1+3SD=",
                                 round(thresholds[i, "trimodal_thr"], digits = 2)), vjust = 2, geom = "label")
        
        print(plot1)
        print(paste0("Marker ", i, " (trimodal) ", " (#", which(names(fittings)==i), ") out of ", length(names(fittings))))
        output_table <- modality_decision(marker_name=i, 
                                thresholds_value=thresholds[i, "trimodal_thr"], 
                                modality_value=k_list[[i]][j], 
                                output_table)
        
      } else if (k_list[[i]][j] == 2) {
        x <- fittings[[i]][[j]]$x
        df <- data.frame(x = x,
                         density1 = dnorm(x, mean = min(fittings[[i]][[j]]$mu), 
                                          sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))],
                         density2 = dnorm(x, mean = max(fittings[[i]][[j]]$mu), 
                                          sd = fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))]) * fittings[[i]][[j]]$lambda[which(fittings[[i]][[j]]$mu == max(fittings[[i]][[j]]$mu))])
        
        plot2 <- ggplot(df, aes(x = x)) +
          geom_density(aes(x = x, y = after_stat(density)), adjust = 1) +
          ggtitle(paste("Bimodal Fitting for", i)) +
          geom_line(aes(y = density1, color = "C1")) +
          geom_line(aes(y = density2, color = "C2")) +
          labs(x = paste(i), y = "Density", colour = "Components") +
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
                panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
          geom_vline(xintercept = thresholds[i, "bimodal_thr"],
                     linetype = "dashed",
                     color = "black", linewidth = 0.5) +
          annotate(x = thresholds[i, "bimodal_thr"], y = +Inf,
                   label = paste("M1+3SD=",
                                 round(thresholds[i, "bimodal_thr"], digits = 2)), vjust = 2, geom = "label")
        
        print(plot2)
        print(paste0("Marker ", i, " (bimodal) ", " (#", which(names(fittings)==i), ") out of ", length(names(fittings))))
        output_table <- modality_decision(marker_name=i, 
                                          thresholds_value=thresholds[i, "bimodal_thr"], 
                                          modality_value=k_list[[i]][j], 
                                          output_table)
        
      }
    }}
  }
  return(output_table)
}



