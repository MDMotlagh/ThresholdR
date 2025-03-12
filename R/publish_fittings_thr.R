#' Save the fitted plots for each possible k in surface marker distributions
#'
#' @param fittings A list of fitted parameters for each distribution (component) under each k across all markers.
#' @param k_list A list with set of possible k values for each marker as an element.
#' @param thresholds A data frame with calculated thresholds.
#' @param path A path to the directory in which the plots will be saved.
#' @returns A data frame with calculated thresholds: bi_thr1 defines the mean+3SD of the first component and bi_cut refers to the cutpoint on x axis where the two components meet in a bimodal (k=2) distribution. tri_thr1 and tri_thr2 indicate mean+3SD of the first and second components in a trimodal (k=3) distribution. Subsequently, tri_cut1 and tri_cut2 represent the cutpoints between first and second component, and second and third components, repectively, in a trimodal distribution. The user should choose which thereshold to use the recomended value is mean+3SD of the first component in k=2.
#' @examples
#' publish_fittings_thr(fittings = fittings, k_list = k.list, thresholds = thresholds, path = '/Users/oliaei/Projects/Packages/ThresholdR/Vignette/02Plots/', seed = 42)



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
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
                  geom_vline(xintercept = thresholds[i,"tri_thr1"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"tri_thr1"],y=+Inf,
                           label=paste("thr1=",
                                       round(thresholds[i,"tri_thr1"], digits =2)),vjust=2,geom="label")+
                  geom_vline(xintercept = thresholds[i,"tri_cut1"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"tri_cut1"],y=max(density(fittings[[i]][[j]]$x)$y)*0.8,
                           label=paste("cut1=",
                                       round(thresholds[i,"tri_cut1"], digits =2)),vjust=2,geom="label")+
                  geom_vline(xintercept = thresholds[i,"tri_thr2"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"tri_thr2"],y=max(density(fittings[[i]][[j]]$x)$y)*0.6,
                           label=paste("thr2=",
                                       round(thresholds[i,"tri_thr2"], digits =2)),vjust=2,geom="label")+
                  geom_vline(xintercept = thresholds[i,"tri_cut2"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"tri_cut2"],y=max(density(fittings[[i]][[j]]$x)$y)*0.4,
                           label=paste("cut2=",
                                       round(thresholds[i,"tri_cut2"], digits =2)),vjust=2,geom="label")
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
                          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
                    geom_vline(xintercept = thresholds[i,"bi_thr"],
                               linetype="dashed",
                               color = "black", linewidth=0.5)+
                    annotate(x= thresholds[i,"bi_thr"],y=+Inf,
                             label=paste("M1+3SD=",
                                         round(thresholds[i,"bi_thr"], digits =2)),vjust=2,geom="label")

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
                        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
                  geom_vline(xintercept = thresholds[i,"tri_thr1"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"tri_thr1"],y=+Inf,
                           label=paste("thr1=",
                                       round(thresholds[i,"tri_thr1"], digits =2)),vjust=2,geom="label")+
                  geom_vline(xintercept = thresholds[i,"tri_cut1"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"tri_cut1"],y=max(density(fittings[[i]][[j]]$x)$y)*0.8,
                           label=paste("cut1=",
                                       round(thresholds[i,"tri_cut1"], digits =2)),vjust=2,geom="label")+
                  geom_vline(xintercept = thresholds[i,"tri_thr2"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"tri_thr2"],y=max(density(fittings[[i]][[j]]$x)$y)*0.6,
                           label=paste("thr2=",
                                       round(thresholds[i,"tri_thr2"], digits =2)),vjust=2,geom="label")+
                  geom_vline(xintercept = thresholds[i,"tri_cut2"],
                             linetype="dashed",
                             color = "black", linewidth=0.5)+
                  annotate(x= thresholds[i,"tri_cut2"],y=max(density(fittings[[i]][[j]]$x)$y)*0.4,
                           label=paste("cut2=",
                                       round(thresholds[i,"tri_cut2"], digits =2)),vjust=2,geom="label")
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
                          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))+
                    geom_vline(xintercept = thresholds[i,"bi_thr"],
                               linetype="dashed",
                               color = "black", linewidth=0.5)+
                    annotate(x= thresholds[i,"bi_thr"],y=+Inf,
                             label=paste("M1+3SD=",
                                         round(thresholds[i,"bi_thr"], digits =2)),vjust=2,geom="label")+
                    geom_vline(xintercept = thresholds[i,"bi_cut"],
                               linetype="dashed",
                               color = "black", linewidth=0.5)+
                    annotate(x= thresholds[i,"bi_cut"],y=max(density(fittings[[i]][[j]]$x)$y)*0.5,
                             label=paste("cut=",
                                         round(thresholds[i,"bi_cut"], digits =2)),vjust=2,geom="label")
            )
            dev.off()
          }
        }

      }
    }

  }
}
