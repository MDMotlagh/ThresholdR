#' Publish Fittings with Assigned Thresholds
#'
#' This function publishes bimodal and trimodal fittings with assigned thresholds in the location of interest.
#'
#'
#' @param fittings a list of fittings for each protein, each containg fitting for all possible fits
#' @param k_list a list of vectors, each vector contains values for possible components fitted best between k=1 to k=k for normalized distribution of each protein
#' @param thresholds a dataframe of thresholds which has proteins or surface markers as rows and bimodal and/or trimodal thresholds as columns.
#'
#' @param path the path you would like to save the fitting density plots.
#'
#' @return It prints the the plots as pdf files in the path provided for path argument.
#'
#'
#' @examples
#' data("adt.df")
#' publish_fittings_thr(fittings = fittings.test, k_list =k_list.test ,thresholds = thresholds.test, path = "~/test_run/figures/")

#'

publish_fittings_thr <- function(fittings, k_list,thresholds, path=path) {
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
        pdf(paste0(path,"/", i, "_Trimodal.pdf"), width = 7, height = 5, unit="in")
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
