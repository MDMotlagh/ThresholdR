source("/Users/ataraskin/Git/DataAnalysisProjects/CITE-Seq/Melissa_Human_CITE-Seq_2023/Modified_pipeline/05_ThresholdeR/functions.R") 

testdf<- data.frame(matrix(NA, nrow = 12, ncol = 3))
for (n in 1:3) {
  colnames(testdf)[n] <- paste0("bic.", n, "-", n+1) 
}


rownames(testdf) <- c("up_up_pos",
                      "up_down_pos",
                      "down_up_pos",
                      "down_down_pos",
                      "up_up_neg",
                      "up_down_neg",
                      "down_up_neg",
                      "down_down_neg",
                      "up_up_mix",
                      "up_down_mix",
                      "down_up_mix",
                      "down_down_mix")

testdf["up_up_pos", ] <- c(2, 5, 6)
testdf["up_down_pos", ] <- c(3, 2, 4)
testdf["down_up_pos", ] <- c(4, 2, 3)
testdf["down_down_pos", ] <- c(6, 5, 2)
testdf["up_up_neg", ] <- c(-6, -5, -2)
testdf["up_down_neg", ] <- c(-6, -5, -9)
testdf["down_up_neg", ] <- c(-5, -7, -1)
testdf["down_down_neg", ] <- c(-3, -5, -6)
testdf["up_up_mix", ] <- c(-2, 1, 3)
testdf["up_down_mix", ] <- c(1, 2, -1)
testdf["down_up_mix", ] <- c(-1, -2, 1)
testdf["down_down_mix", ] <- c(2, 1, -1)

testls <- list()
testls[["up_up_pos"]] <- c(2, 5, 6, 10)
testls[["up_down_pos"]] <- c(3, 2, 4, 7)
testls[["down_up_pos"]] <- c(4, 2, 3, 6)
testls[["down_down_pos"]] <- c(6, 5, 2, 1)
testls[["up_up_neg"]] <- c(-6, -5, -2, -1)
testls[["up_down_neg"]] <- c(-6, -5, -9, -10)
testls[["down_up_neg"]] <- c(-5, -7, -4, -1)
testls[["down_down_neg"]] <- c(-3, -5, -6, -10)
testls[["up_up_mix"]] <- c(-2, 1, 3, 5)
testls[["up_down_mix"]] <- c(1, 2, -1, -3)
testls[["down_up_mix"]] <- c(-1, -2, 1, 3)
testls[["down_down_mix"]] <- c(2, 1, -1, -4)
testls[["zero"]] <- c(1, 4, 4, 5)
testls[["infinity"]] <- c(-4, -1, 0, 1)

list_names <- c("up_up_pos",
                   "up_down_pos",
                   "down_up_pos",
                   "down_down_pos",
                   "up_up_neg",
                   "up_down_neg",
                   "down_up_neg",
                   "down_down_neg",
                   "up_up_mix",
                   "up_down_mix",
                   "down_up_mix",
                   "down_down_mix",
                   "zero",
                   "infinity")

all_margins <- defineMargin(list_names, testls)
out2 <- GetAllPossibleFits(testls, margin_thr=0.1)

out3 <- lapply(out2, function(x) x[!is.na(x)])
out4 <- lapply(out3, function(x) x[1])

cleaned_list_of_lists <- lapply(list_of_lists, remove_na)

output <- data.frame(matrix(nrow = 1, ncol = 3))
colnames(output) <- paste0("bic.", 1:3, "-", 2:(3 + 1))



fittings <- fitting
names(fittings)
k_list

names(bic_values)
bic_values[1]

thresholds <- data.frame(row.names = names(fittings),
                         bimodal_thr = rep(NA, length(names(fittings))),
                         trimodal_thr = rep(NA, length(names(fittings))))
for (i in names(fittings)) {
  for (j in 1:length(k_list[[i]])) {
    if (k_list[[i]][j]==2) {
      thresholds[i,"bimodal_thr"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
    }else {
      if (k_list[[i]][j]==3) {
        thresholds[i,"trimodal_thr"] <- min(fittings[[i]][[j]]$mu)+3*fittings[[i]][[j]]$sigma[which(fittings[[i]][[j]]$mu == min(fittings[[i]][[j]]$mu))]
      }
    }
  }
}




