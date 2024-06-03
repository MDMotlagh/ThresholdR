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
testls[["up_up_pos"]] <- c(2, 5, 6)
testls[["up_down_pos"]] <- c(3, 2, 4)
testls[["down_up_pos"]] <- c(4, 2, 3)
testls[["down_down_pos"]] <- c(6, 5, 2)
testls[["up_up_neg"]] <- c(-6, -5, -2)
testls[["up_down_neg"]] <- c(-6, -5, -9)
testls[["down_up_neg"]] <- c(-5, -7, -1)
testls[["down_down_neg"]] <- c(-3, -5, -6)
testls[["up_up_mix"]] <- c(-2, 1, 3)
testls[["up_down_mix"]] <- c(1, 2, -1)
testls[["down_up_mix"]] <- c(-1, -2, 1)
testls[["down_down_mix"]] <- c(2, 1, -1)


