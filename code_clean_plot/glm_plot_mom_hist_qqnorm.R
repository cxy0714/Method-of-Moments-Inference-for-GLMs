rm(list = ls())

library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)

indices_first_1 <- 1:10
indices_middle_2 <- 100:109
indices_selected <- c(indices_first_1, indices_middle_2)
# 定义三元组
triplets <- c( "0_0_0", "0_0_1")

# c("0_0_0", "0_0_1", "1_0_0", "1_0_1", "1_1_1")
# 遍历每个三元组
for (triplet in triplets) {
  file_path <- file.path("data/data_glm_mom_after_cluster", triplet)
  
  file_names <- list.files(file_path, pattern = "\\.Rda$", full.names = TRUE)
  
  # 提取n后的数字并找到最大的
  n_values <- sapply(file_names, function(x) as.numeric(sub(".*_n_([0-9]+)_.*", "\\1", x)))
  
  max_n_index <- which.max(n_values)
  # 加载最大的n值对应的文件
  load(file_names[max_n_index])
  

  # 提取三元组的值
  triplet_parts <- strsplit(triplet, "_")[[1]]
  Is_sparse <- as.numeric(triplet_parts[1])
  Is_sparse_only_one <- as.numeric(triplet_parts[2])
  Is_Rad <- as.numeric(triplet_parts[3])
  
  # 计算 p_over_n 并保留两位小数
  p_over_n <- round(p / n, 2)
  # 定义图表布局和标题
  pdf_filename <- paste0(
    "plots/plots_glm_mom_hist/glm_mom_norm_sparse_",
    Is_sparse,
    "_one_",
    Is_sparse_only_one,
    "_Ra_",
    Is_Rad,
    "_p_over_n_",
    p_over_n,
    ".pdf"
  )
  
  pdf(pdf_filename, height = 20, width = 16)
  par(mfrow = c(4, 4))
  
  # 标题
  titles <- c(
    expression(paste('E[A]')),
    expression(paste('E[', X^T, ']', Sigma^-1, 'E[X]', sep = '')),
    expression(paste('E[A ', X^T, ']', Sigma^-1, 'E[X]', sep = '')),
    expression(paste('E[A ', X^T, ']', Sigma^-1, 'E[A X]', sep = ''))
  )

  # 绘图
  for (i in 1:4) {
    hist(moments_em[, i], breaks = 20, main = titles[i], xlab = '', freq = FALSE)
    lines(c(moments_truth[i], moments_truth[i]), c(0, 50), col = 'red', lwd = 2)
    qqnorm(moments_em[, i] - moments_truth[i], main = titles[i], pch = 20)
    qqline(moments_em[, i] - moments_truth[i], col = 'red', lwd = 2)
  }
  
  # 标题
  titles <- c(
    expression(paste(alpha^T, 'E[X]')),
    expression(paste(alpha^T, Sigma, alpha, sep = ''))
  )

  # 绘图
  for (i in 1:2) {
    hist(sols_em[, i], breaks = 20, main = titles[i], xlab = '', freq = FALSE)
    lines(c(par_truth[i], par_truth[i]), c(0, 50), col = 'red', lwd = 2)
    qqnorm(sols_em[, i] - par_truth[i], main = titles[i], pch = 20)
    qqline(sols_em[, i] - par_truth[i], col = 'red', lwd = 2)
  }
  
  # 标题
  titles <- c(
    expression(alpha[1] ~ "(nonzero)"),
    expression(alpha[100] ~ "(nonzero)")
  )
  
  alpha <- alpha[indices_selected]
  # 绘图
  for (i in c(1,2)) {
    hist(alpha_est_N[, (i-1)*10+ 1], breaks = 20, main = titles[i], xlab = '', freq = FALSE)
    lines(c(alpha[(i-1)*10+ 1], alpha[(i-1)*10+ 1]), c(0, 50), col = 'red', lwd = 2)
    qqnorm(alpha_est_N[, (i-1)*10+ 1] - alpha[(i-1)*10+ 1], main = titles[i], pch = 20)
    qqline(alpha_est_N[, (i-1)*10+ 1] - alpha[(i-1)*10+ 1], col = 'red', lwd = 2)
  }
  dev.off()
  
}