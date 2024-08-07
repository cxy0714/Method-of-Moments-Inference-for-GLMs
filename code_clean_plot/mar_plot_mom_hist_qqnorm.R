rm(list = ls())

library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)


# 定义三元组
triplets <- c("1_0","0_1","0_0")

# 遍历每个三元组
for (triplet in triplets) {
  file_path <- file.path("data/data_mar_mom_hist/", triplet)
  
  file_names <- list.files(file_path, pattern = "\\.Rda$", full.names = TRUE)
  
  # 提取n后的数字并找到最大的
  n_values <- sapply(file_names, function(x) as.numeric(sub(".*n_([0-9]+)_.*", "\\1", x)))
 
  max_n_index <- which.max(n_values)
  # 加载最大的n值对应的文件
    load(file_names[max_n_index])
  


  # 提取三元组的值
  triplet_parts <- strsplit(triplet, "_")[[1]]
  Is_sparse <- as.numeric(triplet_parts[1])
  Is_Rad <- as.numeric(triplet_parts[2])
  
  # 计算 p_over_n 并保留两位小数
  p_over_n <- round(p / n, 2)
  
  # 定义图表布局和标题
  pdf_filename_parameters <- paste0(
    "plots/plots_mar_hist_qqnorm/glm_mar_para_norm_sparse_",
    Is_sparse,
    "_Ra_",
    Is_Rad,
    "_p_over_n_",
    p_over_n,
    ".pdf"
  )


# 定义图表布局和标题
pdf(pdf_filename_parameters, height = 20, width = 16)
par(mfrow = c(5, 4))

# 标题
titles <- c(
  substitute(paste('E[A]')),
  substitute(paste('E[', X^T, ']', Sigma^-1, 'E[X]', sep = '')),
  substitute(paste('E[A ', X^T, ']', Sigma^-1, 'E[X]', sep = '')),
  substitute(paste('E[A ', X^T, ']', Sigma^-1, 'E[A X]', sep = '')),
  substitute(paste('E[A Y]')),
  substitute(paste('E[A Y ', X^T, ']', Sigma^-1, 'E[X]', sep = '')),
  substitute(paste('E[A Y ', X^T, ']', Sigma^-1, 'E[ X A]', sep = ''))
)

# 绘图
for (i in 1:6) {
  hist(moments_em[, i], breaks = 20, main = titles[i], xlab = '', freq = FALSE)
  lines(c( moments_truth[i], moments_truth[i]), c(0, 50), col = 'red', lwd = 2)
  qqnorm(moments_em[, i] - moments_truth[i], main = titles[i], pch = 20)
  qqline(moments_em[, i] - moments_truth[i], col = 'red', lwd = 2)
}

# 标题
titles <- c(
  substitute(paste(alpha^T, 'E[X]')),
  substitute(paste(alpha^T, Sigma, alpha, sep = '')),
  substitute(paste(alpha^T, Sigma, beta, sep = '')),
  substitute(paste('E[Y]'))
)

# 绘图
for (i in 1:4) {
  hist(sols_em[, i], breaks = 20, main = titles[i], xlab = '', freq = FALSE)
  lines(c(par_truth[i], par_truth[i]), c(0, 50), col = 'red', lwd = 2)
  qqnorm(sols_em[, i] - par_truth[i], main = titles[i], pch = 20)
  qqline(sols_em[, i] - par_truth[i], col = 'red', lwd = 2)
}
dev.off()


}



