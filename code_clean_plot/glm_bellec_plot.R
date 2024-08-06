rm(list = ls())

library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)

file_path <- "Data_glm_sparse_Rad"
file_names <- list.files(file_path, pattern = "\\.RData$", full.names = TRUE)

ns <- ps <- rep(NA,length(file_names))
alpha_bellec_p2_list <- list()
alpha_MoM_list <- list()
alpha_L2_bellec_list <- list()
alpha_L2_MoM_list <- list()
for (j in seq_along(file_names)) {
  load(file_names[j])
  alpha_bellec_p2_list[[j]] <- do.call(rbind, lapply(alpha_de_p2_T, function(x) x[1, ])) - alpha[1]
  alpha_MoM_list[[j]] <- alpha_est_N[,1]- alpha[1]
  alpha_L2_bellec_list[[j]] <- do.call(rbind,alpha_L2_est_B_p2)- par_truth[2]
  alpha_L2_MoM_list[[j]] <- alpha_L2_est_N - par_truth[2]
  ns[j] <- n
  ps[j] <- p
}
# 创建一个空的数据框来保存计算结果
results <- data.frame(n = integer(), lambda = numeric(), bias = numeric(), variance = numeric(), mse = numeric(), method = character())

# 计算alpha_bellec_p2_list的统计量
for (i in seq_along(alpha_bellec_p2_list)) {
  for (j in seq_along(lambda_value)) {
    simulations <- alpha_bellec_p2_list[[i]][, j]
    bias <- mean(simulations, na.rm = TRUE) 
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations )^2, na.rm = TRUE)  
    
    results <- rbind(results, data.frame(n = ns[i], lambda = lambda_value[j], bias = bias, variance = variance, mse = mse, method = "Bellec-Ridge"))
  }
}
# 计算alpha_bellec_p2_list的统计量
for (i in seq_along(alpha_bellec_p2_list)) {
  for (j in seq_along(lambda_value)) {
    simulations <- alpha_L2_bellec_list[[i]][, j]
    bias <- mean(simulations, na.rm = TRUE)
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations )^2, na.rm = TRUE) 
    
    results <- rbind(results, data.frame(n = ns[i], lambda = lambda_value[j], bias = bias, variance = variance, mse = mse, method = "L2-Bellec-Ridge"))
  }
}

# 计算alpha_MoM_list的统计量
for (i in seq_along(alpha_MoM_list)) {
  simulations <- alpha_MoM_list[[i]]
  bias <- mean(simulations, na.rm = TRUE) 
  variance <- var(simulations, na.rm = TRUE)
  mse <- mean((simulations)^2, na.rm = TRUE)
  
  results <- rbind(results, data.frame(n = ns[i], lambda = NA, bias = bias, variance = variance, mse = mse, method = "MoM"))
}
# 计算alpha_MoM_list的统计量
for (i in seq_along(alpha_MoM_list)) {
  simulations <- alpha_L2_MoM_list[[i]]
  bias <- mean(simulations, na.rm = TRUE) 
  variance <- var(simulations, na.rm = TRUE)
  mse <- mean((simulations )^2, na.rm = TRUE)
  
  results <- rbind(results, data.frame(n = ns[i], lambda = NA, bias = bias, variance = variance, mse = mse, method = "L2-MoM"))
}

# 按统计量转换数据框结构
results_long <- results %>%
  pivot_longer(cols = c(bias, variance, mse), names_to = "metric", values_to = "value")
results_long_mom <- results_long[results_long $method == "MoM",]
results_long_bellec <- results_long[results_long $method == "Bellec-Ridge",]

results_long_bellec_L2 <- results_long[results_long $method == "L2-Bellec-Ridge",]
results_long_mom_L2 <- results_long[results_long $method == "L2-MoM",]
library(gridExtra)  # For arranging multiple plots

# 数据已经是长格式
metrics <- c("bias" = "Bias", "variance" = "Variance", "mse" = "Mean Square Error")

# 准备所有图表的列表
plots <- list()
plot_index <- 1
for (met in names(metrics)) {
  extreme_data <- results_long_bellec %>%
    filter(metric == met) %>%
    filter(value > quantile(value, 0.90))  # 筛选90%分位数以上的值
  max_y <- results_long_bellec %>%
    filter(metric == met) %>%
    filter(value <= quantile(value, 0.90)) %>%  # 排除极端值
    summarise(max_value = max(value)) %>%
    pull(max_value)
  # 绘制Bellec-Ridge和MoM的图表
  plot <- ggplot() +
    geom_line(
      data = results_long_bellec %>% filter(metric == met),
      aes(
        x = log10(n),
        y = value,
        color = lambda,
        group = lambda
      ),
      alpha = 0.5
    ) +
    scale_color_gradient(name = expression(lambda),
                         low = "blue",
                         high = "red") +
    geom_point(
      data = results_long_mom %>% filter(metric == met),
      aes(x = log10(n), y = value),
      color = "black",
      shape = 21,
      size = 3,
      fill = "white"
    ) +
    geom_line(
      data = results_long_mom %>% filter(metric == met),
      aes(x = log10(n), y = value),
      color = "black",
      linetype = "dashed",
      size = 1
    ) +
    # geom_point(
    #   data = extreme_data,
    #   aes(x = log10(n), y = value),
    #   color = "red",
    #   shape = 8,  # 星形标记
    #   size = 4
    # ) +
    labs(
      title = expression(alpha[1]),
      x = "n",
      y = metrics[met]
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # 设置标题居中
      legend.position = "right"
    ) +
    scale_x_continuous(
      labels = function(x) 10^x,
      breaks = log10(c(1000, 2000, 4000))  # 设置显示特定的n值
    ) +
    coord_cartesian(ylim = c(NA, max_y))  # 限制y轴的范围
  
  # 存储每个指标的图表
  plots[[2 * plot_index - 1]] <- plot
  plot_index <- plot_index + 1
}

plot_index <- 1

for (met in names(metrics)) {
  extreme_data <- results_long_bellec_L2 %>%
    filter(metric == met) %>%
    filter(value > quantile(value, 0.90))  # 筛选90%分位数以上的值
  max_y <- results_long_bellec_L2 %>%
    filter(metric == met) %>%
    filter(value <= quantile(value, 0.90)) %>%  # 排除极端值
    summarise(max_value = max(value)) %>%
    pull(max_value)
  
  # 绘制Bellec-Ridge和MoM的图表
  plot <- ggplot() +
    geom_line(
      data = results_long_bellec_L2 %>% filter(metric == met),
      aes(
        x = log10(n),
        y = value,
        color = lambda,
        group = lambda
      ),
      alpha = 0.5
    ) +
    scale_color_gradient(name = expression(lambda),
                         low = "blue",
                         high = "red") +
    geom_point(
      data = results_long_mom_L2 %>% filter(metric == met),
      aes(x = log10(n), y = value),
      color = "black",
      shape = 21,
      size = 3,
      fill = "white"
    ) +
    geom_line(
      data = results_long_mom_L2 %>% filter(metric == met),
      aes(x = log10(n), y = value),
      color = "black",
      linetype = "dashed",
      size = 1
    )   + 
    # geom_point(
    #   data = extreme_data,
    #   aes(x = log10(n), y = value),
    #   color = "red",
    #   shape = 8,  # 星形标记
    #   size = 4
    # ) 
    labs(
      title = substitute(paste(alpha^T, Sigma, alpha, sep = '')),
      x = "n",
      y = metrics[met]
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # 设置标题居中
      legend.position = "right"
    ) +
    scale_x_continuous(
      labels = function(x) 10^x,
      breaks = log10(c(1000, 2000, 4000))  # 设置显示特定的n值
    ) +
    coord_cartesian(ylim = c(NA, max_y))  # 限制y轴的范围
  
  # 存储每个指标的图表
  plots[[2 * plot_index ]] <- plot
  plot_index <- plot_index + 1
}

# 合并图表为单列布局
combined_plot <- do.call(grid.arrange, c(plots, ncol = 2))
print(combined_plot)

# 打开PDF设备
pdf("glm_lambda_sparse_rademacher.pdf", width = 8.85, height = 6.6375)

# 使用grid.draw()来绘制gtable对象

grid.draw(combined_plot)

# 关闭PDF设备
dev.off()
