rm(list = ls())

library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)

# 定义三元组
triplets <- c("0_0_0", "0_0_1", "1_0_0", "1_0_1", "1_1_1")
# c("0_0_0", "0_0_1", "1_0_0", "1_0_1", "1_1_1")
# 遍历每个三元组
for (triplet in triplets) {
  file_path <- file.path("data_norm_1.2_clean_2", triplet)
  
  file_names <- list.files(file_path, pattern = "\\.Rda$", full.names = TRUE)
  
  ns <- ps <- rep(0, length(file_names))
  alpha_bellec_p2_list_nonzero <- list()
  alpha_MoM_list_nonzero <- list()
  alpha_bellec_p2_list_zero <- list()
  alpha_MoM_list_zero <- list()
  alpha_L2_bellec_list <- list()
  alpha_L2_MoM_list <- list()
  
  for (j in seq_along(file_names)) {
    load(file_names[j])
    alpha_MoM_list_zero[[j]] <- alpha_est_N_total[, 11] - alpha[100]
    alpha_MoM_list_nonzero[[j]] <- alpha_est_N_total[, 1] - alpha[1]
    alpha_L2_MoM_list[[j]] <- alpha_L2_est_N_total - par_truth[2]
    ns[j] <- n_value
    ps[j] <- p_value
  }
  
  # 创建一个空的数据框来保存计算结果
  results <- data.frame(n = integer(), lambda = numeric(), bias = numeric(), variance = numeric(), mse = numeric(), method = character())
  
  # 计算alpha_MoM_list的统计量
  for (i in seq_along(alpha_MoM_list_zero)) {
    simulations <- alpha_MoM_list_zero[[i]]
    bias <- mean(simulations, na.rm = TRUE) 
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations)^2, na.rm = TRUE)
    
    results <- rbind(results, data.frame(n = ns[i], lambda = NA, bias = bias, variance = variance, mse = mse, method = "MoM-zero"))
  }
  for (i in seq_along(alpha_MoM_list_nonzero)) {
    simulations <- alpha_MoM_list_nonzero[[i]]
    bias <- mean(simulations, na.rm = TRUE) 
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations)^2, na.rm = TRUE)
    
    results <- rbind(results, data.frame(n = ns[i], lambda = NA, bias = bias, variance = variance, mse = mse, method = "MoM-nonzero"))
  }
  # 计算alpha_MoM_list的统计量
  for (i in seq_along(alpha_MoM_list_zero)) {
    simulations <- alpha_L2_MoM_list[[i]]
    bias <- mean(simulations, na.rm = TRUE) 
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations )^2, na.rm = TRUE)
    
    results <- rbind(results, data.frame(n = ns[i], lambda = NA, bias = bias, variance = variance, mse = mse, method = "L2-MoM"))
  }
  
  # 按统计量转换数据框结构
  results_long <- results %>%
    pivot_longer(cols = c(bias, variance, mse), names_to = "metric", values_to = "value")
  results_long_mom_zero <- results_long[results_long$method == "MoM-zero",]
  results_long_mom_nonzero <- results_long[results_long$method == "MoM-nonzero",]
  results_long_mom_L2 <- results_long[results_long$method == "L2-MoM",]
  
  library(gridExtra)  # For arranging multiple plots
  
  # 数据已经是长格式
  metrics <- c("bias" = "Bias", "variance" = "Variance", "mse" = "Mean Square Error")
  
  # 准备所有图表的列表
  plots <- list()
  plot_index <- 1
  
  for (met in names(metrics)) {
    max_y <- results_long_mom_nonzero %>%
      filter(metric == met) %>%
      filter(value <= quantile(value, 1)) %>%  # 排除极端值
      summarise(max_value = max(value)) %>%
      pull(max_value)
    
    # 绘制MoM的图表
    plot <- ggplot() +
      geom_point(
        data = results_long_mom_nonzero %>% filter(metric == met),
        aes(x = n, y = abs(value)),
        color = "black",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_nonzero %>% filter(metric == met),
        aes(x = n, y = abs(value)),
        color = "black",
        linetype = "dashed",
        size = 1
      ) +
      labs(
        title = expression(alpha[1] ~ "(nonzero)"),
        x = "n",
        y = metrics[met]
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),  # 设置标题居中
        legend.position = "right"
      ) +
      scale_x_continuous(
        breaks = ns  # 设置显示特定的n值
      ) 
    
    # 存储每个指标的图表
    plots[[3 * plot_index - 2]] <- plot
    plot_index <- plot_index + 1
  }
  
  plot_index <- 1
  
  for (met in names(metrics)) {
    max_y <- results_long_mom_zero %>%
      filter(metric == met) %>%
      filter(value <= quantile(value, 1)) %>%  # 排除极端值
      summarise(max_value = max(value)) %>%
      pull(max_value)
    
    # 绘制MoM的图表
    plot <- ggplot() +
      geom_point(
        data = results_long_mom_zero %>% filter(metric == met),
        aes(x = n, y = abs(value)),
        color = "black",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_zero %>% filter(metric == met),
        aes(x = n, y = abs(value)),
        color = "black",
        linetype = "dashed",
        size = 1
      ) +
      labs(
        title = expression(alpha[100] ~ "(zero)" ),
        x = "n",
        y = metrics[met]
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),  # 设置标题居中
        legend.position = "right"
      ) +
      scale_x_continuous(
        breaks = ns  # 设置显示特定的n值
      )
    # 存储每个指标的图表
    plots[[3 * plot_index - 1]] <- plot
    plot_index <- plot_index + 1
  }
  
  plot_index <- 1
  for (met in names(metrics)) {
    max_y <- results_long_mom_L2 %>%
      filter(metric == met) %>%
      filter(value <= quantile(value, 1)) %>%  # 排除极端值
      summarise(max_value = max(value)) %>%
      pull(max_value)
    
    # 绘制MoM的图表
    plot <- ggplot() +
      geom_point(
        data = results_long_mom_L2 %>% filter(metric == met),
        aes(x = n, y = abs(value)),
        color = "black",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_L2 %>% filter(metric == met),
        aes(x = n, y = abs(value)),
        color = "black",
        linetype = "dashed",
        size = 1
      ) +
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
        breaks = ns  # 设置显示特定的n值
      ) 
    # 存储每个指标的图表
    plots[[3 * plot_index]] <- plot
    plot_index <- plot_index + 1
  }
  
  # 合并图表为单列布局
  combined_plot <- do.call(grid.arrange, c(plots, ncol = 3))
  print(combined_plot)
  
  # 提取三元组的值
  triplet_parts <- strsplit(triplet, "_")[[1]]
  Is_sparse <- as.numeric(triplet_parts[1])
  Is_sparse_only_one <- as.numeric(triplet_parts[2])
  Is_Rad <- as.numeric(triplet_parts[3])
  
  # 计算 p_over_n 并保留两位小数
  p_over_n <- round(p_value / n_value, 2)
  
  # 生成 PDF 文件名
  pdf_filename <- paste0(
    "glm_mom_sparse_",
    Is_sparse,
    "_one_",
    Is_sparse_only_one,
    "_Ra_",
    Is_Rad,
    "_p_over_n_",
    p_over_n,
    ".pdf"
  )
  
  # 打开PDF设备
  pdf(pdf_filename, width = 8.85, height = 6.6375)
  grid.draw(combined_plot)
  # 关闭PDF设备
  dev.off()
  
}