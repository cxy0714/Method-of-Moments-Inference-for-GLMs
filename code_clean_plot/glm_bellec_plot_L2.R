rm(list = ls())

library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)

# c("0_0_0", "0_0_1", "1_0_0", "1_0_1", "1_1_1")
# 定义三元组
triplets <- c("1_0_1")
file_path_be <- file.path("data_bellec_1.2_has_L2_clean_100")
file_path_mom <- file.path("data_norm_1.2_clean_2")



indices_first_1 <- 1:10
indices_middle_2 <- 100:109
indices_selected <- c(indices_first_1, indices_middle_2)
lambda_value = exp(seq(log(10) , log(.05) , length.out = 12))

for (triplet in triplets) {
  # bellec -----
  # 将三元组拆分为单独的元素
  triplet_parts <- strsplit(triplet, "_")[[1]]
  
  # 构建文件名模式
  triplet_pattern <- paste0("sparse_",
                            triplet_parts[1],
                            "_one_",
                            triplet_parts[2],
                            "_Rad_",
                            triplet_parts[3])
  
  file_path_be <- file.path("data_bellec_1.2_has_L2_clean_100")
  # 获取所有文件名
  file_names <- list.files(file_path_be, pattern = "\\.Rda$", full.names = TRUE)
  
  # 获取当前三元组的所有文件
  current_files_be <- file_names[grepl(triplet_pattern, basename(file_names))]
  load(current_files_be)
  
  N.simu <- length(bellec_data_total)
  N.lambda <- length(lambda_value)
  
  ns_alpha <- ps_alpha <-   lambda_alpha <-  rep(NA, 2 * N.simu * N.lambda)
  Z_test_p2_TT_Ge <- alpha_de_N_p2_TT_Ge <- alpha_de_p2_TT_Ge <-
    matrix(
      NA,
      nrow = 2 * N.simu * N.lambda,
      ncol = nrow(bellec_data_total[[1]]$alpha_de_p2_TT_Ge[[1]])
    )
  
  ns_alpha_L2 <- ps_alpha_L2 <-   lambda_alpha_L2 <-  rep(0, N.simu * N.lambda)
  alpha_L2_de_N_Ge <- alpha_L2_de_Ge <-  rep(NA,  N.simu * N.lambda)
  
  alpha_de_N_p2_TT_Ge_nonzero <- rep(NA, 2 * N.simu * N.lambda)
  alpha_de_N_p2_TT_Ge_zero <- rep(NA, 2 * N.simu * N.lambda)
  
  alpha_de_p2_TT_Ge_nonzero <- rep(NA, 2 * N.simu * N.lambda)
  alpha_de_p2_TT_Ge_zero <- rep(NA, 2 * N.simu * N.lambda)
  
  for (i in c(1:N.simu)) {
    simulations <- bellec_data_total[[i]]
    alpha_L2_true <- sum(simulations$alpha *simulations$alpha)
    alpha_term <- simulations$alpha[indices_selected]

    for (j in c(1:N.lambda)) {
      ns_alpha_L2[(i - 1) * N.lambda + j] <- simulations$n
      ps_alpha_L2[(i - 1) * N.lambda + j] <- simulations$p
      lambda_alpha_L2[(i - 1) * N.lambda + j] <- lambda_value[j]
  
      alpha_L2_de_Ge[(i - 1) * N.lambda + j] <- simulations$alpha_L2_de_Ge[j] - alpha_L2_true
      alpha_L2_de_N_Ge[(i - 1) * N.lambda + j] <- simulations$alpha_L2_de_N_Ge[j] - alpha_L2_true
      
      for (l in c(1:2)) {
        ns_alpha[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l] <- simulations$n
        ps_alpha[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l] <- simulations$p
        lambda_alpha[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l] <- lambda_value[j]
        
        Z_test_p2_TT_Ge[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l, ] <- simulations$Z_test_p2_TT_Ge[[l]][, j]
        alpha_de_N_p2_TT_Ge[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l, ] <- simulations$alpha_de_N_p2_TT_Ge[[l]][, j]
        alpha_de_p2_TT_Ge[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l, ] <- simulations$alpha_de_p2_TT_Ge[[l]][, j]
        
        alpha_de_N_p2_TT_Ge_nonzero[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l] <- alpha_de_N_p2_TT_Ge[(i -
                                                                                                        1) * 2 * N.lambda + (j - 1) * 2 + l, 1] - alpha_term[1]
        alpha_de_p2_TT_Ge_nonzero[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l] <- alpha_de_p2_TT_Ge[(i -
                                                                                                    1) * 2 * N.lambda + (j - 1) * 2 + l, 1] - alpha_term[1]
        
        alpha_de_N_p2_TT_Ge_zero[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l] <- alpha_de_N_p2_TT_Ge[(i -
                                                                                                     1) * 2 * N.lambda + (j - 1) * 2 + l, 11] - alpha_term[11]
        alpha_de_p2_TT_Ge_zero[(i - 1) * 2 * N.lambda + (j - 1) * 2 + l] <- alpha_de_p2_TT_Ge[(i -
                                                                                                 1) * 2 * N.lambda + (j - 1) * 2 + l, 11] - alpha_term[11]
      }
      
    }
  }
  table_bellec <- data.frame(
    ns_alpha,
    lambda_alpha,
    alpha_de_p2_TT_Ge_nonzero,
    alpha_de_p2_TT_Ge_nonzero,
    alpha_de_N_p2_TT_Ge_zero,
    alpha_de_N_p2_TT_Ge_zero
  )
  
  # 定义每个 n 对应的样本数
  sample_size <- min(sum(ns_alpha == 1000 )/N.lambda , sum(ns_alpha == 5000 )/N.lambda)
  # 对每个 n 进行抽样
  table_bellec <- table_bellec %>%
    group_by(ns_alpha, lambda_alpha) %>%
    sample_n(size = sample_size, replace = FALSE) %>%
    ungroup()
  
  table_bellec_L2 <- data.frame(
    ns_alpha_L2,
    lambda_alpha_L2,
    alpha_L2_de_N_Ge,
    alpha_L2_de_Ge
  )
  sample_size <- min(sum(ns_alpha_L2 == 1000 )/N.lambda , sum(ns_alpha_L2 == 5000 )/N.lambda)
  # 对每个 n 进行抽样
  table_bellec_12 <- table_bellec_L2 %>%
    group_by(ns_alpha_L2, lambda_alpha_L2) %>%
    sample_n(size = sample_size, replace = FALSE) %>%
    ungroup()
  
  # 定义方法名称
  methods <- c(
    "alpha_de_p2_TT_Ge_nonzero",
    "alpha_de_p2_TT_Ge_nonzero.1",
    "alpha_de_N_p2_TT_Ge_zero",
    "alpha_de_p2_TT_Ge_zero"
  )
  
  # 创建一个空的数据框来保存计算结果
  results_bellec <- data.frame(
    n = integer(),
    lambda = numeric(),
    bias = numeric(),
    variance = numeric(),
    mse = numeric(),
    method = character()
  )
  
  # 遍历每个方法，计算并填充数据框
  for (method in methods) {
    temp_df <- table_bellec %>%
      group_by(ns_alpha, lambda_alpha) %>%
      summarise(
        bias = mean(!!sym(method)),
        variance = var(!!sym(method)),
        mse = mean((!!sym(method)) ^ 2),
        method = method,
        .groups = 'drop'
      )
    results_bellec <- rbind(results_bellec, temp_df)
  }
  # 定义方法名称
  methods <- c(
    "alpha_L2_de_N_Ge",
    "alpha_L2_de_Ge"
  )
  
  # 创建一个空的数据框来保存计算结果
  results_bellec_L2 <- data.frame(
    n = integer(),
    lambda = numeric(),
    bias = numeric(),
    variance = numeric(),
    mse = numeric(),
    method = character()
  )
  
  # 遍历每个方法，计算并填充数据框
  for (method in methods) {
    temp_df <- table_bellec_L2 %>%
      group_by(ns_alpha_L2, lambda_alpha_L2) %>%
      summarise(
        bias = mean(!!sym(method)),
        variance = var(!!sym(method)),
        mse = mean((!!sym(method)) ^ 2),
        method = method,
        .groups = 'drop'
      )
    results_bellec_L2 <- rbind(results_bellec_L2, temp_df)
  }
  # mom -----
  
  file_path <-  file.path(file_path_mom, triplet)
  file_names <- list.files(file_path, pattern = "\\.Rda$", full.names = TRUE)
  
  ns <- ps <- rep(0, length(file_names))
  alpha_MoM_list_nonzero <- list()
  alpha_MoM_list_zero <- list()
  alpha_L2_MoM_list <- list()
  
  for (j in seq_along(file_names)) {
    load(file_names[j])

    alpha_MoM_list_zero[[j]] <- alpha_est_N_total[, 11] - alpha[indices_selected[11]]
    alpha_MoM_list_nonzero[[j]] <- alpha_est_N_total[, 1] - alpha[1]
    alpha_L2_MoM_list[[j]] <- alpha_L2_est_N_total - par_truth[2]
    ns[j] <- n_value
    ps[j] <- p_value
  }
  results_mom <- data.frame(
    n = integer(),
    lambda = numeric(),
    bias = numeric(),
    variance = numeric(),
    mse = numeric(),
    method = character()
  )
  
  # 计算alpha_MoM_list的统计量
  for (i in seq_along(alpha_MoM_list_zero)) {
    simulations <- alpha_MoM_list_zero[[i]]
    bias <- mean(simulations, na.rm = TRUE)
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations) ^ 2, na.rm = TRUE)
    
    results_mom <- rbind(
      results_mom,
      data.frame(
        n = ns[i],
        lambda = NA,
        bias = bias,
        variance = variance,
        mse = mse,
        method = "MoM-zero"
      )
    )
  }
  for (i in seq_along(alpha_MoM_list_nonzero)) {
    simulations <- alpha_MoM_list_nonzero[[i]]
    bias <- mean(simulations, na.rm = TRUE)
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations) ^ 2, na.rm = TRUE)
    
    results_mom <- rbind(
      results_mom,
      data.frame(
        n = ns[i],
        lambda = NA,
        bias = bias,
        variance = variance,
        mse = mse,
        method = "MoM-nonzero"
      )
    )
  }
  # 计算alpha_MoM_list的统计量
  for (i in seq_along(alpha_MoM_list_zero)) {
    simulations <- alpha_L2_MoM_list[[i]]
    bias <- mean(simulations, na.rm = TRUE)
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations) ^ 2, na.rm = TRUE)
    
    results_mom <- rbind(
      results_mom,
      data.frame(
        n = ns[i],
        lambda = NA,
        bias = bias,
        variance = variance,
        mse = mse,
        method = "L2-MoM"
      )
    )
  }
  # final data ----

  # 按统计量转换数据框结构
  results_long_mom <- results_mom %>%
    pivot_longer(
      cols = c(bias, variance, mse),
      names_to = "metric",
      values_to = "value"
    )
  results_long_bellec <- results_bellec %>%
    pivot_longer(
      cols = c(bias, variance, mse),
      names_to = "metric",
      values_to = "value"
    )
  names(results_long_bellec) = names(results_long_mom)
  results_long_bellec_L2 <- results_bellec_L2 %>%
    pivot_longer(
      cols = c(bias, variance, mse),
      names_to = "metric",
      values_to = "value"
    )
  names(results_long_bellec_L2) = names(results_long_mom)
  
  
  # Root n
  results_long_bellec <- results_long_bellec %>%
    mutate(value = if_else(metric == "bias", value * sqrt(n), value))
  
  results_long_bellec_L2 <- results_long_bellec_L2 %>%
    mutate(value = if_else(metric == "bias", value * sqrt(n), value))
  
  results_long_mom <- results_long_mom %>%
    mutate(value = if_else(metric == "bias", value * sqrt(n), value))
  
  
  results_long_mom_zero <- results_long_mom[results_long_mom$method == "MoM-zero", ]
  results_long_mom_nonzero <- results_long_mom[results_long_mom$method == "MoM-nonzero", ]
  results_long_mom_L2 <- results_long_mom[results_long_mom$method == "L2-MoM", ]
  
  results_long_bellec_zero <- results_long_bellec[results_long_bellec$method == "alpha_de_p2_TT_Ge_zero", ]
  results_long_bellec_nonzero <- results_long_bellec[results_long_bellec$method == "alpha_de_p2_TT_Ge_nonzero", ]

  # plot -----
  
  library(gridExtra)  # For arranging multiple plots
  
  # 数据已经是长格式
  metrics <- c("bias" = "Root n Bias",
               "variance" = "Variance",
               "mse" = "Mean Square Error")
  
  # 准备所有图表的列表
  plots <- list()
  plot_index <- 1
  for (met in names(metrics)) {
    extreme_data <- results_long_bellec_nonzero %>%
      filter(metric == met) %>%
      filter(value > quantile(value, 0.95))  # 筛选90%分位数以上的值
    max_y <- results_long_bellec_nonzero %>%
      filter(metric == met) %>%
      filter(value <= quantile(value, 0.95)) %>%  # 排除极端值
      summarise(max_value = max(value)) %>%
      pull(max_value)
    max_y <- max(max_y,results_long_mom_nonzero %>%
                   filter(metric == met) %>%summarise(max_value = max(value)) %>%
                   pull(max_value) )
    # 绘制Bellec-Ridge和MoM的图表
    plot <- ggplot() +
      geom_line(
        data = results_long_bellec_nonzero %>% filter(metric == met),
        aes(
          x = n,
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
        data = results_long_mom_nonzero %>% filter(metric == met),
        aes(x = n, y = value),
        color = "black",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_nonzero %>% filter(metric == met),
        aes(x = n, y = value),
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
      theme(plot.title = element_text(hjust = 0.5),
            # 设置标题居中
            legend.position = "right") +
      scale_x_continuous(
        breaks = ns  # 设置显示特定的n值
      ) +
      coord_cartesian(ylim = c(NA, max_y))  # 限制y轴的范围
    
    # 存储每个指标的图表
    plots[[3 * plot_index - 2]] <- plot
    plot_index <- plot_index + 1
  }
  
  plot_index <- 1
  
  for (met in names(metrics)) {
    extreme_data <- results_long_bellec_zero %>%
      filter(metric == met) %>%
      filter(value > quantile(value, 0.95))  # 筛选90%分位数以上的值
    max_y <- results_long_bellec_zero %>%
      filter(metric == met) %>%
      filter(value <= quantile(value, 0.95)) %>%  # 排除极端值
      summarise(max_value = max(value)) %>%
      pull(max_value)
    max_y <- max(max_y,results_long_mom_zero %>%
                   filter(metric == met) %>%summarise(max_value = max(value)) %>%
                   pull(max_value) )
    # 绘制Bellec-Ridge和MoM的图表
    plot <- ggplot() +
      geom_line(
        data = results_long_bellec_zero %>% filter(metric == met),
        aes(
          x = n,
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
        data = results_long_mom_zero %>% filter(metric == met),
        aes(x = n, y = value),
        color = "black",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_zero %>% filter(metric == met),
        aes(x = n, y = value),
        color = "black",
        linetype = "dashed",
        size = 1
      )  +
      labs(
        title = expression(alpha[100] ~ "(nonzero)"),
        x = "n",
        y = metrics[met]
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            # 设置标题居中
            legend.position = "right") +
      scale_x_continuous(
        breaks = ns  # 设置显示特定的n值
      ) +
      coord_cartesian(ylim = c(NA, max_y))  # 限制y轴的范围
    
    # 存储每个指标的图表
    plots[[3 * plot_index - 1 ]] <- plot
    plot_index <- plot_index + 1
  }
  
  
  ## plot L2 ----
  plot_index <- 1
  
  for (met in names(metrics)) {
    extreme_data <- results_long_bellec_L2 %>%
      filter(metric == met) %>%
      filter(value > quantile(value, 0.95))  # 筛选90%分位数以上的值
    max_y <- results_long_bellec_L2 %>%
      filter(metric == met) %>%
      filter(value <= quantile(value, 0.95)) %>%  # 排除极端值
      summarise(max_value = max(value)) %>%
      pull(max_value)
    max_y <- max(max_y,results_long_mom_L2 %>%
                   filter(metric == met) %>%summarise(max_value = max(value)) %>%
                   pull(max_value) )
    # 绘制Bellec-Ridge和MoM的图表
    plot <- ggplot() +
      geom_line(
        data = results_long_bellec_L2 %>% filter(metric == met),
        aes(
          x = n,
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
        aes(x = n, y = value),
        color = "black",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_L2 %>% filter(metric == met),
        aes(x = n, y = value),
        color = "black",
        linetype = "dashed",
        size = 1
      )  +
      labs(
        title = substitute(paste(alpha^T, Sigma, alpha, sep = '')),
        x = "n",
        y = metrics[met]
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            # 设置标题居中
            legend.position = "right") +
      scale_x_continuous(
        breaks = ns  # 设置显示特定的n值
      ) +
      coord_cartesian(ylim = c(NA, max_y))  # 限制y轴的范围
    
    # 存储每个指标的图表
    plots[[3 * plot_index ]] <- plot
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
  p_over_n <- round(ps_alpha[1] /ns_alpha[1], 2)
  
  # 生成 PDF 文件名
  pdf_filename <- paste0(
    "plots_bellec_lambda/L2_rootn_bias_glm_mom_sparse_",
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
  pdf(pdf_filename, width = 8.85*3/2, height = 6.6375)
  grid.draw(combined_plot)
  # 关闭PDF设备
  dev.off()
}