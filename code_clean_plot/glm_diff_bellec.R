rm(list = ls())

library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(grid)

# c("0_0_0", "0_0_1", "1_0_0", "1_0_1", "1_1_1")
# 定义三元组
triplets <- c("0_0_0")
file_path_be <- file.path("data_bellec_1.2_clean_100")
file_path_mom <- file.path("data_glm_diff")


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
  alpha_L2_bellec <-  matrix(
    NA,
    nrow = N.simu * N.lambda,
    ncol = nrow(bellec_data_total[[1]]$alpha_de_p2_TT_Ge[[1]])
  )
  
  alpha_de_N_p2_TT_Ge_nonzero <- rep(NA, 2 * N.simu * N.lambda)
  alpha_de_N_p2_TT_Ge_zero <- rep(NA, 2 * N.simu * N.lambda)
  
  alpha_de_p2_TT_Ge_nonzero <- rep(NA, 2 * N.simu * N.lambda)
  alpha_de_p2_TT_Ge_zero <- rep(NA, 2 * N.simu * N.lambda)
  
  for (i in c(1:N.simu)) {
    simulations <- bellec_data_total[[i]]
    alpha_term <- simulations$alpha[indices_selected]
    
    for (j in c(1:N.lambda)) {
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
    alpha_de_p2_TT_Ge_zero
  )
  
  # 定义每个 n 对应的样本数
  sample_size <- 80
  # 对每个 n 进行抽样
  table_bellec <- table_bellec %>%
    group_by(ns_alpha, lambda_alpha) %>%
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
  # mom -----
  

  file_names <- list.files(file_path_mom, pattern = "\\.RDa$", full.names = TRUE)
  
  ns <- ps <- rep(0, length(file_names))
  alpha_em_New_7_zero <- alpha_em_New_7_nonzero <- alpha_L2_New_7 <-  list()
  alpha_em_New_6_zero <- alpha_em_New_6_nonzero <- alpha_L2_New_6 <-  list()
  alpha_em_New_4_zero <- alpha_em_New_4_nonzero <- alpha_L2_New_4 <-  list()
  alpha_em_Old_4_zero <- alpha_em_Old_4_nonzero <-  list()

  for (j in seq_along(file_names)) {
    load(file_names[j])
    alpha_em_New_7_zero[[j]] <-   alpha_em_New_7[, 100] - alpha[100]
    alpha_em_New_7_nonzero[[j]] <-   alpha_em_New_7[, 1] - alpha[1]
    alpha_L2_New_7[[j]] <- sols_em_7[,2]
    
    alpha_em_New_6_zero[[j]] <-   alpha_em_New_6[, 100] - alpha[100]
    alpha_em_New_6_nonzero[[j]] <-   alpha_em_New_6[, 1] - alpha[1]
    alpha_L2_New_6[[j]] <- sols_em_6[,2]
    
    alpha_em_New_4_zero[[j]] <-   alpha_em_New_4[, 100] - alpha[100]
    alpha_em_New_4_nonzero[[j]] <-   alpha_em_New_4[, 1] - alpha[1]
    alpha_L2_New_4[[j]] <- sols_em_4[,2]
    
    ns[j] <- n
    ps[j] <- p
  }
  results_mom <- data.frame(
    n = integer(),
    lambda = numeric(),
    bias = numeric(),
    variance = numeric(),
    mse = numeric(),
    method = character()
  )
  
  # 加载数据
  for (j in seq_along(file_names)) {
    load(file_names[j])
    alpha_em_New_7_zero[[j]] <- alpha_em_New_7[, 100] - alpha[100]
    alpha_em_New_7_nonzero[[j]] <- alpha_em_New_7[, 1] - alpha[1]
    alpha_L2_New_7[[j]] <- sols_em_7[, 2] - par_truth_7[2]
    
    alpha_em_New_6_zero[[j]] <- alpha_em_New_6[, 100] - alpha[100]
    alpha_em_New_6_nonzero[[j]] <- alpha_em_New_6[, 1] - alpha[1]
    alpha_L2_New_6[[j]] <- sols_em_6[, 2] - par_truth_7[2]
    
    alpha_em_New_4_zero[[j]] <- alpha_em_New_4[, 100] - alpha[100]
    alpha_em_New_4_nonzero[[j]] <- alpha_em_New_4[, 1] - alpha[1]
    alpha_L2_New_4[[j]] <- sols_em_4[, 2] - par_truth_7[2]
    
    ns[j] <- n
    ps[j] <- p
  }
  
  # 初始化结果数据框
  results_mom <- data.frame(
    n = integer(),
    lambda = numeric(),
    bias = numeric(),
    variance = numeric(),
    mse = numeric(),
    method = character()
  )
  
  # 通用的计算统计量的函数
  calculate_stats <- function(simulations, n, method) {
    bias <- mean( simulations, na.rm = TRUE)
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations) ^ 2, na.rm = TRUE)
    
    data.frame(
      n = n,
      lambda = NA,
      bias = bias,
      variance = variance,
      mse = mse,
      method = method
    )
  }
  
  # 定义需要处理的alpha列表和方法名
  alpha_lists <- list(
    list(name = "alpha_em_New_7_zero", method = "em-New-7-zero"),
    list(name = "alpha_em_New_7_nonzero", method = "em-New-7-nonzero"),
    list(name = "alpha_L2_New_7", method = "L2-New-7"),
    list(name = "alpha_em_New_6_zero", method = "em-New-6-zero"),
    list(name = "alpha_em_New_6_nonzero", method = "em-New-6-nonzero"),
    list(name = "alpha_L2_New_6", method = "L2-New-6"),
    list(name = "alpha_em_New_4_zero", method = "em-New-4-zero"),
    list(name = "alpha_em_New_4_nonzero", method = "em-New-4-nonzero"),
    list(name = "alpha_L2_New_4", method = "L2-New-4")
  )
  
  # 计算每个alpha列表的统计量并合并结果
  for (alpha_list in alpha_lists) {
    list_name <- alpha_list$name
    method_name <- alpha_list$method
    for (i in seq_along(get(list_name))) {
      simulations <- get(list_name)[[i]]
      stats <- calculate_stats(simulations, ns[i], method_name)
      results_mom <- rbind(results_mom, stats)
    }
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
  
  
  # 对 results_long_bellec 进行处理
  results_long_bellec <- results_long_bellec %>%
    mutate(value = if_else(metric == "bias", value * sqrt(n), value))
  
  # 对 results_long_mom 进行处理
  results_long_mom <- results_long_mom %>%
    mutate(value = if_else(metric == "bias", value * sqrt(n), value))
  
  
  results_long_mom_4_zero <- results_long_mom[results_long_mom$method == "em-New-4-zero", ]
  results_long_mom_4_nonzero <- results_long_mom[results_long_mom$method == "em-New-4-nonzero", ]
  results_long_mom_4_L2 <- results_long_mom[results_long_mom$method == "L2-New-4", ]
  
  
  results_long_mom_6_zero <- results_long_mom[results_long_mom$method == "em-New-6-zero", ]
  results_long_mom_6_nonzero <- results_long_mom[results_long_mom$method == "em-New-6-nonzero", ]
  results_long_mom_6_L2 <- results_long_mom[results_long_mom$method == "L2-New-6", ]
  
  
  results_long_mom_7_zero <- results_long_mom[results_long_mom$method == "em-New-7-zero", ]
  results_long_mom_7_nonzero <- results_long_mom[results_long_mom$method == "em-New-7-nonzero", ]
  results_long_mom_7_L2 <- results_long_mom[results_long_mom$method == "L2-New-7", ]
  
  
  results_long_bellec_zero <- results_long_bellec[results_long_bellec$method == "alpha_de_p2_TT_Ge_zero", ]
  results_long_bellec_nonzero <- results_long_bellec[results_long_bellec$method == "alpha_de_p2_TT_Ge_nonzero", ]
  # results_long_bellec_L2 <- results_long_bellec[results_long_bellec$method == "L2-MoM", ]
  
  # plot -----
  library(gridExtra)
  library(ggplot2)
  library(dplyr)
  
  # 定义度量指标
  metrics <- c("bias" = "Root n Bias",
               "variance" = "Variance",
               "mse" = "Mean Square Error")
  
  # 准备绘图
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
    max_y <- max(max_y, results_long_mom_4_nonzero %>%
                   filter(metric == met) %>% summarise(max_value = max(value)) %>%
                   pull(max_value))
    
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
      labs(
        title = expression(alpha[1] ~ "(nonzero)"),
        x = "n",
        y = metrics[met]
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") +
      scale_x_continuous(
        breaks = ns  # 设置显示特定的n值
      ) +
      coord_cartesian(ylim = c(NA, max_y))  # 限制y轴的范围
    
    # 后缀4的数据
    plot <- plot +
      geom_point(
        data = results_long_mom_4_nonzero %>% filter(metric == met),
        aes(x = n, y = value),
        color = "black",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_4_nonzero %>% filter(metric == met),
        aes(x = n, y = value, linetype = "Suffix 4"),
        color = "black",
        size = 1.5
      )
    
    # 后缀6的数据
    plot <- plot +
      geom_point(
        data = results_long_mom_6_nonzero %>% filter(metric == met),
        aes(x = n, y = value),
        color = "green",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_6_nonzero %>% filter(metric == met),
        aes(x = n, y = value, linetype = "Suffix 6"),
        color = "green",
        size = 1.5
      )
    
    # 后缀7的数据
    plot <- plot +
      geom_point(
        data = results_long_mom_7_nonzero %>% filter(metric == met),
        aes(x = n, y = value),
        color = "purple",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_7_nonzero %>% filter(metric == met),
        aes(x = n, y = value, linetype = "Suffix 7"),
        color = "purple",
        size = 1.5
      )
    
    # 添加图例
    plot <- plot + scale_linetype_manual(
      values = c("Suffix 4" = "dashed", "Suffix 6" = "dotted", "Suffix 7" = "solid"),
      name = "Suffix"
    )
    
    # 存储每个指标的图表
    plots[[2 * plot_index - 1]] <- plot
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
    max_y <- max(max_y, results_long_mom_4_zero %>%
                   filter(metric == met) %>% summarise(max_value = max(value)) %>%
                   pull(max_value))
    
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
      labs(
        title = expression(alpha[1] ~ "(zero)"),
        x = "n",
        y = metrics[met]
      ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = "right") +
      scale_x_continuous(
        breaks = ns  # 设置显示特定的n值
      ) +
      coord_cartesian(ylim = c(NA, max_y))  # 限制y轴的范围
    
    # 后缀4的数据
    plot <- plot +
      geom_point(
        data = results_long_mom_4_zero %>% filter(metric == met),
        aes(x = n, y = value),
        color = "black",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_4_zero %>% filter(metric == met),
        aes(x = n, y = value, linetype = "Suffix 4"),
        color = "black",
        size = 1.5
      )
    
    # 后缀6的数据
    plot <- plot +
      geom_point(
        data = results_long_mom_6_zero %>% filter(metric == met),
        aes(x = n, y = value),
        color = "green",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_6_zero %>% filter(metric == met),
        aes(x = n, y = value, linetype = "Suffix 6"),
        color = "green",
        size = 1.5
      )
    
    # 后缀7的数据
    plot <- plot +
      geom_point(
        data = results_long_mom_7_zero %>% filter(metric == met),
        aes(x = n, y = value),
        color = "purple",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_7_zero %>% filter(metric == met),
        aes(x = n, y = value, linetype = "Suffix 7"),
        color = "purple",
        size = 1.5
      )
    
    # 添加图例
    plot <- plot + scale_linetype_manual(
      values = c("Suffix 4" = "dashed", "Suffix 6" = "dotted", "Suffix 7" = "solid"),
      name = "Suffix"
    )
    
    # 存储每个指标的图表
    plots[[2 * plot_index]] <- plot
    plot_index <- plot_index + 1
  }
  
  # 合并图表为单列布局
  combined_plot <- do.call(grid.arrange, c(plots, ncol = 2))
  print(combined_plot)
  
  # 提取三元组的值
  triplet_parts <- strsplit(triplet, "_")[[1]]
  Is_sparse <- as.numeric(triplet_parts[1])
  Is_sparse_only_one <- as.numeric(triplet_parts[2])
  Is_Rad <- as.numeric(triplet_parts[3])
  
  # 计算 p_over_n 并保留两位小数
  p_over_n <- round(ps_alpha[1] / ns_alpha[1], 2)
  
  # 生成 PDF 文件名
  pdf_filename <- paste0(
    "plots_bellec_lambda/diff-rootn_bias_glm_mom_sparse_",
    Is_sparse,
    "_one_",
    Is_sparse_only_one,
    "_Ra_",
    Is_Rad,
    "_p_over_n_",
    p_over_n,
    ".pdf"
  )
  
  plot_index <- 1
  for (met in names(metrics)) {
    
    # 后缀4的数据
    plot <- plot +
      geom_point(
        data = results_long_mom_4_L2 %>% filter(metric == met),
        aes(x = n, y = value),
        color = "black",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_4_L2 %>% filter(metric == met),
        aes(x = n, y = value, linetype = "Suffix 4"),
        color = "black",
        size = 1.5
      )
    
    # 后缀6的数据
    plot <- plot +
      geom_point(
        data = results_long_mom_6_L2 %>% filter(metric == met),
        aes(x = n, y = value),
        color = "green",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_6_L2 %>% filter(metric == met),
        aes(x = n, y = value, linetype = "Suffix 6"),
        color = "green",
        size = 1.5
      )
    
    # 后缀7的数据
    plot <- plot +
      geom_point(
        data = results_long_mom_7_L2 %>% filter(metric == met),
        aes(x = n, y = value),
        color = "purple",
        shape = 21,
        size = 3,
        fill = "white"
      ) +
      geom_line(
        data = results_long_mom_7_L2 %>% filter(metric == met),
        aes(x = n, y = value, linetype = "Suffix 7"),
        color = "purple",
        size = 1.5
      )
    
    # 添加图例
    plot <- plot + scale_linetype_manual(
      values = c("Suffix 4" = "dashed", "Suffix 6" = "dotted", "Suffix 7" = "solid"),
      name = "Suffix"
    )
    
    # 存储每个指标的图表
    plots[[3 * plot_index]] <- plot
    plot_index <- plot_index + 1
  }
  
  # 打开PDF设备
  pdf(pdf_filename, width = 8.85, height = 6.6375)
  grid.draw(combined_plot)
  # 关闭PDF设备
  dev.off()
  
  
}