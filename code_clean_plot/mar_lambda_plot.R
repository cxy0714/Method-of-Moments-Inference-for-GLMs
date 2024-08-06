rm(list = ls())


library(tidyverse)
library(dplyr)
library(reshape)
library(viridis)
library(stringr)
library(grid)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(latex2exp)

# 定义三元组
# c("0_0", "0_1", "1_0")
triplets <- c("0_0")
file.path_mom <- "data/data_mar_mom"

# 遍历每个三元组
for (triplet in triplets) {
  
  ## mom data process -----
  file_path <- file.path(file.path_mom, triplet)
  file_names <- list.files(file_path, pattern = "\\.Rda", full.names = TRUE)
  n_n <- p_n <-  rep(NA, 10)
  
  psi_est_list <- list()
  for (j in seq_along(file_names)) {
    # N.replicates <- as.numeric(gsub(".*sim_(\\d+).*", "\\1", file_names[j]))
    load(file_names[j])
    psi_est_list[[j]] <- data.frame(n = rep(n, N.replicates), MoM = sols_em[, 4])
    n_n[j] <- n
    p_n[j] <- p
    
  }
  
  psi_est <- do.call(rbind, psi_est_list)
  # 提取 plot_data 中包含的 estimator 列表
  psi_stats <- psi_est %>%
    group_by(n) %>%
    summarize(
      mean =   mean(sqrt(n) * MoM , na.rm = TRUE),
      var = var(MoM, na.rm = TRUE),
      mse = mean(MoM ^ 2, na.rm = TRUE),
    ) %>%
    ungroup() %>%
    as.data.frame()
  
  # 为绘图准备数据
  plot_data_psi <- bind_rows(
    mutate(psi_stats, metric = "mean", value = mean),
    mutate(psi_stats, metric = "var", value = var),
    mutate(psi_stats, metric = "mse", value = mse)
  )
  
  # 将 psi_est 数据标记为 "PSI Estimator"，以便在图中区分
  plot_data_psi <- mutate(plot_data_psi, estimator = "MoM-Stein", lambda.y = NA)
  plot_data_psi <- plot_data_psi %>% filter(n %in% unique(plot_data_psi$n))
  
  plot_data_psi[plot_data_psi$metric == "mean", ]$metric <- "Root n Bias"
  plot_data_psi[plot_data_psi$metric == "var", ]$metric <- "Variance"
  plot_data_psi[plot_data_psi$metric == "mse", ]$metric <- "Mean Square Error"
  
  ## celentano data process lambda------
  
  # Load  data
  all_data.df <- read.csv(paste0("data_celetano/lambda_sparse_",as.numeric(Is_sparse),"_Rad_",as.numeric(Is_Rad),".csv"))
  # lambda_sparse_ 1 _Rad_ 1 _ 20240720_225529 .csv
  # lambda_sparse_ 0 _Rad_ 1 _ 20240718_000608 .csv
  # lambda_sparse_ 0 _ 20240716_003208 .csv
  # lambda_sparse_ 1 _ 20240713_200744 .csv
  
  head(all_data.df)
  nrow(all_data.df)
  print(unique(all_data.df$n))
  print(unique(all_data.df$lambda.y))
  
  
  
  category_order <- c(
    "Ridge",
    "Ridge (IPW)",
    "Debiased ridge (naive)",
    "Debiased ridge (IPW)",
    "Oracle ASCW",
    "Empirical SCA",
    "MoM-Stein"
  )
  
  # Define the estimator renaming function
  estimator.rename <- function(estimator) {
    if (estimator == "mu.y.hat") {
      return("Ridge")
    } else if (estimator == "mu.y.hat.ipw") {
      return("Ridge (IPW)")
    } else if (estimator == "mu.y.hat.d.naive") {
      return("Debiased ridge (naive)")
    } else if (estimator == "mu.y.hat.d.ipw") {
      return("Debiased ridge (IPW)")
    } else if (estimator == "mu.y.hat.d.cfd") {
      return("Oracle ASCW")
    } else if (estimator == "mu.y.hat.d.emp.moment") {
      return("Empirical SCA")
    } else if (estimator == "mu.y.hat.MoM") {
      return("MoM-Stein")
    } else {
      return(estimator)
    }
  }
  estimator.rename = Vectorize(estimator.rename)
  
  # Process data
  plotting.df <- select(
    all_data.df,
    c(
      n,
      lambda.y,
      mu.y.hat,
      mu.y.hat.ipw,
      mu.y.hat.d.naive,
      mu.y.hat.d.cfd,
      mu.y.hat.d.ipw,
      mu.y.hat.d.emp.moment
    )
  ) %>%
    melt(id = c("n", "lambda.y")) %>%
    group_by(n, lambda.y, variable) %>%
    summarize(
      mean = mean(sqrt(n) * value, na.rm = TRUE),
      var = var(value, na.rm = TRUE),
      mse = mean(value ^ 2, na.rm = TRUE),
      mean_se = sd(value, na.rm = TRUE) / sqrt(n()),
      var_se = sqrt(var((
        value - mean(value)
      ) ^ 2) / (n() - 2)),
      mse_se = NA
    ) %>%
    ungroup() %>%
    as.data.frame() %>%
    dplyr::rename(estimator = variable) %>%
    mutate(estimator = estimator.rename(estimator)) %>%
    melt(id = c("n", "lambda.y", "estimator", "mean_se", "var_se")) %>%
    dplyr::rename(metric = variable, value = value) %>%
    mutate(
      se = if_else(metric == "mean", mean_se, var_se),
      estimator = factor(estimator, levels = category_order)
    ) %>%
    select(c(n, lambda.y, estimator, metric, value, se))
  
  head(plotting.df)
  
  # 从plotting.df中提取所有这些最优lambda对应的数据
  plot_data <- plotting.df %>% filter(estimator %in% c("Oracle ASCW", "Empirical SCA"))
  ## plot ----
  
  
  plot_data <- plot_data %>%
    mutate(metric = recode(
      metric,
      "mean" = "Root n Bias",
      "var" = "Variance",
      "mse" = "Mean Square Error"
    ))
  # 假设plot_data和plot_data_psi已正确加载
  estimators <- unique(plot_data$estimator)
  
  plots <- list()
  plot_index <- 1
  # 为每个 estimator 创建图表
  for (met in c("Root n Bias", "Variance", "Mean Square Error")) {
    for (est in estimators)  {
      y_label <- met
      
      
      plot <- ggplot() +
        geom_line(
          data = plot_data %>% filter(estimator == est, metric == met),
          aes(
            x = log(n),
            y = value,
            color = lambda.y,
            group = lambda.y
          ),
          alpha = 0.5
        ) +
        scale_color_gradient(name = expression(lambda),
                             low = "blue",
                             high = "red") +
        # MoM 数据
        geom_point(
          data = plot_data_psi %>% filter(metric == met),
          aes(x = log(n), y = value),
          color = "black",
          shape = 21,
          size = 3,
          fill = "white"
        ) +
        geom_line(
          data = plot_data_psi %>% filter(metric == met),
          aes(x = log(n), y = value),
          color = "black",
          linetype = "dashed",
          size = 1
        ) +
        labs(
          title = paste(est),
          x = "n",
          y = y_label
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5)  # 设置标题居中
        ) +
        scale_x_continuous(
          labels = function(x) round(exp(x)),
          breaks = log(c(100, 1000, 10000))  # 设置显示特定的n值
        ) +
        guides(
          color = guide_colorbar(title = expression(lambda)),
          shape = guide_legend(title = "MoM-Stein"),
          linetype = guide_legend(title = "MoM-Stein")
        )
      plots[[plot_index]] <- plot
      plot_index <- plot_index + 1
    }
  }
  
  p_over_n <- round(p_value / n, 2)
  pdf_filename_lambda <- paste0(
    "plots_mar_lambda/6-glm_mar_lambda_sparse_",
    as.numeric(Is_sparse),
    "_Ra_",
    as.numeric(Is_Rad),
    "_p_over_n_",
    p_over_n,
    ".pdf"
  )
  
  
  
  combined_plot <- plot_grid(plotlist = plots, ncol = 2)
  
  # 保存为PDF
  pdf(pdf_filename_lambda, width = 8.85, height = 6.6375)
  grid.draw(combined_plot)
  dev.off()
  
}