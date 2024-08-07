library(dplyr)

# 设置文件路径和输出路径
file_path <- "data/data_glm_mom"
output_path <- "data/data_glm_mom_after_cluster"
# 确保输出路径存在
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

# 获取所有文件名
file_names <- list.files(file_path, pattern = "\\.Rda$", full.names = TRUE)

# 提取文件名前缀
prefixes <- unique(gsub("(.*)_\\d{8}_\\d{6}\\.Rda", "\\1", basename(file_names)))

# 遍历每个前缀
for (prefix in prefixes) {
  # 提取 n, p, sparse, omega 的值
  n_value <- as.numeric(gsub(".*n_(\\d+).*", "\\1", prefix))
  p_value <- as.numeric(gsub(".*p_(\\d+).*", "\\1", prefix))
  sparse_value <- as.numeric(gsub(".*sparse_(\\d+).*", "\\1", prefix))
  omega_value <- as.numeric(gsub(".*omegma(\\d+).*", "\\1", prefix))
  Is_Rad <- as.logical(as.numeric(gsub(".*Ra_(\\d).*", "\\1", prefix)))

  # 计算 n_iter * length(current_files)
  n_iter <- as.numeric(gsub(".*iter_(\\d+).*", "\\1", prefix))


  # 判断 Is_sparse 和 Is_sparse_only_one
  Is_sparse_only_one <- sparse_value == 1
  Is_sparse <- sparse_value < p_value

  # 获取当前前缀的所有文件
  current_files <- file_names[grepl(paste0("^", prefix), basename(file_names))]
  total_iter <- n_iter * length(current_files)


  alpha_MoM_list_nonzero <- list()
  alpha_MoM_list_zero <- list()
  alpha_L2_MoM_list <- list()
  par_truth_list <- sols_em_list <- list()
  alpha_est_N_total <- list()
  moments_em_list <- moments_truth_list <- list()
  for (j in seq_along(current_files)) {
    load(current_files[j])

    alpha_est_N_total[[j]] <- alpha_est_N
    alpha_MoM_list_zero[[j]] <- alpha_est_N[, 11] - alpha[indices_selected[11]]
    alpha_MoM_list_nonzero[[j]] <- alpha_est_N[, 1] - alpha[1]
    alpha_L2_MoM_list[[j]] <- alpha_L2_est_N - par_truth[2]
    sols_em_list[[j]] <- sols_em
    moments_em_list[[j]] <- moments_em
    moments_truth_list[[j]] <- matrix(rep(moments_truth, each = nrow(moments_em)), nrow = nrow(moments_em))
    par_truth_list[[j]] <- matrix(rep(par_truth, each = nrow(sols_em)), nrow = nrow(sols_em))
  }
  sols_em_matrix <- do.call(rbind, sols_em_list)
  par_truth_matrix <- do.call(rbind, par_truth_list)
  moments_em_matrix <- do.call(rbind, moments_em_list)
  moments_truth_matrix <- do.call(rbind, moments_truth_list)
  alpha_est_N_total <- do.call(rbind, alpha_est_N_total)

  alpha_MoM_list_zero <- as.vector(unlist(alpha_MoM_list_zero))
  alpha_MoM_list_nonzero <- as.vector(unlist(alpha_MoM_list_nonzero))
  alpha_L2_MoM_list <- as.vector(unlist(alpha_L2_MoM_list))

  results_mom <- data.frame(
    n = integer(),
    lambda = numeric(),
    bias = numeric(),
    variance = numeric(),
    mse = numeric(),
    method = character()
  )

  # 计算alpha_MoM_list的统计量
    simulations <- alpha_MoM_list_zero
    bias <- mean(simulations, na.rm = TRUE)
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations) ^ 2, na.rm = TRUE)

    results_mom <- rbind(
      results_mom,
      data.frame(
        n = n_value,
        lambda = NA,
        bias = bias,
        variance = variance,
        mse = mse,
        method = "MoM-zero"
      )
    )
    simulations <- alpha_MoM_list_nonzero
    bias <- mean(simulations, na.rm = TRUE)
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations) ^ 2, na.rm = TRUE)

    results_mom <- rbind(
      results_mom,
      data.frame(
        n = n_value,
        lambda = NA,
        bias = bias,
        variance = variance,
        mse = mse,
        method = "MoM-nonzero"
      )
    )

    simulations <- alpha_L2_MoM_list
    bias <- mean(simulations, na.rm = TRUE)
    variance <- var(simulations, na.rm = TRUE)
    mse <- mean((simulations) ^ 2, na.rm = TRUE)

    results_mom <- rbind(
      results_mom,
      data.frame(
        n = n_value,
        lambda = NA,
        bias = bias,
        variance = variance,
        mse = mse,
        method = "MoM-L2"
      )
    )


  # 保存合并后的对象
  save(
    n_value,
    p_value,
    alpha,
    sols_em_matrix,
    par_truth_matrix,
    Is_sparse,
    Is_sparse_only_one,
    Is_Rad,
    moments_truth_matrix,
    moments_em_matrix,
    results_mom,
    alpha_est_N_total,
    R.version,
    glmnet.version,
    system.version,

    file = file.path(output_path, paste0(
      "MoM_n_",
      n_value,
      "_p_",
      p_value,
      "_sparse_",
      as.numeric(Is_sparse),
      "_one_",
      as.numeric(Is_sparse_only_one),
      "_Ra_",
      as.numeric(Is_Rad),
      "_sim_",
      total_iter,
      ".Rda"
    ))
  )
}

cat("合并和保存完成。\n")

library(dplyr)
library(purrr)


# 获取所有文件名
file_names <- list.files(output_path, pattern = "\\.Rda$", full.names = TRUE)

# 提取三元组
triplets <- unique(gsub(".*sparse_(\\d)_one_(\\d)_Ra_(\\d).*", "\\1_\\2_\\3", basename(file_names)))

# 遍历每个三元组
for (triplet in triplets) {
  # 创建对应的三元组文件夹
  triplet_folder <- file.path(output_path, triplet)
  if (!dir.exists(triplet_folder)) {
    dir.create(triplet_folder)
  }
  
  # 获取当前三元组的所有文件
  triplet_parts <- strsplit(triplet, "_")[[1]]
  sparse_value <- triplet_parts[1]
  one_value <- triplet_parts[2]
  Ra_value <- triplet_parts[3]
  
  current_files <- file_names[grepl(paste0(".*sparse_", sparse_value, "_one_", one_value, "_Ra_", Ra_value, ".*"), basename(file_names))]
  
  # 移动文件到对应的三元组文件夹
  file.rename(from = current_files, to = file.path(triplet_folder, basename(current_files)))
}

cat("文件分组和移动完成。\n")
