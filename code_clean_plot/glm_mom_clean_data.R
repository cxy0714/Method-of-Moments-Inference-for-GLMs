library(dplyr)

# 设置文件路径和输出路径
file_path <- "data_mom_one_1.2_500"
output_path <- "data_norm_1.2_clean_500"

# 确保输出路径存在
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

# 获取所有文件名
file_names <- list.files(file_path, pattern = "\\.Rda$", full.names = TRUE)

# 提取文件名前缀
prefixes <- unique(gsub("(.*)_\\d{8}_\\d{6}\\.Rda", "\\1", basename(file_names)))
# 遍历每个前缀
# 遍历每个前缀
for (prefix in prefixes) {
  # 提取 n, p, sparse, omega 的值
  n_value <- as.numeric(gsub(".*n_(\\d+).*", "\\1", prefix))
  p_value <- as.numeric(gsub(".*p_(\\d+).*", "\\1", prefix))
  sparse_value <- as.numeric(gsub(".*sparse_(\\d+).*", "\\1", prefix))
  omega_value <- as.numeric(gsub(".*omegma(\\d+).*", "\\1", prefix))
  Is_Rad <- as.logical(as.numeric(gsub(".*Ra_(\\d).*", "\\1", prefix)))
  
  # 判断 Is_sparse 和 Is_sparse_only_one
  Is_sparse_only_one <- sparse_value == 1
  Is_sparse <- sparse_value < p_value
  
  # 获取当前前缀的所有文件
  current_files <- file_names[grepl(paste0("^", prefix), basename(file_names))]
  
  # 初始化合并对象
  mu_alpha_est_total <- numeric()
  alpha_L2_est_N_total <- numeric()
  m_em_total <- matrix(nrow = 0, ncol = 4)
  alpha_est_N_total <- matrix(nrow = 0, ncol = 20)
  seeds_total <- unique_seed_total <-  numeric()

  # 读取并合并文件
  for (file in current_files) {
    load(file)
    seeds_total <-  c(seeds_total,seeds)
    unique_seed_total <- c(unique_seed_total,unique_seed)
    if (exists("mu_alpha_est")) {
      mu_alpha_est_total <- c(mu_alpha_est_total, mu_alpha_est)
      
    }
    if (exists("alpha_L2_est_N")) {
      alpha_L2_est_N_total <- c(alpha_L2_est_N_total, alpha_L2_est_N)
    }
    if (exists("m_em")) {
      m_em_total <- rbind(m_em_total, m_em)
    }
    if (exists("alpha_est_N")) {
      alpha_est_N_total <- rbind(alpha_est_N_total, alpha_est_N)
    }
  }
  
  # 计算 n_iter * length(current_files)
  n_iter <- as.numeric(gsub(".*iter_(\\d+).*", "\\1", prefix))
  total_iter <- n_iter * length(current_files)
  
  # 保存合并后的对象
  save(
    n_value,
    p_value,
    alpha,
    
    Is_sparse,
    Is_sparse_only_one,
    Is_Rad,
    
    m_em_total,
    m_truth,
    
    par_truth,
    mu_alpha_est_total,
    alpha_L2_est_N_total,
    
    alpha_est_N_total,
    
    R.version,
    glmnet.version,
    system.version,

    unique_seed_total,
    seeds_total,
    indices_selected,
    
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