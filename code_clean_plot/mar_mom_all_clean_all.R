library(dplyr)
library(purrr)

# 设置文件路径和输出路径
file_path <- "data_mar_6_all"
output_path <- "data_mar_6_all_class"
# 获取所有文件名
file_names <- list.files(file_path , pattern = "\\.RDa$", full.names = TRUE)

# 提取三元组
triplets <- unique(gsub(".*sparse_(\\d)_Ra_(\\d).*", "\\1_\\2", basename(file_names)))

# 提取文件名前缀
prefixes <- unique(gsub("(.*)_\\d{8}_\\d{6}\\.RDa", "\\1", basename(file_names)))
# 遍历每个前缀
for (prefix in prefixes) {
  # 提取 n, p, sparse, omega 的值
  n_value <- as.numeric(gsub(".*n_(\\d+).*", "\\1", prefix))
  p_value <- as.numeric(gsub(".*p_(\\d+).*", "\\1", prefix))
  sparse_value <- as.numeric(gsub(".*sparse_(\\d+).*", "\\1", prefix))
  Is_sparse <- as.logical(sparse_value)
  Is_Rad <- as.logical(as.numeric(gsub(".*Ra_(\\d).*", "\\1", prefix)))
  
  # 获取当前前缀的所有文件
  current_files <- file_names[grepl(paste0("^", prefix), basename(file_names))]
  
  # 初始化合并对象

  moments_em_total <- matrix(nrow = 0, ncol = 6)
  sols_em_total <- matrix(nrow = 0, ncol = 4)
  N.success_total <- numeric()
  
  # 读取并合并文件
  for (file in current_files) {
    load(file)
    N.success_total <- c(N.success_total,N.success)
      sols_em_total <- rbind(sols_em_total, sols_em)
      moments_em_total <- rbind(moments_em_total, moments_em)
  }
  
  n_iter <- as.numeric(gsub(".*iter_(\\d+).*", "\\1", prefix))
  total_iter <- n_iter * length(current_files)
  
  # 保存合并后的对象
  save(
    n_value,
    p_value,
    alpha,
    beta,
    
    moments_truth,
    par_truth,
    Is_sparse,
    Is_Rad,
    
    sols_em_total,
    moments_em_total,
    
    R.version,
    glmnet.version,
    system.version,
    
    file = file.path(output_path, paste0(
      "n_",
      n_value,
      "_p_",
      p_value,
      "_sparse_",
      as.numeric(Is_sparse),
      "_Ra_",
      as.numeric(Is_Rad),
      "_sim_",
      total_iter,
      ".Rda"
    ))
  )
}
# 

# 获取所有文件名
file_names <- list.files(output_path, pattern = "\\.Rda$", full.names = TRUE)

# 提取三元组
triplets <- unique(gsub(".*sparse_(\\d)_Ra_(\\d).*", "\\1_\\2", basename(file_names)))

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
  Ra_value <- triplet_parts[2]

  current_files <- file_names[grepl(paste0(".*sparse_", sparse_value,"_Ra_", Ra_value, ".*"), basename(file_names))]

  # 移动文件到对应的三元组文件夹
  file.rename(from = current_files, to = file.path(triplet_folder, basename(current_files)))
}

cat("文件分组和移动完成。\n")