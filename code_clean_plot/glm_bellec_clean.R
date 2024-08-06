library(dplyr)

# 设置文件路径和输出路径
file_path <- "data_bellec"
output_path <- "data_bellec_1.2_clean_100"

# 确保输出路径存在
if (!dir.exists(output_path)) {
  dir.create(output_path)
}

# 获取所有文件名
file_names <- list.files(file_path, pattern = "\\.Rda$", full.names = TRUE)

# 提取所有可能的三元组组合
triplets <- unique(gsub(".*(sparse_([01])_one_([01])_Rad_([01])).*", "\\1", basename(file_names)))

# 遍历每个三元组
for (triplet in triplets) {
  # 获取当前三元组的所有文件
  current_files <- file_names[grepl(triplet, basename(file_names))]
  
  if (length(current_files) == 0) {
    next
  }
  
  # 初始化合并对象
  bellec_data_total <- list()
  unique_seed_total <- numeric()
  
  # 读取并合并文件
  for (file in current_files) {
    load(file)
    
    bellec_data_total <- c(bellec_data_total, experiment.data)
    unique_seed_total <- c(unique_seed_total, unique_seed)
  }
  
  # 计算 n_iter * length(current_files)
  total_iter <- 2 * length(bellec_data_total)
  
  # 保存合并后的对象
  save(
    unique_seed_total,
    bellec_data_total,
    lambda_value,
    file = file.path(output_path, paste0(
      "Bellec_",
      triplet,
      "_sim_",
      total_iter,
      ".Rda"
    ))
  )
}

cat("合并和保存完成。\n")
