library(dplyr)
library(purrr)

# 设置文件路径和输出路径
file_path <- "data/data_mar_mom"

# 获取所有文件名
file_names <- list.files(file_path , pattern = "\\.Rda$", full.names = TRUE)

# 提取三元组
triplets <- unique(gsub(".*sparse_(\\d)_Ra_(\\d).*", "\\1_\\2", basename(file_names)))

# 遍历每个三元组
for (triplet in triplets) {
  # 创建对应的三元组文件夹
  triplet_folder <- file.path(file_path, triplet)
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