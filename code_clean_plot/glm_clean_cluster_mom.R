library(dplyr)
library(purrr)

# 设置文件路径和输出路径
file_path <- "data/data_glm_mom"
output_path <- "data/data_glm_mom_after_cluster"

# 获取所有文件名
file_names <- list.files(file_path, pattern = "\\.Rda$", full.names = TRUE)

# 提取三元组
triplets <- unique(gsub(".*sparse_(\\d)_one_(\\d)_Rad_(\\d).*", "\\1_\\2_\\3", basename(file_names)))

# 遍历每个三元组
for (triplet in triplets) {
  # 创建对应的三元组文件夹
  triplet_folder <- file.path(output_path, triplet)
  if (!dir.exists(triplet_folder)) {
    dir.create(triplet_folder)
  }
  
  # 获取当前三元组的所有文件
  triplet_parts <- strsplit(triplet, "_")[[1]]
  Is_sparse <- triplet_parts[1]
  Is_sparse_only_1 <- triplet_parts[2]
  Is_Rad <- triplet_parts[3]
  
  current_files <- file_names[grepl(paste0(".*sparse_", Is_sparse, "_one_", Is_sparse_only_1, "_Rad_", Is_Rad, ".*"), basename(file_names))]
  
  # 移动文件到对应的三元组文件夹
  file.rename(from = current_files, to = file.path(triplet_folder, basename(current_files)))
}

cat("文件分组和移动完成。\n")