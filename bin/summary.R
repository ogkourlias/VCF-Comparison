read_files <- function(file_name) {
  df <- read.csv(file_name, sep = "\t")
  df$file_name <- file_name
  return(df)
}
file_names <- list.files(pattern = "../*.tsv")
df_list <- lapply(file_names, read_files)
chr21 <- df_list[0]
chr21$file_name