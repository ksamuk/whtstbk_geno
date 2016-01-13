# rename structure output

files <- list.files("data/structure/results", full.names = TRUE)

for (i in 1:length(files)){
  file_name_new <- gsub("K[0-9]_","", files[i])
  file.rename(files[i], file_name_new)
}
