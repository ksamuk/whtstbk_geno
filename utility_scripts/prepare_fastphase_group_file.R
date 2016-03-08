# make a fastphase popfile

fp_pop <- read.table("metadata/pop_file_fastphase.txt", header = FALSE)

fp_pop_file <- as.numeric(fp_pop[,2])

write(fp_pop_file, "metadata/fastphase_labels.txt", ncolumns = length(fp_pop_file))
