# the speciation continuum
# is real is it false, who cares! 


library("ggplot2")

fst_df <- read.table("metadata/fst_df.txt", h = T)

fst_df <- read.csv("metadata/fst_df_reduced.csv", h = T)

head(fst_df)


fst_df %>%
  ggplot(aes(x = fst)) +
  geom_histogram(bins = 5, binwidth = 0.02)
