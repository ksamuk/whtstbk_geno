fstat_label_to_columns <- function(x){
  x <- names(x)
  chr <- gsub("\\..*", "", x)
  pos <- gsub("[^0-9]", "", x)
  data.frame(chr, pos)
}