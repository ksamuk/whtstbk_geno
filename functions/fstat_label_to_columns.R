fstat_label_to_columns <- function(x){
  if(!is.null(names(x))){
    x <- names(x)
  }
  chr <- gsub("\\..*", "", x)
  pos <- gsub("[^0-9]", "", x)
  data.frame(chr, pos)
}