
plot_resid_fixed <- function(stem, pop_order, min = -0.009, max = 0.009, cex = 1, usemax = T, wcols = "r"){
  c = read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  m = read.table(gzfile(paste(stem, ".modelcov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  names(c) = rownames(c)
  names(m) = rownames(m)
  #o = read.table(pop_order, as.is = T, comment.char = "", quote = "")
  o = data.frame(pop_order)
  se = read.table(gzfile(paste(stem, ".covse.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  mse = apply(se, 1, mean)
  mse = mean(mse)
  print(mse)	
  c = c[order(names(c)), order(names(c))]
  m = m[order(names(m)), order(names(m))]
  tmp = c -m 
  #tmp = m - c
  #tmp = (m-c)/m
  #print(tmp)
  toplot = data.frame(matrix(nrow = nrow(tmp), ncol = ncol(tmp)))
  for(i in 1:nrow(o)){
    for( j in 1:nrow(o)){
      #print(paste(o[i,1], o[j,1]))
      if (o[i,1] %in% names(tmp) ==F){
        print(paste("not found", o[i,1]))
      }
      if (o[j,1] %in% names(tmp) ==F){
        print(paste("not found", o[j,1]))
      }
      toplot[i, j] = tmp[which(names(tmp)==o[i,1]), which(names(tmp)==o[j,1])]
    }
  }
  #print(toplot)
  if (usemax){
    m1 = max(abs(toplot), na.rm = T)
    max = m1*1.02
    min = -(m1*1.02)	
  }
  print("here")
  names(toplot) = o[,1]
  toreturn = plot_resid_internal(toplot, max = max, min = min, wcols = wcols, mse = mse, o = o, cex = cex)
  return(toreturn)
}