library(viridis)
library(shape)

plot_tree = function(stem, o = NA, cex = 1, disp = 0.003, plus = 0.01, flip = vector(), arrow = 0.05, scale = T, ybar = 0.1, mbar = T, plotmig = T, plotnames = T, xmin = 0, lwd = 1, arrow_lwd = 1, font = 1, spoof_labels = NA, spoof_cols = spoof_cols, use_viridis = TRUE, use_alpha = TRUE, arrow_lty = 1, shadow = 0.5){
  d = paste(stem, ".vertices.gz", sep = "")
  e = paste(stem, ".edges.gz", sep = "")
  se = paste(stem, ".covse.gz", sep = "")
  d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
  e = read.table(gzfile(e), as.is  = T, comment.char = "", quote = "")
  if (!is.na(o)){
    o = read.table(o, as.is = T, comment.char = "", quote = "")
  }
  e[,3] = e[,3]*e[,4]
  e[,3] = e[,3]*e[,4]
  
  se = read.table(gzfile(se), as.is = T, comment.char = "", quote = "")
  m1 = apply(se, 1, mean)
  m = mean(m1)
  #m = 0
  for(i in 1:length(flip)){
    d = flip_node(d, flip[i])
  }
  d$x = "NA"
  d$y = "NA"
  d$ymin = "NA"
  d$ymax = "NA"
  d$x = as.numeric(d$x)
  d$y = as.numeric(d$y)
  d$ymin = as.numeric(d$ymin)
  d$ymax = as.numeric(d$ymax)
  
  d = set_y_coords(d)
  d = set_x_coords(d, e)
  #print(d)
  d = set_mig_coords(d, e)
  
  if(any(!is.na(spoof_labels))){
    
    d[,2][!is.na(d[,2])] <- spoof_labels
  }
  plot_tree_internal(d, e, o = o, cex = cex, xmin = xmin, disp = disp, plus = plus, arrow = arrow, ybar = ybar, mbar = mbar, mse = m, scale = scale, plotmig = plotmig, plotnames = plotnames, lwd = lwd, font = font, arrow_lwd = arrow_lwd, spoof_cols = spoof_cols, use_viridis = use_viridis, use_alpha= use_alpha, arrow_lty = arrow_lty, shadow = shadow)
  #return(list( d= d, e = e))
}

plot_tree_internal = function(d, e, o = NA, cex = 1, disp = 0.005, plus = 0.005, arrow = 0.05, ybar = 0.01, scale = T, mbar = T, mse = 0.01, plotmig = T, plotnames = T, xmin = 0, lwd = 1, arrow_lwd = 1, font = 1, spoof_cols = spoof_cols, use_viridis = TRUE, use_alpha= TRUE, arrow_lty = 1, shadow = 0.5){
  
  plot(d$x, d$y, axes = F, ylab = "", xlab = "Drift parameter", xlim = c(xmin, max(d$x)+plus), pch = "")
  axis(1)
  mw = max(e[e[,5]=="MIG",4])
  
  if (use_viridis){
    mcols = viridis(150)
  }else if (use_alpha){
    mcols = as.character(lapply(0:149, function(x) rgb(0, 0, 0, alpha = (x/149))))
  }else{
    mcols = rev( heat.colors(150) )
  }

  for(i in 1:nrow(e)){
    col = "black"
    if (e[i,5] == "MIG"){
      w = floor(e[i,4]*200)+50
      if (mw > 0.5){
        w = floor(e[i,4]*100)+50
      }
      col = mcols[w]
      if (is.na(col)){
        col = "blue"
      }
    }
    v1 = d[d[,1] == e[i,1],]
    v2 = d[d[,1] == e[i,2],]
    if (e[i,5] == "MIG"){
      if (plotmig){
        #arrows( v1[1,]$x, v1[1,]$y, v2[1,]$x, v2[1,]$y, col = col, length = arrow, lwd = arrow_lwd, lty = arrow_lty)
        Arrows( v1[1,]$x, v1[1,]$y, v2[1,]$x, v2[1,]$y, lcol = col, 
                arr.length = arrow, arr.width = arrow, arr.type = "triangle", arr.lwd = 1, arr.adj = 1,
                lty = 1, lwd = arrow_lwd, segment = FALSE)
        
        Arrows( v1[1,]$x, v1[1,]$y, v2[1,]$x, v2[1,]$y, lcol = col, 
                arr.length = 0, arr.width = 0, arr.type = "triangle", arr.lwd = 1,
                lty = arrow_lty, lwd = arrow_lwd, segment = TRUE)
        
      }
    }
    else{
      lines( c(v1[1,]$x, v2[1,]$x), c(v1[1,]$y, v2[1,]$y), col = col, lwd = lwd)
    }
  }
  tmp = d[d[,5] == "TIP",]
  #print(tmp$x)
  #print(disp)
  if (!is.na(o)){
    for(i in 1:nrow(tmp)){
      tcol = o[o[,1] == tmp[i,2],2]
      if(plotnames){
        text(tmp[i,]$x+disp, tmp[i,]$y, labels = tmp[i,2], adj = 0, cex = cex, col  = tcol, font = font)
      }
    }
  }
  else{
    if (plotnames){
      shadowtext(tmp$x+disp, tmp$y, labels = tmp[,2], adj = 0, cex = cex, font = font, col = spoof_cols,
                 bg = "white", r=shadow)
      #text(tmp$x+disp, tmp$y, labels = tmp[,2], adj = 0, cex = cex, font = font, col = spoof_cols)
    }
  }
  if (scale){
    #print (paste("mse", mse))
    lines(c(0, mse*10), c(ybar, ybar))
    text( 0, ybar - 0.04, lab = "10 s.e.", adj = 0, cex  = 0.8)
    lines( c(0, 0), c( ybar - 0.01, ybar+0.01))
    lines( c(mse*10, mse*10), c(ybar- 0.01, ybar+ 0.01))
  }
  if (mbar){
    if (use_viridis){
      mcols = viridis(150)
    }else if (use_alpha){
      mcols = as.character(lapply(0:149, function(x) rgb(0, 0, 0, alpha = (x/149))))
    }else{
      mcols = rev( heat.colors(150) )
    }
    
    mcols = mcols[50:length(mcols)]
    ymi = ybar+0.15
    yma = ybar+0.35
    l = 0.2
    w = l/100
    xma = max(d$x/20)
    rect( rep(0, 100), ymi+(0:99)*w, rep(xma, 100), ymi+(1:100)*w, col = mcols, border = mcols)
    text(xma + 0.0001, ymi, lab = "0", adj = 0, cex = 0.7)
    if ( mw >0.5){ text(xma+disp, yma, lab = "1", adj = 0, cex = 0.7)}
    else{
      text(xma + 0.0001, yma, lab = "0.5", adj = 0, cex =0.7)
    }
    text(0, yma+0.09, lab = "Migration", adj = 0 , cex = 1)
    text(0, yma+0.06, lab = "weight", adj = 0 , cex = 1)
  }	
}

# from TeachingDemos

shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
         theta= seq(pi/4, 2*pi, length.out=500), r=0.05,... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('A')
  yo <- r*strheight('A')
  
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... ) }