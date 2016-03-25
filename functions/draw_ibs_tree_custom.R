draw_ibs_tree_custom <- function (obj, clust.count = NULL, dend.idx = NULL, type = c("dendrogram", 
                                                             "z-score"), yaxis.height = TRUE, yaxis.kinship = TRUE, y.kinship.baseline = NaN, 
          y.label.kinship = FALSE, outlier.n = NULL, shadow.col = c(rgb(0.5, 
                                                                        0.5, 0.5, 0.25), rgb(0.5, 0.5, 0.5, 0.05)), outlier.col = rgb(1, 
                                                                                                                                      0.5, 0.5, 0.5), leaflab = "none", labels = NULL, y.label = 0.2, 
          ...) 
{
  stopifnot(is.null(dend.idx) | is.numeric(dend.idx))
  type <- match.arg(type)
  stopifnot(is.numeric(y.kinship.baseline))
  if (type == "dendrogram") {
    stopifnot(!is.null(obj$dendrogram))
    stopifnot(is.null(outlier.n) | is.numeric(outlier.n))
    if (is.null(clust.count)) 
      clust.count <- obj$clust.count
    if (is.null(outlier.n)) 
      outlier.n <- obj$outlier.n
    if (!is.null(dend.idx)) {
      den <- obj$dendrogram[[dend.idx]]
      x.offset <- 0
      for (i in 1:length(dend.idx)) {
        if (dend.idx[i] == 2) {
          IX <- dend.idx[1:i]
          IX[i] <- 1
          x.offset <- x.offset + attr(obj$dendrogram[[IX]], 
                                      "member")
        }
      }
    }
    else {
      den <- obj$dendrogram
      x.offset <- 0
    }
    par(mar = c(4, 4, 4, 4))
    oldpar <- par(mgp = c(5, 1, 0))
    plot(den, leaflab = leaflab, axes = FALSE, ...)
    par(oldpar)
    if (yaxis.height) {
      axis(side = 2, line = 0)
      tmp <- list(...)
      if (!is.null(tmp$ylab)) 
        ylab <- tmp$ylab
      else ylab <- "individual dissimilarity"
      mtext(ylab, side = 2, line = 2.5)
    }
    if (yaxis.kinship) {
      if (is.finite(y.kinship.baseline)) {
        y.kinship.baseline <- y.kinship.baseline[1]
      }
      else {
        y.kinship.baseline <- attr(den, "height")
      }
      ym <- pretty(c(0, 1))
      axis(side = 4, (1 - ym) * y.kinship.baseline, ym, 
           line = 0)
      mtext("coancestry coefficient", 4, line = 2.5)
    }
    if (!is.null(clust.count)) {
      m <- c(0, cumsum(clust.count))
      jj <- 1
      k <- 1
      for (i in 1:length(clust.count)) {
        if (clust.count[i] > outlier.n) {
          rect(m[i] + 0.5 - x.offset, par("usr")[3L], 
               m[i + 1] + 0.5 - x.offset, par("usr")[4L], 
               col = shadow.col[jj], border = NA)
          jj <- 3 - jj
          if (!is.null(labels[k])) 
            text((m[i] + m[i + 1])/2 - x.offset, y.label, 
                 labels[k])
          k <- k + 1
        }
        else {
          rect(m[i] + 0.5 - x.offset, par("usr")[3L], 
               m[i + 1] + 0.5 - x.offset, par("usr")[4L], 
               col = outlier.col, border = NA)
        }
      }
    }
    if (yaxis.kinship & y.label.kinship) {
      h1 <- (1 - 0.5) * y.kinship.baseline
      abline(h = h1, lty = 2, col = "gray")
      h2 <- (1 - 0.25) * y.kinship.baseline
      abline(h = h2, lty = 2, col = "gray")
      h3 <- (1 - 1/8) * y.kinship.baseline
      abline(h = h3, lty = 2, col = "gray")
      h4 <- (1 - 1/16) * y.kinship.baseline
      abline(h = h4, lty = 2, col = "gray")
      axis(side = 4, c(h1, h2, h3, h4), c("twins", "PC/FS", 
                                          "DFC/HS", "FC"), tick = FALSE, line = -0.75, 
           las = 2, cex.axis = 0.75, col.axis = "gray25")
    }
  }
  else if (type == "z-score") {
    if (is.null(obj$merge)) 
      stop("There is no Z score in this object.")
    y <- obj$merge[, 1]
    y <- y[order(y, decreasing = TRUE)]
    plot(y, xlab = "the order of Z score", ylab = "Z score", 
         type = "b", pch = "+", log = "x", ...)
    abline(h = 15, col = "gray", lty = 2)
  }
  invisible()
}
