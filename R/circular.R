# circular dendrogram with heat map

plotCircle<- function(
  mat = matrix(rnorm(100*20), nrow = 10, ncol = 50), # each col is a TF
  col_names_mat=  rep("testTF", 50), # TF names of each col
  factors = rep(letters[1], ncol(mat)), # grouping info of each col (class,family...)
  f = colorRamp2(breaks = c(-1, 0, 1), colors = c("white", "steelblue", "red"))
)
{
   pacman::p_load(circlize)

  # dend: as dendrogram object, usually returned by hclust and as.dendrogram
  # maxy: maximum height of the tree
  circos.dendrogram = function(dend, maxy) {
    labels = as.character(labels(dend))
    x = seq_along(labels) - 0.5
    names(x) = labels

    is.leaf = function(object) (is.logical(L <- attr(object, "leaf"))) && L

    # recursive function to draw the tree
    draw.d = function(dend, maxy) {
      leaf = attr(dend, "leaf")
      d1 = dend[[1]]
      d2 = dend[[2]]
      height = attr(dend, 'height')
      midpoint = attr(dend, 'midpoint')

      if(is.leaf(d1)) {
        x1 = x[as.character(attr(d1, "label"))]
      } else {
        x1 = attr(d1, "midpoint") + x[as.character(labels(d1))[1]]
      }
      y1 = attr(d1, "height")

      if(is.leaf(d2)) {
        x2 = x[as.character(attr(d2, "label"))]
      } else {
        x2 = attr(d2, "midpoint") + x[as.character(labels(d2))[1]]
      }
      y2 = attr(d2, "height")

      circos.lines(c(x1, x1), maxy - c(y1, height), straight = TRUE)
      circos.lines(c(x1, x2), maxy - c(height, height))
      circos.lines(c(x2, x2), maxy - c(y2, height), straight = TRUE)

      if(!is.leaf(d1)) {
        draw.d(d1, maxy)
      }
      if(!is.leaf(d2)) {
        draw.d(d2, maxy)
      }
    }

    draw.d(dend, maxy)
  }

  # factors = rep(letters[1:2], 100)


  par(mar = c(1, 1, 1, 1))
  circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5) # gap between 2 sectors
  circos.initialize(factors, xlim = c(0, ncol(mat)))
  maxy = 0 # height from the heatplot




  # heat map
  circos.trackPlotRegion(ylim = c(0, nrow(mat)), bg.border = NA,track.height = 0.2, panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    m = mat[, factors == sector.index]

    dend.col = as.dendrogram(hclust(dist(t(m))))

    maxy = ifelse(maxy > attr(dend.col, "height"), maxy, attr(dend.col, "height"))
    assign("maxy", maxy, envir = parent.env(environment()))

    m2 = m[, labels(dend.col)]
    nr = nrow(m2)
    nc = ncol(m2)
    for(i in 1:nr) {
      for(j in 1:nc) {
        circos.rect(j-1, nr-i, j, nr-i+1, border = f(m2[i, j]), col = f(m2[i, j]))
      }
    }

  })


  # text track
  circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.1, panel.fun = function(x, y) {
    for(i in seq_len(length(col_names_mat))) {
      circos.text(i-0.5, 0, col_names_mat[i], adj = c(1, 0.5),
                  facing = "reverse.clockwise", niceFacing = TRUE,
                  col = 1, cex = 0.7)
    }
  })


  # dendrogram , track.height is the height of the tree
  circos.trackPlotRegion(ylim = c(0, maxy), bg.border = NA, track.height = 0.6, panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    m = mat[, factors == sector.index]

    dend.col = as.dendrogram(hclust(dist(t(m))))

    circos.dendrogram(dend.col, maxy)

  })

  circos.clear()

  ## legend at the center
  # x = seq(-10, 10, length.out=100)/40
  # col =f(seq(-2, 2, length.out = length(x-1)))
  # for(i in seq_along(x)) {
  #   if(i == 1) next
  #   rect(x[i-1], -0.05, x[i], 0.05, col = col[i], border = col[i])
  # }
  #
  # text(x[1], -0.08, "-2", adj = c(0.5, 1), cex = 1.2)
  # text(x[ceiling(length(x)/2)], -0.08, "0", adj = c(0.5, 1), cex = 1.2)
  # text(x[length(x)], -0.08, "2", adj = c(0.5, 1), cex = 1.2)
}
