.onLoad <- function(libname, pkgname){
  # geom settings updated
  ggplot2::update_geom_defaults("point",list(size = 0.6,alpha=0.99))
  ggplot2::theme_set(gg_theme_Publication())
  if(!is.element("SeqTools", installed.packages()[,1])) devtools::install_github("lianos/seqtools/R/pkg")
  if(!is.element("ggseqlogo", installed.packages()[,1])) pacman::p_load(ggseqlogo)
}
