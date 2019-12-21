
# return Vect with corrected baseline
math_sub_baseline <-function(Vect,baseline_hwm=5)
{
  pacman::p_load(baseline)
  baseline.medianWindow(matrix(Vect,nrow = 1),hwm = baseline_hwm)$corrected[1,]
}
math_get_baseline <-function(Vect,baseline_hwm=5)
{
  pacman::p_load(baseline)
  baseline.medianWindow(matrix(Vect,nrow = 1),hwm = baseline_hwm)$baseline[1,]
}


# squish and normalize to 0,1
  squishQuan_normalize <- function(vect,quan=c(0.005,0.995),  normalize=TRUE)
  {
    vect1= vect %>% {scales::squish(.,quantile(.,c(0.005,0.995)))}
    if(normalize) vect1%<>% {.- min(.)} %>% {./max(.)}
    vect1
  }
#!!!-..12


# p value from winflat
winflat<-function(xvalue=100,yvalue=99,sizex=100000,sizey=100000)
{
  cmd= glue::glue("ssh yin@172.25.122.75 winflat -xvalue {xvalue} -yvalue {yvalue} -diff {sizex} {sizey}")
  options(warn=-1)
  result=system(cmd,intern = T)
  options(warn=0)
  result %>% str_extract("(?<==)[^=]*$") %>% as.numeric() %>% min()
}
winflat_v=Vectorize(winflat)

