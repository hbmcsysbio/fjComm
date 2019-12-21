# take txt or gz file
kmerCnt_allgap_file <- function(file, k=2L, minGap=0L, maxGap=0L, asDf=TRUE, diffLen=FALSE)
{
  file_cnt=read_csv(file)
  return(kmerCnt_allgap(file_cnt[[1]], k, maxGap, asDf, diffLen, minGap))
}
