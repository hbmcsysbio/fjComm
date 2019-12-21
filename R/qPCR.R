
# calc ct points
qPCR.get_d2max <- function(qPCR_file, interpolate_points=100, smooth_bdwidth=2, threshold=3) # if all point fluo < threshold, call as "NA"
{
  pacman::p_load(features)

  qData= readr::read_tsv(qPCR_file,col_names = T,skip = 1)
  colnames(qData)= c("well","name","prog","seg","cyc","Time","Temp","FL","notUsed")
  qData= dplyr::filter(qData, prog==2) %>% reshape2::acast(well ~cyc, mean, value.var = "FL") # prog=2 is the qPCR amp, 3 for melting curve
  # results <- apply(qData, 1, fjComm::analyse,  base.line = "slanted", output = "all", silent = TRUE) %>% t # tolerance is important to deal with bad data

  corrOrder= order(stringr::str_sub(rownames(qData),1,1), as.integer(stringr::str_sub(rownames(qData),2,-1))  ) #sort first with A-Z, then 1-11
  qData=qData[corrOrder,]

  d2max=numeric(nrow(qData))
  x= as.integer(colnames(qData))
  for (i in seq_len(nrow(qData)))
  {
    if (max(qData[i,])< threshold) { d2max_tmp=NA }
    else{
      interpolated= approx(x, qData[i,], n = interpolate_points)
      fitted=features(interpolated$x,interpolated$y,bandwidth= smooth_bdwidth) %>% attr("fits")
      d2max_tmp= fitted$x[which.max(fitted$d2)]
    }
    d2max[i]=d2max_tmp
  }
  names(d2max)=rownames(qData)
  return(d2max)
}

qPCR.getData <- function(qPCR_file)
{
  pacman::p_load(features)

  qData= readr::read_tsv(qPCR_file,col_names = T,skip = 1)
  colnames(qData)= c("well","name","prog","seg","cyc","Time","Temp","FL","notUsed")
  qData= dplyr::filter(qData, prog==2) %>% reshape2::acast(well ~cyc, mean, value.var = "FL") # prog=2 is the qPCR amp, 3 for melting curve
  # results <- apply(qData, 1, fjComm::analyse,  base.line = "slanted", output = "all", silent = TRUE) %>% t # tolerance is important to deal with bad data

  corrOrder= order(stringr::str_sub(rownames(qData),1,1), as.integer(stringr::str_sub(rownames(qData),2,-1))  ) #sort first with A-Z, then 1-11
  qData=qData[corrOrder,]
  return(qData)
}
