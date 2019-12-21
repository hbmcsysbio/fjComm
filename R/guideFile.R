# df should with "pos" in the 1st column
addColtoGuide<-function(guidefile,df)
{
  guide=read.csv(guidefile,sep = "\t")

  # first remove empty row in guide
  guide=guide[(!guide$pos=="") & (!is.na(guide$pos)),] # & is elementwise while && only calc 1st element of a vector

  # order guide
  corrOrderGuide= order(str_sub(guide$pos,1,1), as.numeric(str_sub(guide$pos,2,-1))  ) #sort first with A-Z, then 1-11
  guide=guide[corrOrderGuide,]

  # order df
  corrOrder= order(str_sub(df$pos,1,1), as.numeric(str_sub(df$pos,2,-1))  ) #sort first with A-Z, then 1-11
  df=df[corrOrder,] #order and remove the pos column

  # guide=cbind(guide,df[,-1])
  guide=cbind(guide,df[!(colnames(df)=="pos")])
  colnames(guide)[(length(colnames(guide))-ncol(df)+2):length(colnames(guide)) ]= colnames(df[!(colnames(df)=="pos")]) #colnames(df)[2:length(colnames(df))] # correct name
  write.table(guide, file = guidefile, sep = "\t", quote = F, row.names = F)
  return("edit guide success")
}


guide.get_Info_by_well <- function(wellVect, guideFilePath, fields=c("TF","family","class","pos"))
{
  guide=readr::read_tsv(guideFilePath)
  return(guide[match(wellVect,guide$pos),][fields])
}
