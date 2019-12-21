# calc duplication rate of sequences
sr_dup_rate<-function(ShortReadObj,isfqfileName=FALSE)
{
  if(isfqfileName) ShortReadObj=ShortRead::readFastq(ShortReadObj)
  ShortRead::srduplicated(ShortReadObj) %>% table() %>% {.["TRUE"]/(.["TRUE"]+.["FALSE"])} %>% set_names(NULL) %>% {sprintf("%.2f",.*100)} %>% as.numeric()
}
