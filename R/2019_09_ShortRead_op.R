# calc duplication rate of sequences
sr_dup_rate<-function(ShortReadObj,isfqfileName=FALSE)
{
  if(isfqfileName) ShortReadObj=ShortRead::readFastq(ShortReadObj)
  ShortRead::srduplicated(ShortReadObj) %>% table() %>% {.["TRUE"]/(.["TRUE"]+.["FALSE"])} %>% set_names(NULL) %>% {sprintf("%.2f",.*100)} %>% as.numeric()
}



getSeq_fqfachrFile<-function(file)
  #getSeq from .fq and .fa or seqonly files
{
  pacman::p_load(ShortRead)
  head_line=readLines(file,n = 1)
  if (grepl("@",head_line)) seqs=readFastq(file) %>% {.@sread} %>% as.character() else {
    if(grepl(">",head_line))
    {
      seqs=readFasta(file) %>% {.@sread} %>% as.character()
    }else seqs=readLines(file)
  }
  seqs
}

