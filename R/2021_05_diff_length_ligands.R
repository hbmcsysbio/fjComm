
length_adjust<-function(seqs,output_length=40L,seq_in_middle=FALSE)
  # adjust all ligands to the same length, out put seq only
{
  seqs %<>% str_sub(1,output_length)
  Ns=output_length-nchar(seqs)
  topaste=strrep("N",Ns)
  if(seq_in_middle){
    leftNs=(Ns/2) %>% as.integer()
    rightNs=Ns-leftNs
    leftNs=strrep("N",leftNs)
    rightNs=strrep("N",rightNs)
    seqs=paste0(leftNs,seqs,rightNs)
  }else seqs=paste0(seqs,topaste)
  seqs
}
