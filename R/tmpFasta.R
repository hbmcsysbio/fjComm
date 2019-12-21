# convert Just reads to fasta and store in ./tmp
# source(paste(RlibDir,"tmpFasta.R",sep = ""))

#return path of tmp fasta file
genTmpFasta <- function (JustReadsFile_full,lines=NULL,gz=NULL,tmpFileName=NULL)
{
  cat1 = ifelse(!is.null(gz), "zgrep . ", "cat")
  headn= ifelse(!is.null(lines), paste(" | head -n ", lines, sep = ""), "")
  if(!is.null(tmpFileName)) {tmpFastaName=tmpFileName} else {tmpFastaName= paste("./tmp_",basename(JustReadsFile_full),".fa",sep = "")}
  system(paste(cat1," ",JustReadsFile_full, " ", headn, " | perl -lpe '",'$_=">\n".',"$_' > ",tmpFastaName, sep = "" ) )
  return(tmpFastaName)
}

rmTmpFasta <- function (tmpFastaName)
{
  system( paste("rm ",tmpFastaName, sep = "") )
}

kmer_table_JustReads <- function(JustReadsFile,lines=NULL,gz=NULL,tmpFileName=NULL,kmerLen=4L,aggregate=FALSE, posRng=NULL)
{
  tempFasta=genTmpFasta(JustReadsFile,lines,gz,tmpFileName)

  library(qrqc); library(dplyr)
  allseq<- readSeqFile(tempFasta,type = "fasta",k = kmerLen, hash.prop = 1)
  kmer=slot(allseq,"kmer")
  kmer=data.table(kmer) %>% dplyr::filter(!grepl("N",kmer))
  if(!is.null(posRng)) {kmer=subset(kmer, position %in% posRng)}
  if (aggregate==T) { kmer= aggregate(count~kmer, data = kmer, sum) }

  rmTmpFasta(tempFasta)
  return(kmer)
}
