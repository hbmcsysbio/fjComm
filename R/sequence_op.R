# read and write files

seq_gen_rnd <- function(seqLen, seqNum, prob=c("A"=0.25,"C"=0.25,"G"=0.25,"T"=0.25) )
{
  replicate(seqNum, paste(sample(c("A","T","C","G"),seqLen,replace = T, prob = prob),collapse = "" ) )
}


seq_dimer_sig_gen <- function(half_site_len=5, spacing=3,seqNum, dirct=c("rc","ori"))
{
  firstHalf= seq_gen_rnd(half_site_len,seqNum)
  if (dirct=="rc") return(paste0(firstHalf,seq_gen_rnd(spacing,seqNum),revComp(firstHalf)))
  if (dirct=="ori") return(paste0(firstHalf,seq_gen_rnd(spacing,seqNum),firstHalf))
}

seq_rnd_substitute <- function(seq= seq_gen_rnd(seqLen = 101,seqNum = 100000), sub= seq_dimer_sig_gen(4,3,100000), seqNum=100000)
{
  seqLen= nchar(seq[1])
  subLen= nchar(sub[1])
  subPos=sample(1:(seqLen-subLen+1),seqNum,replace = T)
  seq %>% stringr::`str_sub<-`(subPos,subPos+subLen-1,sub)
}


# seq=seq_gen_rnd(seqLen = 101,seqNum = 100000)
# sub=seq_dimer_sig_gen(half_site_len = 7,spacing = 3,100000,dirct = "rc")
# seq=seq_rnd_substitute(seq,sub)




        # getKmerTable <- function (fastaFile, kmerLength)
        # {
        #   # old kmer counter using qrqc
        #   library(pacman); pacman::p_load(qrqc)
        #   allseq<- readSeqFile(fastaFile,type = "fasta",k = kmerLength, hash.prop = 1)
        #   kmer=slot(allseq,"kmer")
        #   kmer=kmer[!grepl("N",kmer$kmer),]
        #   kmer
        # }
        #
        # # rm duplicated seq in file, using perl  //deprecated
        # dedup <- function (file_full,tmpDir=paste0(outDir,"/tmp/"),keepFile=FALSE)
        # {
        #     # rmdup_3point
        #     system(paste("mkdir -p ",tmpDir,sep = ""));
        #     cmd=(paste( "rmdup_3point.pl -f ",file_full," -k 20 -o ",tmpDir , sep = "")); system(cmd);
        #     file_base=basename(file_full)
        #     newName=sub("(.*)(\\..{1,})$", "\\1_u\\2", file_base)
        #     tmpFile= paste0(tmpDir,newName);
        #     if (!keepFile)
        #     {
        #       allreads= read_table(tmpFile,col_names = F)
        #       system(paste("rm ",tmpFile,sep = ""))
        #       return(allreads)
        #     }else
        #     {
        #        return(tmpFile)
        #     }
        # }

# reads in df, specify col containing reads, 2 20bp-wins as default, 4-5 best
rmdup <- function(seq_df, seq_col_num=1, winsize=20, winNo=2)
{
  # file and strings also OK
  if (typeof(seq_df)=="character"){ if (length(seq_df)==1){ seq_df=read_csv(seq_df,col_names = FALSE)} else {seq_df=data.frame(seq_df,stringsAsFactors = FALSE)}  }

  nrowIni=nrow(seq_df)
  seqlen=nchar(seq_df[[seq_col_num]][1])
  shift= floor(seqlen/winNo)
  for (i in (0:(winNo-1)))
  {
    start= shift*i+1; end=shift*i+winsize
    seq_df$sub= stringr::str_sub(string = seq_df[[seq_col_num]],start,end)
    seq_df= distinct(seq_df,sub,.keep_all = T)
  }

  nrowEnd=nrow(seq_df)
  print(paste0(nrowIni," input reads",", duplication rate is ",prettyNum((1-nrowEnd/nrowIni)*100, digits=4),"%"))
  seq_df$sub=NULL
  seq_df
}

rmdup_headtail<-function(seq_df, winsize=18)
{
  # file and strings also OK
  if (typeof(seq_df)=="character"){ if (length(seq_df)==1){ seq_df=read_csv(seq_df,col_names = FALSE)} else {seq_df=data.frame(seq_df,stringsAsFactors = FALSE)}  }

  nrowIni=nrow(seq_df)

  seq_df$sub= stringr::str_sub(seq_df[[1]],1,winsize)
  seq_df= distinct(seq_df,sub,.keep_all = T)

  seqlen=nchar(seq_df[[1]])
  seq_df$sub= stringr::str_sub(string = seq_df[[1]],-winsize,seqlen)
  seq_df= distinct(seq_df,sub,.keep_all = T)

  nrowEnd=nrow(seq_df)
  print(paste0(nrowIni," input reads",", duplication rate is ",prettyNum((1-nrowEnd/nrowIni)*100, digits=4),"%"))
  seq_df$sub=NULL
  seq_df
}

revComp <- function(str_array){Biostrings::DNAStringSet(str_array) %>% Biostrings::reverseComplement() %>% as.character()}

# only 1 instance for fw and rev, can remove palindromic also
fw_rev_dedup <- function(strArr,gapArr=NULL,rmPalin=F,return_Boolean=F)
{
  if (rmPalin){
    rev_k=revComp(strArr)
    strArr= strArr[strArr!=rev_k] # remove palindromic kmer
  }

  rev_k=revComp(strArr)
  ori_k_val=chartr("ACTG","0124",strArr) %>% as.integer() # !!!! sum of 2 should not equal to the other base, 0123 NG
  rev_k_val=chartr("ACTG","0124",rev_k) %>% as.integer()
  if(!is.null(gapArr)){dup= ifelse(ori_k_val<rev_k_val,strArr,rev_k) %>% paste0(gapArr) %>% duplicated()}
    else{dup= ifelse(ori_k_val<rev_k_val,strArr,rev_k) %>% duplicated()}
  if (return_Boolean) return(dup)
  strArr= strArr[!dup]
  strArr
}


# scoring RNA for Jilin
score_with_mono_di_pwm <- function(kcnt_df, monopwm, dipwm, stem=11, loop=0 )
{
  winNum= motifLen-kmerLen+1
  noPairLen_endpos= stem+loop
  noPairLen_startpos= stem+1
  tmpBkMat=matrix(log((0.25+0.00001)/bkfreq),nrow(kcnt),kmerLen)
  kcnt$scoreBk=rowSums(tmpBkMat)
  bkscoreVect= tmpBkMat[,1]

  kmerMask= rep(1,kmerLen)
  names(kmerMask)= paste0("k_",1:kmerLen)
  loopMask= c(rep(0,stem),rep(1,loop),rep(0,stem))
  names(loopMask)= c(paste0("s1_",1:stem),if(loop==0) character(0) else paste0("l_",1:loop), paste0("s2_",stem:1)) %>% na.omit()
  stemMask= !loopMask

  for (startpos in 1:winNum)
  {

    kmerMaskPos= c(rep(0,startpos-1),kmerMask,rep(0,motifLen-(startpos-1+kmerLen)))
    scoreCurrStartPos=numeric(nrow(kcnt)) # 0s, sum of the following
    loopScore=numeric(nrow(kcnt))
    unpairedSocre=numeric(nrow(kcnt))
    pairedScore=numeric(nrow(kcnt))

    # if (startpos==7) browser()

    if (any(loopMask & kmerMaskPos))
    {
      loop_pos_model= which(loopMask & kmerMaskPos) %>% names %>% stringr::str_split("_",simplify = T) %>% .[,2] %>% as.integer()
      loop_pos_kmer= which(kmerMaskPos & loopMask) %>% names %>% stringr::str_split("_",simplify = T) %>% .[,2] %>% as.integer()
      for(i in seq_along(loop_pos_model))
      {
        loopScore= loopScore + monopwm[match(kmerMat[,loop_pos_kmer[i]],rownames(monopwm)), loop_pos_model[i]]
      }
    }

    # tt= data.frame(kcnt$kmer,loopScore)

    # if (is.na(mean(loopScore))) browser()

    # stem, unpaired
    stem_model_left=(stemMask & kmerMaskPos) %>% .[1:stem]
    stem_model_right=(stemMask & kmerMaskPos) %>% .[-(1:(stem+loop))] %>% rev
    unpaired_bp_num= which(xor(stem_model_left,stem_model_right)) %>% length()
    unpairedSocre= unpairedSocre + bkscoreVect*unpaired_bp_num

    stem_paired_left= stem_model_left & stem_model_right
    stem_paired_right= stem_model_right & stem_model_left
    if(any(stem_paired_left))
    {
      paired_stem_mask= stemMask
      paired_stem_mask[1:stem]= stem_paired_left
      paired_stem_mask[-(1:(stem+loop))]= rev( stem_paired_right)

      stemKmerPos= names(which(kmerMaskPos & paired_stem_mask))
      stemModelPos= names(which(paired_stem_mask & kmerMaskPos))
      dict= stemKmerPos %>% stringr::str_split(.,"_",simplify = T)%>% .[,2] %>% as.integer()
      names(dict)= stemModelPos

      all_model_pos= names(dict) %>% stringr::str_split(.,"_",simplify = T)%>% .[,2] %>% unique() %>% as.integer()
      for (i in all_model_pos)
      {
        kmerPos1= dict[paste0("s1_",i)]
        kmerPos2= dict[paste0("s2_",i)]
        kmerPair= paste0(kmerMat[,kmerPos1],kmerMat[,kmerPos2])
        pairedScore= pairedScore + dipwm[match(kmerPair,rownames(dipwm)),i]
      }
    }

    print(paste(startpos,mean(loopScore),mean(unpairedSocre),mean(pairedScore) ))
    scoreCurrStartPos= loopScore + unpairedSocre + pairedScore
    kcnt[[paste0("pos",startpos)]]= scoreCurrStartPos
  }

  kcnt$maxScore= do.call(pmax, kcnt[,-c(1,2,3)])
  p=ggplot(data = kcnt,aes(log(counts),maxScore,kmer=kmer)) + geom_point() + theme_bw()
  print(p)
  kcnt
}

# count kmer form file
kmerCntBit_file<-function(file, k=2, diffLen=FALSE, collapse=FALSE, asDf= TRUE, all_possible_k=FALSE, pseudo=0,   rmdup=TRUE, barcode=NA, bkfile=NA)
{
  seq=read_csv(file,col_names = FALSE);
  if(!is.na(barcode)){ barLen=nchar(barcode); seq$bar=str_sub(seq[[1]],1,4); seq[[1]]=str_sub(seq[[1]],5,-1); seq=seq %>% dplyr::filter(bar==barcode) }
  if(rmdup) seq=rmdup(seq,1)
  result=kmerCntBit(seq[[1]], k=k, diffLen=diffLen, collapse=collapse, asDf= asDf, all_possible_k=all_possible_k, pseudo=pseudo)

  if(!is.na(bkfile))
  {
    if (!collapse) {collapse=TRUE; print("collapse because bkfile exist")}
    seq=read_csv(bkfile,col_names = FALSE);
    if(!is.na(barcode)){ barLen=nchar(barcode); seq$bar=str_sub(seq[[1]],1,4); seq[[1]]=str_sub(seq[[1]],5,-1); seq=seq %>% dplyr::filter(bar==barcode) }
    if(rmdup) seq=rmdup(seq,1)
    result1=kmerCntBit(seq[[1]], k=k, diffLen=diffLen, collapse=collapse, asDf= asDf, all_possible_k=all_possible_k, pseudo=pseudo)
    result %<>% mutate(bkcounts=result1$counts, log2chn=log2(counts/bkcounts))
  }
  result
}


gkmerCntBit_file<-function(file, gapNo=2, k=2, gapMins = c(0,0), gapMaxs = c(3,3), pseudo=5, diffLen=FALSE, posInfo=FALSE, all_possible_k=TRUE,    rmdup=TRUE, barcode=NA, melt_result=FALSE,  rc_combine=FALSE, rc_uniq=FALSE)
{
  # required logic relationship
  if (rc_uniq) rc_combine=TRUE
  if (rc_combine) melt_result=TRUE

  seq=read_csv(file,col_names = FALSE);
  if(!is.na(barcode)){ barLen=nchar(barcode); seq$bar=str_sub(seq[[1]],1,4); seq[[1]]=str_sub(seq[[1]],5,-1); seq=seq %>% dplyr::filter(bar==barcode) }
  if(rmdup) seq=rmdup(seq,1)
  result=gkmerCntBit(seq[[1]], gapNo=gapNo, k=k, gapMins = gapMins, gapMaxs = gapMaxs, pseudo=pseudo, diffLen=diffLen, posInfo=posInfo, all_possible_k=all_possible_k)
  if(melt_result) result=reshape2::melt(result) %>% set_colnames(c("kmer","gap","counts")) %>% mutate_if(is.factor,as.character)

  if(rc_combine)
  {
    rc_result= result %>% mutate(kmer=revComp(kmer))
    result %<>% mutate(counts= counts + rc_result[match(paste0(kmer,gap),paste0(rc_result$kmer,rc_result$gap)),]$counts)
    rm(rc_result)
  }

  if(rc_uniq)
  {
      dupped= fw_rev_dedup(result$kmer,result$gap,return_Boolean = T)
    result %<>% dplyr::filter(!dupped)
  }
  result
}



