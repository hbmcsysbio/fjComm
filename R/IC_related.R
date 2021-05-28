recoverDf_from_return_list <- function(returnList)
{
  colNum=returnList[[1]] %>% length()
  lapply(returnList,function(x){x[[1]]})
  rearr_list= vector("list",4)
  names(rearr_list)= returnList[[1]] %>% names
  for (i in seq_len(colNum)){rearr_list[[i]]=lapply(returnList,function(x){x[[i]]}) %>% unlist}
  rearr_list %>% as.data.frame()
}

ic_get_2vect_cnt_table <- function(V1,V2,pseudo)
{
  counts=  ( data.table(V1,V2) %>% .[,.N,by=list(V1,V2)] %>% with(., xtabs(N~V1+V2)) ) + pseudo
  counts[!grepl("N",rownames(counts)),!grepl("N",colnames(counts)), drop=F]
}




    # all arrange according to foldchn, but return diff value
    maxBiascalc_elemental <- function (seqMat,col1,col2,topNo,pseudo,type=c("freq","foldchn","topMI"),dimer=FALSE,head_to_tail=TRUE,head_to_head=TRUE, directed=FALSE) # only hh true for RNA
    {
      if(length(type)>1 || (!(type %in% c("freq","foldchn","topMI"))) ) stop("maxBiascalc_elemental: please specify type of calc!!!")

      actCnt= ic_get_2vect_cnt_table(V1 = seqMat[,col1], V2 = seqMat[,col2], pseudo = pseudo)
      # browser()
        if(length(actCnt)==0) return(0)

        actFreq= actCnt %>% prop.table

        expFreq= margin.table(actFreq,1) %o% margin.table(actFreq,2) %>% as.matrix() %>% melt
        actFreq= melt(actFreq)#; actCnt= melt(actCnt)

        if (dimer)
        {
          e=parent.env(parent.frame())
          if(is.null(e$dimerRows)) # figure out dimer rows
          {
            kmer1=actFreq$V1 %>% as.character()
            kmer2=actFreq$V2 %>% as.character()
            if(head_to_head && head_to_tail) {e$dimerRows= (kmer1==kmer2) | (kmer1==revComp(kmer2))}
            else if(head_to_tail) {e$dimerRows= (kmer1==kmer2)}
            else if(head_to_head) {e$dimerRows= (kmer1==revComp(kmer2))}
            else {stop("error, should specify valid dimer mode in maxBiascalc_elemental()")}
          }
          actFreq= actFreq[e$dimerRows,]
          expFreq= expFreq[e$dimerRows,]
        }

        actFreq= actFreq %>% mutate(foldchn= actFreq$value/ expFreq$value )%>% arrange(foldchn) # only foldchn can be less stable if enrichment is high (value*)
        if (type=="foldchn")
        {
          bottomDevMean= head(actFreq,topNo)$foldchn %>% mean
          topDevMean= tail(actFreq,topNo)$foldchn %>% mean
          topk1=tail(actFreq,1)$V1 %>% as.character()
          topk2=tail(actFreq,1)$V2 %>% as.character()
          return(list(foldchnTop=topDevMean,foldchnEnd=bottomDevMean,topk1=topk1,topk2=topk2))
        }

        if (type=="freq")
        {
          bottomDevMean= head(actFreq,topNo)$value %>% mean
          topDevMean= tail(actFreq,topNo)$value %>% mean
          topk1=tail(actFreq,1)$V1 %>% as.character()
          topk2=tail(actFreq,1)$V2 %>% as.character()
          return(list(freqTop=topDevMean,freqEnd=bottomDevMean,topk1=topk1,topk2=topk2))
        }

        if (type=="topMI")
        {
          actFreq=rbind(head(actFreq,topNo),tail(actFreq,topNo))
          actFreq= actFreq %>% mutate(MI= actFreq$value * log2(foldchn) ) #%>% arrange(MI)
          # browser()
            # only take into acc MI from 1 direction of the kmer
            if(directed){
              allk=with(actFreq, paste0(V1,V2))
              allkrc=revComp(allk)
              k_filter=allkrc<=allk
              actFreq=actFreq %>% mutate(MI=ifelse(k_filter,MI,0))
            }

          bottomMIsum= head(actFreq,topNo)$MI %>% sum
          topMIsum= tail(actFreq,topNo)$MI %>% sum
          topk1=tail(actFreq,1)$V1 %>% as.character()
          topk2=tail(actFreq,1)$V2 %>% as.character()
          return(list(topMIsum=topMIsum,bottomMIsum=bottomMIsum,topk1=topk1,topk2=topk2))
        }
    }


    MIcalc_elemental <- function (seqMat,col1,col2,punish,pseudo)
    {
      counts= ic_get_2vect_cnt_table(V1 = seqMat[,col1], V2 = seqMat[,col2], pseudo = pseudo)
      if(length(counts)==0) return(0+punish)

      mi.empirical(counts,unit ="log2")+punish
    }

    # if filter for spacing is T, only cols with the specified "spacing" is retained, spacing = 0 for neighbouring kmers
ic_related_calc <- function (seqs=character(0), kmerLen=2L, filter_for_spacing=TRUE, spacing=c(0L,1L,2L,3L), verbose=F, pseudo=10L, type=c("MI","maxBias","dimer","RNA","direct_maxBias"),
                             maxBias_dimer_Params=list(type="topMI",topNo=5L) ) # "foldchn", "topMI", "freq"
{
  if(length(type)>1 || (!(type %in% c("MI","maxBias","dimer","RNA","direct_maxBias"))) ) stop("ic_related_calc: please specify type of calc!!!")
  pacman::p_load(entropy)
  seqMat= seqFregments(seqs,k=kmerLen)

  # gen pseudo cnt mat
  eachCol= gtools::permutations(4,kmerLen,c("A","T","C","G"),repeats.allowed = T) %>% apply(1,function(x){paste0(x,collapse="")})
  # pseudoMat= gtools::permutations(4^kmerLen, 2, eachCol, repeats.allowed = T)
  # pseudoMat= table(pseudoMat[,1],pseudoMat[,2]) * pseudo

  # gen all position combinations, with empty row as score
  kmerNo=ncol(seqMat)
  resultDf= utils::combn(c(1:kmerNo),2) %>% t %>% as.data.frame()
  resultDf=resultDf %>% dplyr::filter(!abs(V2-V1)< kmerLen)
  if (filter_for_spacing) resultDf= resultDf %>%  dplyr::filter(abs(V2-V1) %in% (kmerLen+spacing)) #== (kmerLen+spacing))
  colnames(resultDf)=c("pos1","pos2")
  punish=0 # -(4^(kmerLen*2)-1)/(2*log(2)*nrow(seqMat)) + 2* ((4^(kmerLen)-1)/(2*log(2)*nrow(seqMat)))
  seqMat= seqMat %>% as.data.frame() %>% mutate_all(list(factor),levels=eachCol) # to factor, add all possible levels

    # prog bar
    pb <- utils::txtProgressBar(min = 0, max = nrow(resultDf), style = 3); pbCnt=0;
    updateProgressBar <- function(pb){ pbCnt<<- pbCnt+1; Sys.sleep(0.0001); setTxtProgressBar(pb, pbCnt)}
  if (type=="MI")
  {
    resultDf$MI= apply(resultDf,1,function(x){   # calc pair-wise
      if(verbose) {cat(paste0(x[1],"-",x[2],"  "))}
      updateProgressBar(pb)
      MIcalc_elemental(seqMat, x[1],x[2],punish, pseudo)
    })
    return(resultDf)
  }

  if (type=="maxBias")
  {
    calcResult= apply(resultDf,1,function(x){   # calc pair-wise
      if(verbose) { cat(paste0(x[1],"-",x[2],"  "))}
      updateProgressBar(pb)
      maxBiascalc_elemental(seqMat, x[1],x[2], topNo = maxBias_dimer_Params$topNo, pseudo = pseudo, type = maxBias_dimer_Params$type)
    })
    calcResult=recoverDf_from_return_list(calcResult)
    resultDf= cbind(resultDf,calcResult)
    return(resultDf)
  }

  if (type=="direct_maxBias")
  {
    calcResult= apply(resultDf,1,function(x){   # calc pair-wise
      if(verbose) { cat(paste0(x[1],"-",x[2],"  "))}
      updateProgressBar(pb)
      maxBiascalc_elemental(seqMat, x[1],x[2], topNo = maxBias_dimer_Params$topNo, pseudo = pseudo, type = maxBias_dimer_Params$type, directed=TRUE)
    })
    calcResult=recoverDf_from_return_list(calcResult)
    resultDf= cbind(resultDf,calcResult)
    return(resultDf)
  }

  if (type=="dimer")
  {
    dimerRows= NULL
    calcResult= apply(resultDf,1,function(x){   # calc pair-wise
      if(verbose) { cat(paste0(x[1],"-",x[2],"  "))}
      updateProgressBar(pb)
      maxBiascalc_elemental(seqMat, x[1],x[2], topNo = maxBias_dimer_Params$topNo, pseudo = pseudo, type = maxBias_dimer_Params$type, dimer=TRUE)
    })
    calcResult=recoverDf_from_return_list(calcResult)
    resultDf= cbind(resultDf,calcResult)
    return(resultDf)
  }

  if (type=="RNA")
  {
    dimerRows= NULL
    calcResult= apply(resultDf,1,function(x){   # calc pair-wise
      if(verbose) { cat(paste0(x[1],"-",x[2],"  "))}
      updateProgressBar(pb)
      maxBiascalc_elemental(seqMat, x[1],x[2], topNo = maxBias_dimer_Params$topNo, pseudo = pseudo, type = maxBias_dimer_Params$type, dimer=TRUE, head_to_tail =FALSE ,head_to_head =TRUE )
    })
    calcResult=recoverDf_from_return_list(calcResult)
    resultDf= cbind(resultDf,calcResult)
    return(resultDf)
  }

  stop("calc type not valid: ic_related_calc()")

}



IC <- function (vector,punish_single)
{
  cntTab=table(vector[!grepl("N",vector)])
  propTab=prop.table(cntTab)
  2*k_1posInfo+sum(propTab*(log2(propTab)))-punish_single
}

len4_vect_IC <- function(cntVect=c(1,2,2,4))
  # table of "A","C","T","G"
{
  cntVect=c(cntVect,rep(0.001,4-length(cntVect)))
  cntVect %>% prop.table() %>% {sum(.*log2(./0.25))}
}





# file="~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_200/Trulig200v1IIIFJ4-4-CAP2-P10-TF377-TIIIc4_S754_L002_R1_001.peared_trimmed.fq_u.gz"
#
# seq= readr::read_csv(file,col_names = F)#,n_max = 100000)
# seq=rmdup(seq,1)
# result_topMI=ic_related_calc(seqs = seq[[1]],kmerLen = 3,filter_for_spacing = T,spacing = 0:10,type = "maxBias",maxBias_dimer_Params=list(type="topMI",topNo=10L))
# print(gg_heat2D_MI(result_topMI))

