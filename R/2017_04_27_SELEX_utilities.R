# gen 384 well lable
# gen_384_label <- function() expand.grid(letters[1:16] %>% toupper(), 1:24) %>% arrange(Var1,Var2) %>% {paste0(.[[1]],.[[2]])}
gen_96_label <- function(byrow=T) expand.grid(letters[1:8] %>% toupper(), 1:12) %>% {paste0(.[[1]],.[[2]])} %>% matrix(nrow=12, byrow = byrow) %>% as.character()
gen_384_label <- function(byrow=T) expand.grid(letters[1:16] %>% toupper(), 1:24) %>% {paste0(.[[1]],.[[2]])} %>% matrix(nrow=24, byrow = byrow) %>% as.character()

# if test=T do not rmdup, and take head 1000 line
xyplot <- function(xfile, yfile, test=F, kmerLen=3L, rm_dup=F, gapped=T, gapNo = 2, gapMins = c(0,0), gapMaxs = c(5,5), topKmersEach= 3000,alpha=0.2,size=0.2)
{

  if (!gapped)
  {
    y_selex=  {y_selex= SELEXFile(yfile); y_selex$getSeq(dup_rm = rm_dup); if(exists("test") && test)y_selex$seq= head(y_selex$seq,1000); y_selex$count_k(k = kmerLen,collapse = T,asDf = T,pseudo = 10); y_selex}
    x_selex=  {x_selex= SELEXFile(xfile); x_selex$getSeq(dup_rm = rm_dup); if(exists("test") && test)x_selex$seq= head(x_selex$seq,1000); x_selex$count_k(k = kmerLen,collapse = T,asDf = T,pseudo = 10); x_selex}
    kktable=y_selex$kmerCnt %>% rename(ycounts=counts) %>% mutate(xcounts=x_selex$kmerCnt$counts); rm(x_selex,y_selex)
    kktable=dplyr::union(kktable %>% top_n(topKmersEach*10,ycounts), kktable %>% top_n(topKmersEach*10,xcounts))
    topkmers= kktable %>% top_n(topKmersEach, rank(ycounts)+rank(xcounts))

      topkmer_y= kktable %>% top_n(topKmersEach, rank(ycounts))
      topkmer_x= kktable %>% top_n(topKmersEach, rank(xcounts))
    topkmers=dplyr::union(topkmers,topkmer_y,topkmer_x)

    corr=cor(topkmers$ycounts,topkmers$xcounts,method = "pearson")%>% prettyNum(digit=2)
    topkmers %<>% mutate(cntA=str_count(kmer,"A"), cntT=str_count(kmer,"T"),cntC=str_count(kmer,"C"),cntG=str_count(kmer,"G"))
    maxCnt=topkmers$kmer[1] %>% nchar()
    topkmers %<>% mutate_at(c("cntA","cntT","cntC","cntG"),function(x) floor(x/maxCnt*256))
    topkmers %<>% mutate(Rval=(cntT+cntG) ,Gval=(cntA+cntG) ,Bval=(cntC) ,pcolor=rgb(red = Rval,green = Gval,blue = Bval,maxColorValue = 256))

    topkmers=topkmers[sample(1:nrow(topkmers)),]
    xyplot=ggplot() + geom_point(data=topkmers,aes(xcounts,ycounts,kmer=kmer,color=I(pcolor)),alpha=alpha,size=size)+ gg_theme_Publication() + gg_anno_grob(paste0("pearson corr: ", corr))+
      gg_anno_grob("A","darkgreen",0.95,0.5)+
      gg_anno_grob("C","blue",0.95,0.45)+
      gg_anno_grob("G","darkgoldenrod",0.95,0.4)+
      gg_anno_grob("T","red",0.95,0.35)


  }else{
    y_selex= {
      y_selex= SELEXFile(yfile); y_selex$getSeq(dup_rm = rm_dup );
      if(test)y_selex$seq= head(y_selex$seq,1000);
      y_selex$count_gpk(k = kmerLen, gapNo = gapNo,gapMins = gapMins, gapMaxs = gapMaxs);
      y_selex$kmerCntGap= y_selex$kmerCntGap %>% mutate(kmer=paste0(kmer,"_",gapLen)) %>% select(-gapLen)
      y_selex
    }
    x_selex= {
      x_selex= SELEXFile(xfile); x_selex$getSeq(dup_rm = rm_dup);
      if(exists("test") && test)x_selex$seq= head(x_selex$seq,1000);
      x_selex$count_gpk(k = kmerLen,  gapNo = gapNo,gapMins = gapMins, gapMaxs = gapMaxs);
      x_selex$kmerCntGap= x_selex$kmerCntGap %>% mutate(kmer=paste0(kmer,"_",gapLen)) %>% select(-gapLen)
      x_selex
    }

    kktable=y_selex$kmerCntGap %>% dplyr::rename(ycounts=counts); kktable$xcounts=x_selex$kmerCntGap$counts; rm(x_selex,y_selex)
    kktable=dplyr::union(kktable %>% top_n(topKmersEach*10,ycounts), kktable %>% top_n(topKmersEach*10,xcounts))
    topkmers_both= kktable %>% top_n(topKmersEach, rank(ycounts)+rank(xcounts))

      topkmer_y= kktable %>% top_n(topKmersEach, rank(ycounts))
      topkmer_x= kktable %>% top_n(topKmersEach, rank(xcounts))
    topkmers=dplyr::union(topkmer_y,topkmer_x) %>% dplyr::union(topkmers_both)
      # topkmers=dplyr::union(topkmer_y,topkmer_x)

    topkmers= topkmers %>% mutate(gapLen= (lapply(kmer %>%  stringr::str_split("_"), function(x){x[2]}) %>% unlist))

    corr=cor(topkmers$ycounts,topkmers$xcounts,method = "pearson") %>% prettyNum(digit=2)
    kmerLen_t=kmerLen*(gapNo+1)
    topkmers %<>% mutate(gapLen=(gapLen %>% str_extract("[^n]*n$") %>% parse_number())+(gapLen %>% parse_number()),Ghomo=((topkmers$kmer %>% str_count("G"))>=(kmerLen_t-2))|((topkmers$kmer %>% str_count("C"))>=(kmerLen_t-2)))
    topkmers %<>% mutate(pcolor= ifelse(Ghomo,rgb(165,165,165,maxColorValue = 256),rgb(red = gapLen/max(gapLen),green = 0.2,blue = 1-gapLen/max(gapLen))))
    xyplot=ggplot(data=topkmers,aes(xcounts,ycounts,kmer=kmer)) + geom_point(aes(color=I(pcolor)),alpha=0.3)+ gg_theme_Publication() + gg_anno_grob(paste0("pearson corr: ", corr))

  }

  return(xyplot)
}

# p=xyplot("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-A1-BARHL1oooNIIIc4-CCCGCGCT_S1_L003_R1_001.peared_trimmed.fq.gz",
#          "~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4-A10-VSX1IIIc4-GCTCCGAT_S10_L003_R1_001.peared_trimmed.fq.gz",
#          test = F)
# print(p)

xyplot_RDSs<-function(filex,filey,xlab=NA,ylab=NA,alpha=0.2,size=0.2,quantile_cut=NA)
{
  if(is.na(xlab))xlab=rlang::eval_tidy(filex) %>% names
  if(is.na(ylab))ylab=rlang::eval_tidy(filey) %>% names
  topkmers=cbind( filex %>% readRDS() %>% rename(xcounts="counts"),filey %>% readRDS()%>% rename(ycounts="counts") %>% dplyr::select(ycounts))
  if(!is.na(quantile_cut))  dplyr::filter(xcounts+ycounts>quantile(xcounts+ycounts,quantile_cut))

  corr=cor(topkmers$ycounts,topkmers$xcounts,method = "pearson")%>% prettyNum(digit=2)
  topkmers %<>% mutate(cntA=str_count(kmer,"A"), cntT=str_count(kmer,"T"),cntC=str_count(kmer,"C"),cntG=str_count(kmer,"G"))
  maxCnt=topkmers$kmer[1] %>% nchar()
  topkmers %<>% mutate_at(c("cntA","cntT","cntC","cntG"),function(x) floor(x/maxCnt*256))
  topkmers %<>% mutate(Rval=(cntT+cntG) ,Gval=(cntA+cntG) ,Bval=(cntC) ,pcolor=rgb(red = Rval,green = Gval,blue = Bval,maxColorValue = 256))

  topkmers=topkmers[sample(1:nrow(topkmers)),]
  xyplot=ggplot() + geom_point(data=topkmers,aes(xcounts,ycounts,kmer=kmer,color=I(pcolor)),alpha=alpha,size=size)+ gg_theme_Publication() + gg_anno_grob(paste0("pearson corr: ", corr))+
    gg_anno_grob("A","darkgreen",0.95,0.5)+
    gg_anno_grob("C","blue",0.95,0.45)+
    gg_anno_grob("G","darkgoldenrod",0.95,0.4)+
    gg_anno_grob("T","red",0.95,0.35)+xlab(xlab)+ylab(ylab)
  return(xyplot)
}
