
# TF="OLIG2"
# cond="30deg"
# sigFiles="~/Nut_zhuData/seqFiles2/FJ10tempSolVar/allreads/c3c4b2/" %>% paste0(cond,"/Trulig147*",TF,"*IIIc4*") %>% Sys.glob()
# bkFiles="~/Nut_zhuData/seqFiles2/FJ10tempSolVar/allreads/c3c4b2/0deg/Trulig147*" %>% paste0(TF,"*IIIc4*") %>% Sys.glob()
# topKmersEach= 3000
# kmer_l=4
# local_max_Num=100
# sigcnt=local_max_file(sigFiles, local_max_Num, kmer_l)
# bkcnt=local_max_file(bkFiles, local_max_Num, kmer_l)
# xyplot=xyplot_gklm(sigcnt, bkcnt, topKmersEach)
# plotly::ggplotly(xyplot)


hudd_kmer_gen <-function(kmer, gap=0, rc=TRUE)
{
  ## gen all kmers with hudd dist <1 (incl. rc kmer with hudd_dist 0), or with frame shift 1, incl hudd1 kmers for rc
  kmerlen=nchar(kmer)
  kmers_gapv=kmer
  gaps_gapv=if(gap==0) 1 else c(gap-1,gap+1)

  kmers_kmerv= {
    c(kmer %>% stringi::stri_sub_replace(1:kmerlen,1:kmerlen,value = "A"),
      kmer %>% stringi::stri_sub_replace(1:kmerlen,1:kmerlen,value = "C"),
      kmer %>% stringi::stri_sub_replace(1:kmerlen,1:kmerlen,value = "G"),
      kmer %>% stringi::stri_sub_replace(1:kmerlen,1:kmerlen,value = "T"))}
  gaps_kmerv=gap

  ## gen frame shift kmers
  halflen=kmerlen/2
  kmerh=kmer %>% str_sub(1,kmerlen/2)
  kmert=kmer %>% str_sub(kmerlen/2+1,kmerlen)
  ACGT=qw("A C G T")

  kmerh_shift_l= paste0(ACGT, kmerh %>% str_sub(1,halflen-1))
  kmert_shift_l= paste0(ACGT, kmert %>% str_sub(1,halflen-1))
  kmerh_shift_r= paste0(kmerh %>% str_sub(2,halflen),ACGT)
  kmert_shift_r= paste0(kmert %>% str_sub(2,halflen),ACGT)

  kmers_shift_whole= c(paste0(outer(kmerh_shift_l,kmert_shift_l,stringi::stri_c),gap),  paste0(outer(kmerh_shift_r,kmert_shift_r,stringi::stri_c),gap),
                       paste0(outer(kmerh_shift_r,kmert_shift_l,stringi::stri_c),gap-2),  paste0(outer(kmerh_shift_l,kmert_shift_r,stringi::stri_c),gap+2))
  kmers_shift_h= c(paste0(outer(kmerh_shift_l,kmert,stringi::stri_c),gap+1), paste0(outer(kmerh_shift_r,kmert,stringi::stri_c),gap-1))
  kmers_shift_t= c(paste0(outer(kmerh,kmert_shift_r,stringi::stri_c),gap+1), paste0(outer(kmerh,kmert_shift_l,stringi::stri_c),gap-1))

  hudd_kmers=c(paste0(kmers_gapv,gaps_gapv), paste0(kmers_kmerv,gaps_kmerv), kmers_shift_whole, kmers_shift_h, kmers_shift_t)
  if(rc & kmer!=revComp(kmer)) hudd_kmers=c(hudd_kmers,hudd_kmer_gen(kmer = revComp(kmer),gap = gap,rc = FALSE),
                                            paste0(revComp(kmer),gap))
  hudd_kmers %<>% setdiff(paste0(kmer,gap))

  hudd_kmers
}


local_max_gkdf<- function(gkcnt, local_max_Num=100, col=3)
{
  # pass a melted gkdf, annotate the localMax kmers
  # id(pasted kmer and gap e.g. AAAAACCC7)
  gkcnt %<>% mutate_if(is.factor,as.character) %>% mutate(gap=parse_number(gap) %>% as.integer()) %>% mutate(id=paste0(kmer,gap))
  gkbk=gkcnt

  # not considering G/C polymer
    kmerLen_t= gkcnt$kmer[[1]] %>% nchar
  gkcnt %<>% arrange(desc(.[[3]])) %>% dplyr::filter( !(((kmer %>% str_count("G"))>=(kmerLen_t-2))|((kmer %>% str_count("C"))>=(kmerLen_t-2))) )
  for (i in 1:local_max_Num)
  {

    kmer=gkcnt$kmer[i]
    # if(gkcnt$id[i]=="AACCATAT0")
    # View(gkcnt); browser()
    #ATATGGTT sig
    gap=gkcnt$gap[i]
    id_filter= hudd_kmer_gen(kmer,gap,TRUE)
    gkcnt %<>% dplyr::filter(!(id %in% id_filter))
    if(gkcnt$kmer[i]!=kmer) stop("removed current kmer!!")
  }
  lm_kmers=gkcnt$id %>% head(local_max_Num)
  gkbk %>% mutate(localMax=ifelse(id %in% c(lm_kmers),TRUE, FALSE))
}

local_max_file<- function(file, local_max_Num=100, kmerLen=4, rc_combine=T, rmdup=T)
{
  # pass a filename, gen localMax kmers
  sigFiles=file
  sigcnt= gkmerCntBit_file(sigFiles,1,kmerLen,gapMins = 0,gapMaxs = 10,rmdup = rmdup,rc_combine =rc_combine ,rc_uniq = F)
  sigcnt %>% local_max_gkdf(local_max_Num,3)
}

xyplot_gklm<- function(sig_lm_df,bk_lm_df,topKmersEach= 3000)
{
  ## xyplot indicating local max, triangle: both, filled: sig, non-filled: bk
  half_k_len=nchar(sig_lm_df$kmer[[1]])/2
  sig_lm_df %<>% mutate(localMax1=ifelse(bk_lm_df$localMax,"bk", "none"),localMax1=ifelse(localMax,"sig", localMax1),localMax1=ifelse(localMax & bk_lm_df$localMax,"both", localMax1)) %>% mutate(kmer_=paste0(str_sub(kmer,1,half_k_len),paste0(gap,"n"),str_sub(kmer,half_k_len+1,half_k_len*2)))
  sig_lm_df %<>% select(-localMax) %>% rename(localMax="localMax1",cnt_sig="counts") %>% mutate(cnt_bk=bk_lm_df$counts)

  gkcnt_sig=sig_lm_df %>% top_n(topKmersEach,cnt_sig/cnt_bk)
  gkcnt_bk=sig_lm_df %>% top_n(topKmersEach,cnt_bk/cnt_sig)
  gkcnt_sigc=sig_lm_df %>% top_n(topKmersEach,cnt_sig)
  gkcnt_bkc=sig_lm_df %>% top_n(topKmersEach,cnt_bk)
  gkcnt_bksigc=sig_lm_df %>% top_n(topKmersEach,cnt_bk+cnt_sig)
  topkmers=dplyr::union(gkcnt_sig,gkcnt_bk) %>% dplyr::union(gkcnt_sigc) %>% dplyr::union(gkcnt_bkc) %>% dplyr::union(gkcnt_bksigc)


  corr=cor(topkmers$cnt_sig,topkmers$cnt_bk,method = "pearson") %>% prettyNum(digit=2)
  kmerLen_t=nchar(topkmers$kmer[[1]])
  topkmers %<>% mutate(Ghomo=((topkmers$kmer %>% str_count("G"))>=(kmerLen_t-2))|((topkmers$kmer %>% str_count("C"))>=(kmerLen_t-2)) )
  topkmers %<>% mutate(pcolor= ifelse(Ghomo,rgb(165,165,165,maxColorValue = 256),rgb(red = gap/max(gap),green = 0.2,blue = 1-gap/max(gap))))
  # topkmers %<>% mutate(shape_=case_when(localMax == "sig" ~ 115, localMax == "bk" ~ 98,TRUE ~ 1), size_=ifelse(localMax=="none",0.2,1.5))
  xyplot=
    ggplot()+ scale_shape_identity()+
    geom_point(aes(cnt_sig,cnt_bk,kmer=kmer_,color=I(pcolor)),data=topkmers %>% dplyr::filter(localMax=="none"),size=0.2,alpha=0.3)+
    geom_point(aes(cnt_sig,cnt_bk,kmer=kmer_,color=I(pcolor),lm="bk"),data=topkmers %>% dplyr::filter(localMax=="bk"),size=1.8,alpha=0.2,shape=1)+
    geom_point(aes(cnt_sig,cnt_bk,kmer=kmer_,color=I(pcolor),lm="sig"),data=topkmers %>% dplyr::filter(localMax=="sig"),size=1.8,alpha=0.2)+
    geom_point(aes(cnt_sig,cnt_bk,kmer=kmer_,color=I(pcolor),lm="both"),data=topkmers %>% dplyr::filter(localMax=="both"),size=2,alpha=0.2,shape=17)+
    gg_theme_Publication() + gg_anno_grob(paste0("pearson corr: ", corr))

  return(xyplot)
}
