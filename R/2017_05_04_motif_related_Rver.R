## pure R version for motif mining

# motif from seed and seed from motif
            di_nt_to_IUPAC<-function(str_Arr)
              # c("CT","AT") %>% di_nt_to_IUPAC()
              # [1] "Y" "W"
            {
              str_Arr %>%  str_replace_all("AT|TA","W") %>%
                str_replace_all("GC|CG","S") %>%
                str_replace_all("AG|GA","R") %>%
                str_replace_all("CT|TC","Y") %>%
                str_replace_all("AC|CA","M") %>%
                str_replace_all("GT|TG","K")
            }

            IUPAC_to_di_nt_regex <- function(str,asArr=FALSE,target_pos_as_group=FALSE,gen_1_mismatch=TRUE)
            {
              # only interp to gen consensus
              if(!gen_1_mismatch) {return(str %>% str_replace_all("W","[AT]") %>%
                                            str_replace_all("S","[GC]") %>%
                                            str_replace_all("R","[AG]") %>%
                                            str_replace_all("Y","[CT]") %>%
                                            str_replace_all("M","[AC]") %>%
                                            str_replace_all("K","[GT]") %>%
                                            str_replace_all("N","."))}

              # gen all 1-mismatch seqs
              str %<>% str_replace_all("N",".")
              replace_pos=str %>% str_locate_all("[^\\.]") %>% .[[1]]  # get start & end pos of all individual char in string
              str%<>% `str_sub<-`(replace_pos[,1],replace_pos[,2],value="(.)")  %>%
                str_replace_all("W","[AT]") %>%
                str_replace_all("S","[GC]") %>%
                str_replace_all("R","[AG]") %>%
                str_replace_all("Y","[CT]") %>%
                str_replace_all("M","[AC]") %>%
                str_replace_all("K","[GT]")
              if(!asArr) str %<>% paste0(collapse="|")
              # if(target_pos_as_group) str %<>% str_replace_all("\\.","(.)")
              str
            }

            lamda_calc<-function(fg_reads,bg_reads)
            {
              median_8mer_frac <-function(reads)
              {
                median_8mer_frac= kmerCntBit(reads,8,collapse = T,pseudo = 0.5) %>% arrange(desc(counts))
                rows=nrow(median_8mer_frac); per25= ceiling(rows*0.25); per75= ceiling(rows*0.75)
                median_8mer_frac=  (median_8mer_frac[per25:per75,]$counts %>% sum())  /  (median_8mer_frac$counts %>% sum())
              }
              lamda= median_8mer_frac(fg_reads)/median_8mer_frac(bg_reads); lamda=ifelse(lamda>1,1,lamda); #attr(pfm_fg,'lamda (max 1)')=lamda;
              print("lamda = " %>% paste0(lamda))
              normalize_factor= length(fg_reads) / length(bg_reads);
              attr(lamda,"normalize_factor")=normalize_factor
              lamda
            }

            pfm_bk_corr_lamda <- function(PFM,PFM_bk,lamda)
              # lamda is with attr "normalize_factor" = length(fg_reads) / length(bg_reads)
            {
                normalize_factor=attr(lamda,"normalize_factor")
                PFM_bk_corr=PFM_bk * lamda *normalize_factor;
                PFM= PFM-PFM_bk_corr;
                PFM[PFM<0]=1
                PFM
            }

            pfm_correct_direction <- function(pfm,ic=FALSE,bk=0.3,return_logic=FALSE)
              # prefer more characters in A C G T order
            {
              pfm_topseq_sum<-function(pfm,ic=TRUE)
                {
                  values=apply(pfm,2,which.max) %>% c("A","C","G","T")[.] %>% str_replace("A","1") %>% str_replace("C","2") %>% str_replace("G","3") %>% str_replace("T","4") %>% as.integer()
                  if(ic) {icVal= apply(pfm+0.001,2,len4_vect_IC)+bk; values=values*icVal}
                  values %>% sum()
                }
              oriSum= pfm_topseq_sum(pfm, ic=ic)
              pfm_rc=fjComm::matRevComp(pfm)
              rcSum= pfm_topseq_sum(pfm_rc, ic=ic)

                if(return_logic){need_rc=if(oriSum<=rcSum) FALSE else TRUE; return(need_rc)}
              return_pfm= if(oriSum<=rcSum) pfm else pfm_rc %>% setattr("rc_used: ","yes")

              return_pfm %>% set_rownames(c("A","C","G","T"))
            }


pfm_from_seed_R<-function(fg_reads,seed,bg_reads=NULL,direc_corr=FALSE,rc=FALSE)
{
  if(is.na(seed)|(str_count(seed,"[^N]")<5)||is.null(seed)) print("pfm_from_seed_R: bad seed, return NA") & return(NA)
  if (rc) fg_reads=c(fg_reads,revComp(fg_reads))
  seed_expanded= seed %>% IUPAC_to_di_nt_regex(asArr = F,target_pos_as_group = T)
  allhits=str_match_all(fg_reads,seed_expanded) %>% {do.call(rbind,.)}

  # consensus only match to the 1st regex, but should match to all regex to count correctly, so treated separately
  cons= seed %>% IUPAC_to_di_nt_regex(gen_1_mismatch = F)
    cons_filter=allhits[,1] %>% str_detect(cons)
    consensus=allhits[cons_filter,][,1]
  cons_pfm=kmerCntBit(consensus,1L,F,pseudo = 1) %>% t

    mismatch=allhits[!cons_filter,]
  mis_pfm=apply(mismatch[,-1], 2, function(x){table(x) %>%  .[c("A", "C", "G", "T")] %>% as.integer()}) %>% {.[is.na(.)]=0; .}
  PFM=cons_pfm
  if( (seed %>% str_count("N"))>0)
  {
    Npos=seed %>% str_locate_all("N") %>% .[[1]] %>% .[,1]
    mis_pfm_Npos=lapply(Npos,function(x){mismatch[,1] %>% str_sub(x,x) %>% table() %>%  .[c("A", "C", "G", "T")] %>% as.integer()}) %>% as.data.frame() %>% as.matrix()
    PFM[,Npos]=PFM[,Npos]+mis_pfm_Npos; PFM[,-Npos]=PFM[,-Npos]+mis_pfm;
  }else{PFM=PFM+mis_pfm}

  if(!is.null(bg_reads))
  {
    PFM_bk=pfm_from_seed_R(bg_reads,seed)
    lamda=lamda_calc(fg_reads,bg_reads)
    PFM=pfm_bk_corr_lamda(PFM,PFM_bk,lamda)
  }
  if(direc_corr) PFM %<>% pfm_correct_direction()
  PFM
}

      # ##### wrong ver, fixed bases cannot be counted
      # pfm_from_seed_R<-function(fg_reads,seed,bg_reads=NULL,direc_corr=FALSE)
      # {
      #   seed_expanded= seed %>% IUPAC_to_di_nt_regex()
      #   allhits=str_extract(fg_reads,seed_expanded) %>% {.[!is.na(.)]}
      #   PFM=fjComm::kmerCntBit(allhits,1L,F) %>% t
      #     if(!is.null(bg_reads))
      #     {
      #       PFM_bk=pfm_from_seed_R(bg_reads,seed)
      #       lamda=lamda_calc(fg_reads,bg_reads)
      #       PFM=pfm_bk_corr_lamda(PFM,PFM_bk,lamda)
      #     }
      #     if(direc_corr) PFM %<>% pfm_correct_direction()
      #   PFM
      # }


seed_from_pfm_R <- function(pfm, weak_pos_cutoff=0.5, use_2nt_IUPAC_cutoff=c(0.5,2), include_flank_cutoff=2,get_consensus=FALSE)
  # 1. include flanking positions where the ratio between the most and least frequent bases was > 2 (include_flank_cutoff)
  # 2. for pos where max base freq < 0.5 (weak_pos_cutoff ). "N" subsitution,
  # 3. use IUPAC dint code if 2nd_nt_freq / 1st_nt_freq > 0.5 (use_2nt_IUPAC_cutoff, should be 0.75)      &      2nd_nt_freq / 3rd_nt_freq > 2 (use_2nt_IUPAC_cutoff)
{
  if(is.na(pfm)[1]) print("seed_from_pfm_R: pfm is NA, return NA") & return(NA)
  freq_pfm=pfm %>% sweep(2,colSums(.),FUN = "/")
  keepfilter=(colMaxs(freq_pfm)/colMins(freq_pfm))>=include_flank_cutoff
    if ((which(keepfilter) %>% length()) <=2)
    { print("Fun: seed_from_pfm_R     most positions has too low IC and discarded")
      if(get_consensus) { col_rank= apply(-freq_pfm,2,function(x)rank(x,ties.method = "random"))
                          newseed_from_pfm=c("A","C","G","T")[apply(col_rank,2,function(x)which(x==1))]
                          return(newseed_from_pfm %>% paste0(collapse = ""))} else { return(NA) }
     }
  keepfilter[which(keepfilter) %>% {min(.):max(.)}]=TRUE # only trim the flanking columns
  freq_pfm=freq_pfm[,keepfilter]
  weak_pos_filter=colMaxs(freq_pfm)<weak_pos_cutoff

  col_rank= apply(-freq_pfm,2,function(x)rank(x,ties.method = "random"))
  newseed_from_pfm=c("A","C","G","T")[apply(col_rank,2,function(x)which(x==1))]
  if(get_consensus) return(newseed_from_pfm %>% paste0(collapse = ""))
  # corr
    use2ntSym_filter1= (freq_pfm[col_rank==2]/freq_pfm[col_rank==3]>=use_2nt_IUPAC_cutoff[2])
    use2ntSym_filter2= (freq_pfm[col_rank==2]/freq_pfm[col_rank==1]>=use_2nt_IUPAC_cutoff[1])
    use2ntSym_filter=use2ntSym_filter1&use2ntSym_filter2
    # browser()
        # not corr
        # use2ntSym_filter= (freq_pfm[col_rank==2]/freq_pfm[col_rank==3]>=use_2nt_IUPAC_cutoff[2]) & freq_pfm[col_rank==2]>=0.06
  second_nt_letters=rownames(freq_pfm)[apply(col_rank,2,function(x)which(x==2))]
  top2nt_target_pos=paste0(newseed_from_pfm,second_nt_letters)[use2ntSym_filter] %>% di_nt_to_IUPAC()

  newseed_from_pfm[weak_pos_filter]="N"
  newseed_from_pfm[use2ntSym_filter]=top2nt_target_pos
  newseed_from_pfm %>% paste0(collapse="")
}

# seed_from_pfm_R <- function(pfm, weak_pos_cutoff=0.5, use_2nt_IUPAC_cutoff=2, include_flank_cutoff=2)
#   # 1. include flanking positions where the ratio between the most and least frequent bases was > 2 (include_flank_cutoff)
#   # 2. for pos where max base freq < 0.5 (weak_pos_cutoff ). "N" subsitution,
#   # 3. use IUPAC dint code if 2nd_nt_freq / 3rd_nt_freq > 2 (use_2nt_IUPAC_cutoff)
# {
#   freq_pfm=pfm %>% sweep(2,colSums(.),FUN = "/")
#   keepfilter=(colMaxs(freq_pfm)/colMins(freq_pfm))>=include_flank_cutoff
#   keepfilter[which(keepfilter) %>% {min(.):max(.)}]=TRUE # only trim the flanking columns
#   freq_pfm=freq_pfm[,keepfilter]
#   weak_pos_filter=colMaxs(freq_pfm)<weak_pos_cutoff
#
#   newseed_from_pfm= rownames(freq_pfm)[apply(freq_pfm,2,which.max)]
#   col_rank= apply(-freq_pfm,2,function(x)rank(x,ties.method = "random"))
#   use2ntSym_filter= (freq_pfm[col_rank==2]/freq_pfm[col_rank==3]>=use_2nt_IUPAC_cutoff) & weak_pos_filter
#   # use2ntSym_filter= (freq_pfm[col_rank==2]/freq_pfm[col_rank==1]>=use_2nt_IUPAC_cutoff) #& weak_pos_filter
#   second_nt_letters=rownames(freq_pfm)[apply(col_rank,2,function(x)which(x==2))]
#   top2nt_target_pos=paste0(newseed_from_pfm,second_nt_letters)[use2ntSym_filter] %>% di_nt_to_IUPAC()
#
#   newseed_from_pfm[weak_pos_filter]="N"
#   newseed_from_pfm[use2ntSym_filter]=top2nt_target_pos
#   newseed_from_pfm %>% paste0(collapse="")
# }


get_top_seed <- function(c4file, c0file, flankLen=5, half_site_Len=5,  single_string_seed=TRUE, oneOffname_for_bk=NULL, oneOffname_for_fg=NULL)
  # no revcomp
{
  if(!is.null(oneOffname_for_bk))
    bkcnt=oneOffCalc("oneoff/" %>% paste0(oneOffname_for_bk),calcFun = fjComm::gkmerCntBit_file,calcFunParamList = list(file = c0file,gapNo = 1,k = half_site_Len,gapMins = 0,gapMaxs = 10,pseudo = 10,melt_result = T,all_possible_k = T,rmdup = F),useScriptPath = T)
  else
    bkcnt=fjComm::gkmerCntBit_file(file = c0file,gapNo = 1,k = half_site_Len,gapMins = 0,gapMaxs = 10,pseudo = 10,melt_result = T,all_possible_k = T,rmdup = F)
  if(!is.null(oneOffname_for_fg))
    c4cnt=oneOffCalc("oneoff/" %>% paste0(oneOffname_for_fg),calcFun = fjComm::gkmerCntBit_file,calcFunParamList = list(file = c4file,gapNo = 1,k = half_site_Len,gapMins = 0,gapMaxs = 10,pseudo = 10,melt_result = T,all_possible_k = T),useScriptPath = T)
  else
    c4cnt=fjComm::gkmerCntBit_file(file = c4file,gapNo = 1,k = half_site_Len,gapMins = 0,gapMaxs = 10,pseudo = 10,melt_result = T,all_possible_k = T)

  c4cnt %<>% mutate(foldchn= (counts/sum(counts)) / (bkcnt$counts/ sum(bkcnt$counts))) %>% arrange(desc(foldchn)) %>% mutate(kmer=as.character(kmer))

  # use kmers with >=3 non homo nt
  topk_ind=1; cutoff=1.2
  while( (c4cnt$kmer[topk_ind] %>% str_split("") %>% .[[1]] %>% table() %>% len4_vect_IC())> cutoff ){topk_ind=topk_ind+1}

  seedLen=half_site_Len
  topkmer= c4cnt$kmer[topk_ind]
  seed1=str_sub(topkmer,1,seedLen); seed2=str_sub(topkmer,seedLen+1,-1)
  gapLen=c4cnt$gap[topk_ind] %>% str_sub(1,-2) %>% as.integer()
  flankLen=flankLen

  if(!single_string_seed){return(list(seed1=seed1,seed2=seed2,gapLen=gapLen,flankLen=flankLen))}else # for old ver
  {seed_used=paste0(strrep("N",flankLen),seed1,strrep("N",gapLen),seed2,strrep("N",flankLen)); return(seed_used)}
}





iterate_seed <-function(seed,fg_reads,bg_reads=NULL, weak_pos_cutoff=0.5, use_2nt_IUPAC_cutoff=c(0.5,2),include_flank_cutoff=2,maxIter=5,cnter=0,lamda=NULL,debugmode=F,add_N_flank=T,rc=FALSE, direct_corr=TRUE)
  # no rev comp
{

  if(is.na(seed)|(str_count(seed,"[^N]")<7)||is.null(seed)) print("iter_seed: bad seed, return NA") & return(list(pfm=NA,seed=seed,iterates=cnter,lamda=lamda))
  print(paste0("iteration ",cnter,", seed:"));
  print(seed);
  if (add_N_flank) {oriseed=seed; seed=paste0("NN",seed,"NN")}

  if(rc){fg_reads=c(fg_reads,revComp(fg_reads)); if(!is.null(bg_reads)){bg_reads=c(bg_reads,revComp(bg_reads))} }
  # # allow some flanking
  #   oriseed=seed
  #   seed %<>% str_replace_all("^N?(?=[^N])|(?<=[^N])N?$","NN")

  # calc lamda if not exist
  if(!is.null(bg_reads)) bkcorr=TRUE else bkcorr=FALSE

  PFM= pfm_from_seed_R(fg_reads,seed)

  # bk corr
  if(bkcorr)
  {
    if(is.null(lamda)) lamda=lamda_calc(fg_reads,bg_reads)
    PFM_bk= pfm_from_seed_R(bg_reads,seed)
    PFM= pfm_bk_corr_lamda(PFM,PFM_bk,lamda)
  }
  if(debugmode) {browser()}

  newSeed=seed_from_pfm_R(PFM,weak_pos_cutoff=weak_pos_cutoff, use_2nt_IUPAC_cutoff=use_2nt_IUPAC_cutoff, include_flank_cutoff=include_flank_cutoff)
  if (add_N_flank) seed=oriseed

  #//use old seed if signal is too weak, force newSeed==seed to end loop
  if(is.na(newSeed)|(str_count(newSeed,"[^N]")<7)) {newSeed=seed; print("too weak signal, use old seed")}

  if(newSeed==seed|cnter>=maxIter){
    # newSeed %<>% str_replace_all("^NN+|NN+$","NN");
    if(direct_corr)
    {
      PFM=PFM %>% pfm_correct_direction();
      if(PFM %>% pfm_correct_direction(return_logic = TRUE)) seed=revComp(seed)
    }

    if(cnter>=maxIter) print("No convergence after" %>% paste0(maxIter," cyc of Iteration")) else print('Pfm converged after ' %>% paste0(cnter," cycs of iter"))
    return(list(pfm=PFM,seed=seed,iterates=cnter,lamda=lamda))}
  else{iterate_seed(newSeed,fg_reads,bg_reads,weak_pos_cutoff=weak_pos_cutoff, use_2nt_IUPAC_cutoff=use_2nt_IUPAC_cutoff,include_flank_cutoff=include_flank_cutoff,maxIter=maxIter, cnter=cnter+1,lamda=lamda,debugmode = debugmode,add_N_flank=add_N_flank)}
}





