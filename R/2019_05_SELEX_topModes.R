top_n_models<-function(file,topNo=3, local_max_Num=100, kmer_l=4, rc=TRUE, maxIter = 4)
{
  # automatically pick the top n models from a given SELEX file
  # get localMax kmers >> get models (extended from seed) >> rm duplicated localMax kmers by matching to consensus of curr model
  # !!! (1-mis, and gap+-1 not yet removed)

  # local max kmers
  sigcnt=local_max_file(file, local_max_Num, kmer_l)

  # prep reads and lm kmers
  sg_reads=file %>% read_lines() %>%  rmdup() %>% .[[1]]
  if(rc) sg_reads %>% {c(.,revComp(.))}
  sig_lms=sigcnt %>% dplyr::filter(localMax) %>% arrange(desc(counts))

  # get models for top kmers (extend by iteration), once a model obtained, remove localMaxes that match consensus of the model (1-mis, and gap+-1 not yet removed)
  models=vector("list",topNo)
  for (i in 1:topNo)
  {
    sig_lm_kmers=sig_lms$kmer; sig_lm_kmers=paste0(str_sub(sig_lm_kmers,1,kmer_l),strrep("N",sig_lms$gap),str_sub(sig_lm_kmers,kmer_l+1,2*kmer_l))
    sig_lm_kmers_regex= sig_lm_kmers %>% str_replace_all("N",".")
          if(is.na(sig_lm_kmers[i])) break()
    mode1=iterate_seed(seed = sig_lm_kmers[i],sg_reads,maxIter = maxIter)
    mode1_cons=mode1$pfm %>% seed_from_pfm_R(get_consensus = T)

      mode1$consensus=mode1_cons
      mode1$lm_kmer=sig_lm_kmers[i]
      mode1$lm_kmer_id=sig_lms$id[i]
      models[[i]]=mode1

      kmers_to_filter=str_detect(mode1_cons,sig_lm_kmers_regex) %>% {.[1:i]=FALSE;.}
      kmers_to_filter1=str_detect(mode1_cons %>% revComp(),sig_lm_kmers_regex) %>% {.[1:i]=FALSE;.}
    kmers_to_filter= kmers_to_filter | kmers_to_filter1
    sig_lms %<>% dplyr::filter(!kmers_to_filter)
  }
  models
}




top_mode_cmp <- function(sigFiles, bkFiles, kmer_l=4, topNo=5, local_max_Num=100,rc=TRUE,maxIter=6, outFileName="top_mode_cmp.pdf",test=FALSE)
{
  # get top n seeds and models from sig and bk files
  # combine sig seeds and bk seeds together
  # for each seed, get sig model and bk model
  # calculate the difference between PWM models of sig and bk
    ## 1. freq of each base divided by colMax of respective model
    ## 2. normalized_freq1 - normalized_freq2
  # plot all models
  # !!! need to find a way to save plot, and figure out better return value


  if(test){topNo=2; local_max_Num=10; rc=FALSE; maxIter=1}
  if((outFileName %>% dirname())!=".") system("mkdir -p " %>% paste0(outFileName %>% dirname()))

  # bkFiles= bkFiles %>% {.[str_detect(.,paste0("-",well,"-"))]}

  if(length(bkFiles)>1){stop("bkFile has multiple matches: " %>% paste0(bkFiles,collapse = "######"))}

    sig_models=top_n_models(sigFiles,topNo=topNo, local_max_Num=local_max_Num, kmer_l=kmer_l, rc=rc, maxIter = maxIter)
    bk_models=top_n_models(bkFiles,topNo=topNo, local_max_Num=local_max_Num, kmer_l=kmer_l, rc=rc, maxIter = maxIter)

    sig_reads=sigFiles %>% read_lines() %>% rmdup() %>% .[[1]] %>% {c(.,revComp(.))}
    bk_reads=bkFiles %>% read_lines() %>% rmdup() %>% .[[1]] %>% {c(.,revComp(.))}

    allseeds=c(sig_models,bk_models) %>% map(~.$seed);
    sig_pfms=map(allseeds, ~pfm_from_seed_R(sig_reads,.))%>% set_names(allseeds %>% as.character() %>% paste0("_",1:length(allseeds)))
    bk_pfms=map(allseeds, ~pfm_from_seed_R(bk_reads,.))%>% set_names(allseeds %>% as.character() %>% paste0("_",1:length(allseeds)))
    sig_pfms %<>% map(function(x){if(is.na(x)[1]) return(matrix(rep(1,16),nrow = 4)); x })
    bk_pfms %<>% map(function(x){if(is.na(x)[1]) return(matrix(rep(1,16),nrow = 4)); x })

    tofreq<-function(pfm) sweep(pfm,2,colMaxs(pfm),"/")
    diff_pfm=map(1:(topNo*2),function(i){ (sig_pfms[[i]] %>% tofreq ) - (bk_pfms[[i]] %>% tofreq ) +0.0000001 })

    plot_pfmlist_alternative(sig_pfms, bk_pfms, diff_pfm, pfm_per_row=8, totalpfmNum=length(allseeds), outFileName)
    list(sig_models=sig_models,bk_models=bk_models,allseeds=allseeds,sig_pfms=sig_pfms,bk_pfms=bk_pfms,diff_pfm=diff_pfm)

}



              plot_pfmlist_alternative <- function(pfmlist1,pfmlist2,pfmlist3, pfm_per_row=8, totalpfmNum, outFileName)
              {
                ## to plot pfm lists alternatively according to pfm_per_row
                sig_pfms=pfmlist1; bk_pfms=pfmlist2; diff_pfm=pfmlist3

                rows=( totalpfmNum/pfm_per_row) %>% base::ceiling()
                sig_plot=vector("list",rows); bk_plot=vector("list",rows); diff_plot=vector("list",rows)

                fillNull<- function(x,val=0.0000001)
                {
                  NAnum=is.na(names(x)) %>% which() %>% length(); replacement=rnorm(NAnum)
                  names(x)[is.na(names(x)) %>% which()]=replacement %>% prettyNum(digit=8)
                  x %>% map(function(x1){ if(is.null(x1)) matrix(rep(val,16),nrow = 4) else x1  })
                }
                for(i in 1:rows)
                {
                  start= pfm_per_row*(i-1)+1;  end= pfm_per_row*i
                  sig_plot[[i]]= sig_pfms[start:end] %>% fillNull(1) %>% ggseqlogo_lab_list(ncol=pfm_per_row,method="probability")+ylab(cond)
                  bk_plot[[i]]=bk_pfms[start:end] %>% fillNull(1) %>% ggseqlogo_lab_list(ncol=pfm_per_row,method="probability")+ylab("bk_freq")
                  diff_plot[[i]]=diff_pfm[start:end] %>% fillNull(0.0000001) %>% set_names(NULL) %>% ggseqlogo_lab_list(method='custom',ncol=pfm_per_row)+ylab("diff (freq/colMax(freq))")
                }

                allplots=vector("list",3*rows)
                for(i in 1:rows)
                {
                  allplots[[(i-1)*3+1]]=sig_plot[[i]]
                  allplots[[(i-1)*3+2]]=bk_plot[[i]]
                  allplots[[(i-1)*3+3]]=diff_plot[[i]]
                }


                pdf(outFileName, width = 30,height = 8*rows)
                # gg_multiplot(sig_plot,bk_plot,diff_plot)
                gg_multiplot(plotlist = allplots)
                dev.off()
              }



top_mode_cmp_1seed <- function(sigFiles, bkFiles, seed="RTGGAAANWN",rc=TRUE,rmdup=TRUE)
{
  # get models from sig and bk files for seed
  # calculate the difference between PWM models of sig and bk
  ## 1. freq of each base divided by colMax of respective model
  ## 2. normalized_freq1 - normalized_freq2
  # plot all models

  sg_reads=sigFiles %>% read_lines()
  bk_reads=bkFiles %>% read_lines()
  if (rmdup) {sg_reads %<>% rmdup() %>% .[[1]]; bk_reads %<>% rmdup() %>% .[[1]]}
  if (rc) {sg_reads %<>% {c(.,revComp(.))}; bk_reads %<>% {c(.,revComp(.))} }

  maxIter=0
  sig_pfm=iterate_seed(seed,sg_reads,maxIter = maxIter,direct_corr = F)
  bk_pfm=iterate_seed(seed,bk_reads,maxIter = maxIter,direct_corr = F)
  tofreq<-function(pfm) sweep(pfm,2,colMaxs(pfm),"/")
  diff_pfm= (sig_pfm$pfm %>% tofreq ) - (bk_pfm$pfm %>% tofreq ) +0.0000001

  pfmlist=list(sig_pfm=sig_pfm$pfm, bk_pfm=bk_pfm$pfm, diff_pfm=diff_pfm) #%>% ggseqlogo_lab_list(ncol=1) %>%ggtitle(seed)
  plotlist= pfmlist %>% purrr::map(~ggseqlogo_lab(.,method = "probability")) %>% set_names(c("sig_plot","bk_plot","diff_plot"))
  plotlist[[1]]= plotlist[[1]] + ylab("sig_pfm")
  plotlist[[2]]= plotlist[[2]] + ylab("bg_pfm")
  plotlist[[3]]= ggseqlogo_lab(pfmlist$diff_pfm,method = "custom") +ylab("diff \n(freq/colMax(freq))")
  print(gg_multi_ggoutput(plotlist))
  return(list(pfmlist=pfmlist,plotlist=plotlist))
}



