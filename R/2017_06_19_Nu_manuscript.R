nusel_binders<-function()
{
  readRDS("~/Nut_zhuData/Analysis/Analysis2/1702_MI_diag_classify/combined_147_TFMI/different_binders_E-MI_penetration/TF_classification_as_binders.Robj")
}

# source("R/TF_related.R")
# TF curations
nusel_get_147 <- function()
{
  batch_4_5_guide= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig147_curate_v1.txt") #%>% dplyr::filter(use==1 & order==1)
  batch_4_4_guide= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig147_curate_v1.txt") #%>% dplyr::filter(use==1 & order==1)
  lig147= rbind(with(batch_4_5_guide,data.frame(pos=pos,TF=TF,batch="FJ45",family=family,class=class,tag=Tag,use=use,order=order)),with(batch_4_4_guide,data.frame(pos=pos,TF=TF,batch="FJ44",family=family,class=class,tag="Flag",use=use,order=order))) %>% mutate(TF=stringr::str_replace(TF,"ooo.*","")) %>%
    mutate(TF=stringr::str_replace(TF,"(NKX\\d*)\\.(\\d*)","\\1-\\2")) %>% mutate(pos=as.character(pos))
  lig147$stdFamily= TF_classify(lig147$TF)

  binders= readRDS("~/Nut_zhuData/Analysis/Analysis2/1702_MI_diag_classify/combined_147_TFMI/different_binders_E-MI_penetration/TF_classification_as_binders.Robj")
  lig147$binderClass_147= binders[match(lig147$TF,names(binders))]
  lig147 %<>% mutate(c4_file=NA,TFctrl_file=NA,c5Nu_file=NA,c5Free_file=NA,c4_EMI_file=NA,c4_EMI_TFctrl_file=NA,c5_EMI_Nu_file=NA,c5_EMI_Free_file=NA)
      EMIfile_FJ44= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/analysis/2d_IC_related_pub/data/*_Trulig147v1IIIFJ4-4-CAP2-*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
      TFctrl_EMIfile_FJ44= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/analysis/2d_IC_related_pub/data/*_Trulig147v1IIIFJ4-4-TFctrl-*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
      EMIfile_FJ45= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*_Trulig147v1IIIFJ4-5-CAP2r-147c4-*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
      TFctrl_EMIfile_FJ45= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*_Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4-*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
      Nu_EMIfile_FJ45= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*_Trulig147v1IIIFJ4-5-CAP2r-147c5Nu*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
      Free_EMIfile_FJ45= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*_Trulig147v1IIIFJ4-5-CAP2r-147c5Free*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
    c4_files_FJ44= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_147/*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    c4_files_FJ45= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c4*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    TFctrl_files_FJ44= file=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_TFctrl/*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    TFctrl_files_FJ45= file=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-TFctrl147c4*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    c5Nu_files_FJ45= file=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c5Nu*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    c5Free_files_FJ45= file=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-5-CAP2r-147c5Free*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
  lig147[lig147$batch=="FJ44",] %<>% mutate(c4_file= c4_files_FJ44[pos],TFctrl_file= TFctrl_files_FJ44[pos],c4_EMI_file=EMIfile_FJ44[pos],c4_EMI_TFctrl_file=TFctrl_EMIfile_FJ44[pos])
  lig147[lig147$batch=="FJ45",] %<>% mutate(c4_file= c4_files_FJ45[pos],TFctrl_file= TFctrl_files_FJ45[pos],c5Nu_file=c5Nu_files_FJ45[pos],c5Free_file= c5Free_files_FJ45[pos],c4_EMI_file=EMIfile_FJ45[pos],c4_EMI_TFctrl_file=TFctrl_EMIfile_FJ45[pos],
                                            c5_EMI_Nu_file=Nu_EMIfile_FJ45[pos],c5_EMI_Free_file=Free_EMIfile_FJ45[pos])

  return(lig147 %>% distinct() %>% mutate_at(vars(contains("file")),funs(str_replace(.,"/Users/zhu","~"))))
}

nusel_get_200 <- function()
{
  batch_4_5_guide= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.5CAP2r_lig200_curate_v1.txt") #%>% dplyr::filter(use==1 & order==1)
  batch_4_4_guide= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.4CAP2_lig200_curate_v1.txt") #%>% dplyr::filter(use==1 & order==1)
  batch_4_6_FL_guide= read_tsv("~/Nut_zhuData/SELEXphp/guide/FJ4.6CAPFL_lig200_Fangj201709.txt") %>% dplyr::filter(class=="2.3 C2H2 zinc finger") %>% dplyr::filter(!str_detect(TF,"Mecp"))
  lig200= rbind(with(batch_4_5_guide,data.frame(pos=pos,TF=TF,batch="FJ45",family=family,class=class,use=use,order=order)),with(batch_4_4_guide,data.frame(pos=pos,TF=TF,batch="FJ44",family=family,class=class,use=use,order=order)), with(batch_4_6_FL_guide,data.frame(pos=pos,TF=TF,batch="FJ46",family=family,class=class,use=use,order=order))) %>% mutate(TF=stringr::str_replace(TF,"ooo.*","")) %>%
    mutate(TF=stringr::str_replace(TF,"(NKX\\d*)\\.(\\d*)","\\1-\\2")) %>% mutate(pos=as.character(pos))
  lig200$stdFamily= TF_classify(lig200$TF)

  binders= readRDS("~/Nut_zhuData/Analysis/Analysis2/1702_MI_diag_classify/combined_147_TFMI/different_binders_E-MI_penetration/TF_classification_as_binders.Robj")
  lig200$binderClass_147= binders[match(lig200$TF,names(binders))]
  lig200 %<>% mutate(c4_file=NA,TFctrl_file=NA,c5Nu_file=NA,c5Free_file=NA,c4_EMI_file=NA,c4_EMI_TFctrl_file=NA)
    c4_files_FJ44= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_200/*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    c4_files_FJ45= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-200c4*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    TFctrl_files_FJ44= file=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/rmdup_TFctrl/*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    TFctrl_files_FJ45= file=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-TFctrl200c4*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    c5Nu_files_FJ45= file=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-200c5Nu*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    c5Free_files_FJ45= file=Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/allreads/Adpt_Trimmed_Reads/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-5-CAP2r-200c5Free*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    EMIfile_FJ44= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/analysis/2d_IC_related_pub/data/*_Trulig200v1IIIFJ4-4-CAP2-*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
    EMIfile_FJ45= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*_Trulig200v1IIIFJ4-5-CAP2r-200c4-*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
    TFctrl_EMIfile_FJ45= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.5CAP2r_PE/analysis/2d_IC_related_pub/data/*_Trulig200v1IIIFJ4-5-CAP2r-TFctrl200c4-*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
  lig200[lig200$batch=="FJ44",] %<>% mutate(c4_file= c4_files_FJ44[pos],TFctrl_file= TFctrl_files_FJ44[pos],c4_EMI_file=EMIfile_FJ44[pos])
  lig200[lig200$batch=="FJ45",] %<>% mutate(c4_file= c4_files_FJ45[pos],TFctrl_file= TFctrl_files_FJ45[pos],c5Nu_file=c5Nu_files_FJ45[pos],c5Free_file= c5Free_files_FJ45[pos],c4_EMI_file=EMIfile_FJ45[pos],c4_EMI_TFctrl_file=TFctrl_EMIfile_FJ45[pos])
    c4_files_FJ46= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.6CAPFL/allreads/Adpt_Trimmed_Reads_FJ4.6-CAPFL/2_Trim_adaptor_kmer/Trulig200v1IIIFJ4-6-CAPFL-200c4*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    TFctrl_files_FJ46= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.6CAPFL/allreads/Adpt_Trimmed_Reads_FJ4.6-CAPFL/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-6-CAPFL-147c4TFctrl*.gz") %>% set_names(lapply(., extract_well) %>% unlist)
    EMIfile_FJ46= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.6CAPFL/analysis/2d_IC_related_pub/data/*_Trulig200v1IIIFJ4-6-CAPFL-200c4-*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)
    # TFctrl_EMIfile_FJ46= Sys.glob("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/analysis/2d_IC_related_pub/data/*_Trulig147v1IIIFJ4-4-TFctrl-*_3mer_foldchn_maxBias.csv") %>% set_names(lapply(., extract_well) %>% unlist)

  lig200[lig200$batch=="FJ46",] %<>% mutate(c4_file= c4_files_FJ46[pos],TFctrl_file= TFctrl_files_FJ46[pos],c4_EMI_file=EMIfile_FJ46[pos])
  return(lig200 %>% distinct() %>% mutate_at(vars(contains("file")),funs(str_replace(.,"/Users/zhu","~"))))
}


# TF_std_name("BHLHB2")
TF_std_name <- function(TF_non_std_name)
{
  pacman::p_load(limma,org.Hs.eg.db)
  alias2Symbol(TF_non_std_name)
}

