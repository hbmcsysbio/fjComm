ggseqlogo_lab_list<-function(mat_list,...)
{
  # method=c("probability","bits")
  # scales = c("free_x","free_y")
  mat_list=tryCatch(mat_list %>% map(~.x@motif), error = function(e) return(mat_list))
  mat_list=tryCatch(mat_list %>% map(~.x@profileMatrix), error = function(e) return(mat_list))
  p=mat_list %>% purrr::map(~{set_rownames(as.matrix(.),c("A","C","G","T"))}) %>% ggseqlogo(...)
  p+gg_theme_Publication()+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
}


ggseqlogo_save_A4_pdf <- function(mat_list, filePrefix="motifplots_p", motifs_per_page=60,...)
{
  #name147=lapply(mat_list,function(x)attr(x,"name")) %>% unlist(); names(mat_list)=name147 %>% str_replace("_.*$","") %>% str_replace("ooo.*$","")
  mat_list %<>% purrr::map(~{set_rownames(as.matrix(.),c("A","C","G","T"))})
  length=length(mat_list); pages=(length/motifs_per_page) %>% as.integer()+1
  for (page in 1:pages)
  {
    start_=(motifs_per_page*(page-1)+1)
    end_=motifs_per_page*page; if(end_>length) end_=length
    # if(page==4)browser()
    rows=(((start_:end_) %>% length )/5) %>% base::ceiling()
    p147=mat_list[start_:end_] %>% ggseqlogo( ncol=5,...)+ scale_y_continuous(breaks = c(0,2)) +gg_theme_Publication()+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank() ) #theme(axis.text.x = element_blank(),title = )
    gg_save_pdf(p147,21,29.7/12*rows,filename = filePrefix %>% paste0(page))
  }
}

ggseqlogo_save_single_motifs_in_list <- function(named_mat_list,outDir="motif_vis/",width = 7,height = 3,outNamesVect=NA,titlesVect=NA,titleAdd="NCAP, ",directCorr=FALSE)
{
  # item_names= if(is.na(namesVect)) names(named_mat_list) else namesVect
  if (is.na(outNamesVect)) outNamesVect=names(named_mat_list)
  for(i in 1:length(named_mat_list))
  {
      if(directCorr)
      {   # corr direction
        named_mat_list[[i]]=named_mat_list[[i]] %>% as.matrix()
        named_mat_list_tmp=named_mat_list[[i]] %>% pfm_correct_direction(ic = T,bk=0.3)
        if (all_equal(named_mat_list_tmp,named_mat_list[[i]])[1]!=TRUE) {named_mat_list[[i]]=named_mat_list_tmp; titlesVect[i]=revComp(titlesVect[i])}
      }
    p= named_mat_list[[i]] %>% as.matrix() %>% set_rownames(qw("A C G T")) %>% ggseqlogo(method="prob")
    p= p+theme(axis.title = element_blank(),axis.text = element_blank(),plot.margin = margin(.1,.1,.1,.1, "cm"))
    if (!is.na(titlesVect[1])) p=p+ggtitle(paste0(titleAdd, titlesVect[i]))+theme(plot.title = element_text(size=9,hjust = 0.5))
    # browser()
    gg_save_pdf(p,width = width,height = height,path = outDir,filename = outNamesVect[i])
    saveRDS(p,outDir %>% paste0(outNamesVect[i],".Robj"))
  }
}
