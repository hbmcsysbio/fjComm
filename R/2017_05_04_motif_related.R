allGuidesResource= c("~/Nut_zhuData/SELEXphp/guide/all-human-version4.txt" ,
             "~/Nut_zhuData/SELEXphp/guide/Human_IndividualTF_Cell_eLife_Nature_YimengsMethPaper.txt",
             "~/Nut_zhuData/SELEXphp/guide/MethylSELEX_v0.7_FinalCall_without_fail.txt",
             # "~/Nut_zhuData/SELEXphp/guide/MethylSELEX_v0.7_FinalCall.txt",
             "~/Nut_zhuData/SELEXphp/guide/Methylation_SELEX_Version_0.5.txt"
)

# # matrix reverse complement //C version used in seed_to_pfm.cpp
# matRevComp <- function(x) "[<-"(x, rev(x))

motif_getpfm_from_TFname <- function(TFname, guidefiles=allGuidesResource, monomericFirst=TRUE, exactMatchFirst=TRUE, head=20)
{
  seqobj=seqFile("")
  seqobj$getPfmFiles(guidefiles = guidefiles,TFname_full = TFname,head = head,monomericFirst = monomericFirst, exactMatchFirst = exactMatchFirst)
  return(seqobj$pfmFiles)
}


# match reads with pfm file, and plot
motif_plot_moods_singlePWM_hit <- function(readsfile, pfmfile, dup_rm=TRUE, p_moods="0.0001", colorGrad= c("#ceffce","red"), strand_filter=c("+","-"), combine_2strand_hits=FALSE)
{
  selex= SELEXFile(readsfile)
  selex$getSeq(dup_rm = dup_rm)
  selex$moodsMap(pfmfile, p=p_moods )
  selex$moodsResult= selex$moodsResult %>% group_by(pos,strand) %>% summarise(count=n()) %>% ungroup() # bin counts with the same strand and pos

  if (length(strand_filter)==2 & (!combine_2strand_hits))
  {
    selex$moodsResult= selex$moodsResult %>% group_by(strand) %>% mutate(ncount= prop.table(count)) %>% ungroup()
    plot=ggplot(selex$moodsResult)+ geom_tile(aes(pos+1,strand,fill=ncount))+ scale_fill_gradientn(colors =colorGrad)+ guides(fill=F)  +ylab("Strand")+
      gg_axis_x_noExp(breaks=gg_breaks(selex$moodsResult$pos+1))+ scale_y_discrete(expand = c(0,0)) + gg_theme_Publication() + theme(axis.line= element_blank(),axis.ticks.y = element_blank())
  }else{
    if (combine_2strand_hits) strand_filter=c("+","-")
    if (length(strand_filter)==1) {selex$moodsResult= selex$moodsResult %>% dplyr::filter(strand %in% strand_filter) %>% mutate(ncount= prop.table(count))}
    else if (length(strand_filter)==2){selex$moodsResult= selex$moodsResult %>% group_by(pos) %>% summarise(count=sum(count)) %>% ungroup() %>% mutate(ncount= prop.table(count)) }
    else {stop("strand filter length is not correct !!")}
    plot=ggplot(selex$moodsResult)+ geom_tile(aes(pos+1,1,fill=ncount))+ scale_fill_gradientn(colors =colorGrad)+ guides(fill=F) +
      gg_axis_x_noExp(breaks=gg_breaks(selex$moodsResult$pos+1,end_plus_middle = T))+gg_axis_y_noExp()+ gg_theme_Publication_diag
  }
  return(plot+xlab("(bp)")+scale_x_continuous(expand = c(0, 0),breaks =  gg_breaks(selex$moodsResult$pos+1,end_plus_middle = T)))
}

# # old density ver
# motif_plot_moods_singlePWM_hit <- function(readsfile, pfmfile, dup_rm=TRUE, p_moods="0.0001", colorGrad= c("#ceffce","red"), strand_filter=c("+","-"), combine_2strand_hits=FALSE)
# {
#   selex= SELEXFile(readsfile)
#   selex$getSeq(dup_rm = dup_rm)
#   selex$moodsMap(pfmfile, p=p_moods )
#   if (length(strand_filter)==2 & (!combine_2strand_hits))
#   {
#     plot=ggplot(selex$moodsResult)+ geom_tile(aes(pos+1,strand,fill=..density..),stat = "density")+ scale_fill_gradientn(colors =colorGrad)+ guides(fill=F)  +
#       gg_axis_x_noExp()+ scale_y_discrete(expand = c(0,0)) + gg_theme_Publication() + theme(axis.line= element_blank(),axis.ticks.y = element_blank())
#   }else{
#     if (combine_2strand_hits) strand_filter=c("+","-")
#     plot=ggplot(selex$moodsResult %>% dplyr::filter(strand %in% strand_filter))+ geom_tile(aes(pos+1,1,fill=..density..),stat = "density")+ scale_fill_gradientn(colors =colorGrad)+ guides(fill=F) +
#       gg_axis_x_noExp()+gg_axis_y_noExp()+ gg_theme_Publication_diag
#   }
#   return(plot+xlab("position (bp)")+scale_x_continuous(expand = c(0, 0),breaks =  gg_axis_breaks(selex$moodsResult$pos,interval = 20)))
# }

### motifStack related ####
    plotMotif_pfmFile <-function(pfm_file,ic.scale=FALSE,title=NULL,full_path_on_fig=FALSE,TFname=NULL)
    {
      suppressPackageStartupMessages(library(motifStack))
      # pcm <- read.table(file.path(find.package("motifStack"),
      #                             "extdata", "bin_SOLEXA.pcm"))
      # pcm <- pcm[,3:ncol(pcm)]
      file=pfm_file
      pcm <- read.table(file)
      # pcm <- read.table("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/analysis/T_analysis/T_motif/CAC76NGTGmul_n.txt")
      rownames(pcm) <- c("A","C","G","T")
      motif <- new("pcm", mat=as.matrix(pcm), name=
                     ifelse(is.null(title), ifelse(full_path_on_fig, paste0(dirname(file)%>% str_replace("^/Users/zhu","~"),"/\n",basename(file),"\n",TFname ), basename(file) ), title)
                   )
      ##pfm object
      #motif <- pcm2pfm(pcm)
      #motif <- new("pfm", mat=motif, name="bin_SOLEXA")
    # opar<-par(mfrow=c(1,1))
      #plot the logo with same height
      plot(motif,ic.scale=ic.scale, ylab="IC")
      # #try a different font
      # plot(motif, font="mono,Courier")
      # #try a different font and a different color group
      # motif@color <- colorset(colorScheme='basepairing')
      # plot(motif,font="Times")
    }

    plotMotif_pfmMat <- function(pfm_mat,name="default",ic.scale=TRUE,ylab="IC",yaxis=T,xaxis=T)
    {
      suppressPackageStartupMessages(library(motifStack))
      rownames(pfm_mat) <- c("A","C","G","T")
      motif <- new("pcm", mat=as.matrix(pfm_mat), name=name)
      plot(motif, ylab=ylab,xlab="",ic.scale=ic.scale,yaxis=yaxis,xaxis=xaxis)
    }

    ggseqlogo_lab <- function(pfm,method=c("probability","bits"),...)
    {
      pacman::p_load(ggseqlogo)
      if (length(method)==2) method="bits"
      ggseqlogo(pfm %>% as.matrix() %>% set_rownames(c("A","C","G","T")),method=method,...)+theme(axis.text.x = element_blank())}

    # get motif obj from pfm_files
    getMotif <- function(pfm_file)
    {
      suppressPackageStartupMessages(library(motifStack))
      pcm <- read.table(pfm_file)
      rownames(pcm) <- c("A","C","G","T")
      motif= new("pcm", mat=as.matrix(pcm), name=basename(pfm_file))
      return(motif)
    }

    plotMotif_all_tops_by_TFname <- function(TFname, guidefiles=allGuidesResource, monomericFirst=TRUE, exactMatchFirst=TRUE, head=20, pdf_path_to_save=NULL)
    {
      top_motif_files=motif_getpfm_from_TFname(TFname, guidefiles, monomericFirst, exactMatchFirst, head)
        tempfileName=tempfile(pattern = "file", tmpdir = tempdir()) %>% paste0(".pdf")
        filename=ifelse(is.null(pdf_path_to_save),tempfileName,pdf_path_to_save)
      pdf(file=filename,height = 3)
      for (i in seq_len(nrow(top_motif_files))) {plotMotif_pfmFile(top_motif_files[i,2],ic.scale = T,full_path_on_fig = T,TFname =top_motif_files[i,1] )}
      # lapply(top_motif_files,function(x)(plotMotif_pfmFile(x,ic.scale = T,full_path_on_fig = T)))
      dev.off()
      if(is.null(pdf_path_to_save)) {system(paste0('open "', filename, '"')); system(paste0("rm ",filename))}
    }


### motif from raw reads using seed ######
      # element fun for 1 file
      pfm_from_seed_element <- function(seqs, seed1="AAA", gapLen=1, two_strands=TRUE, seed2="", flankLen=4,     all_start_with_specified_gap=T, seed1_start= c(2,3,4),     rmdup=TRUE,  ns_estimate=FALSE)
      {
        # // pfm=pfm_from_seed(seqs = seq,seed1 = "CAC", gapLen = 77, seed2 = "GTG", seed1_start = c(1,2,5),flankLen = 4,all_start_with_specified_gap=T,two_strands = T)
        # // pfm=pfm_from_seed(seqs = seq,seed1 = "TAA", gapLen = 0, seed2 = "CACCT",flankLen = 0, all_start_with_specified_gap=T,two_strands = T)


            from_file=FALSE
          if(length(seqs)==1 && typeof(seqs)=="character") {from_file=TRUE; print (paste0("treated as file: ",basename(seqs)))}
          if(from_file) {seqs=read_csv(seqs,col_names = FALSE); if(rmdup) seqs= rmdup(seqs); seqs=seqs[[1]]}

          if( (!from_file) && typeof(seqs)!="list" && rmdup) {seqs= rmdup(seqs)[[1]]}
          if((!from_file) && typeof(seqs)=="list") {  if(rmdup) seqs= rmdup(seqs)[[1]]} # so data frame also work

        readsNum=length(seqs)
        pfm = pfm_from_seed_notUsed(seqs = seqs, seed1= seed1, gapLen = gapLen, seed2 = seed2, seed1_start=seed1_start, flankLen=flankLen, all_start_with_specified_gap=all_start_with_specified_gap)
        if(two_strands)
        {
          pfm_rc=pfm_from_seed_notUsed(seqs = seqs %>% revComp(), seed1= seed1, gapLen = gapLen, seed2 = seed2, seed1_start=seed1_start, flankLen=flankLen, all_start_with_specified_gap=all_start_with_specified_gap)
          pfm= pfm+pfm_rc
          # pfm_rc=pfm_from_seed_notUsed(seqs = seqs, seed1= seed2 %>% revComp(), gapLen = gapLen, seed2 = seed1 %>% revComp(), seed1_start=seed1_start, flankLen=flankLen, all_start_with_specified_gap=all_start_with_specified_gap)
          # pfm= pfm+(pfm_rc %>% matRevComp())
        }

        if(ns_estimate) # estimate ns carryover as in genome research
        {
          median_8mer_frac= kmerCntBit(seqs,8,collapse = T,pseudo = 0.5) %>% arrange(desc(counts))
          rows=nrow(median_8mer_frac); per25= ceiling(rows*0.25); per75= ceiling(rows*0.75)
          median_8mer_frac=  (median_8mer_frac[per25:per75,]$counts %>% sum())  /  (median_8mer_frac$counts %>% sum())
          attr(pfm, "median_8mer_frac")=median_8mer_frac
        }
        attr(pfm,"readsNum")=readsNum
        return(pfm)
      }

# figure out motif from reads and seed
# seed_to_pfm.cpp
pfm_from_seed <- function(seqs_or_file, seqs_or_file_bg=NA, seed1="AAA", gapLen=1, two_strands=TRUE, seed2="", flankLen=4,     all_start_with_specified_gap=T, seed1_start= c(2,3,4),     rmdup_fg=TRUE, rmdup_bg=TRUE)
{
    # 1. fish out bg and fg PWM first
    # 2. estimate ns carryover using 25-75% 8mer fraction of k and k+1 cycle
    # 3. calc lamda, if >0.8 then =0.8
    # 4. substract pfm_bg from pfm_fg
    ## example ##
    # pfm_from_seed(seqs_or_file = after1file, seqs_or_file_bg = beforefile,seed1 = after1_TAAT$kmer[[1]],gapLen = 0,flankLen = 4, all_start_with_specified_gap = T,two_strands = T) %>% fjComm::plotMotif_pfmMat()

    if(seed2==""){gapLen=0; seed2=str_sub(seed1,-1,-1); seed1=str_sub(seed1,1,-2); print(seed1); print(seed2)}
    ns_estimate=FALSE; if( length(seqs_or_file_bg)>1 || (!is.na(seqs_or_file_bg))) ns_estimate=TRUE

  pfm_fg=pfm_from_seed_element(seqs=seqs_or_file, seed1=seed1, gapLen=gapLen, two_strands=two_strands, seed2=seed2, flankLen=flankLen,     all_start_with_specified_gap=all_start_with_specified_gap, seed1_start= seed1_start,     rmdup=rmdup_fg,ns_estimate=ns_estimate)

  # corr for bg if available
  if( length(seqs_or_file_bg)>1 || (!is.na(seqs_or_file_bg)))
  {
    pfm_bg=pfm_from_seed_element(seqs=seqs_or_file_bg, seed1=seed1, gapLen=gapLen, two_strands=two_strands, seed2=seed2, flankLen=flankLen,     all_start_with_specified_gap=all_start_with_specified_gap, seed1_start= seed1_start,     rmdup=rmdup_bg, ns_estimate=ns_estimate)
    lamda= attr(pfm_fg,"median_8mer_frac") / attr(pfm_bg,"median_8mer_frac"); lamda=ifelse(lamda>1,1,lamda); attr(pfm_fg,'lamda (max 1)')=lamda;
    normalize_factor= attr(pfm_fg,"readsNum") / attr(pfm_bg,"readsNum"); #normalize_factor= colSums(pfm_fg) / colSums(pfm_bg)
    pfm_bg_corr=pfm_bg * lamda *normalize_factor  ;#pfm_bg_corr=sweep(pfm_bg,2,normalize_factor,FUN="*") * lamda  ;
      pfm_fg= pfm_fg-pfm_bg_corr;
      pfm_fg[pfm_fg<0]=1
  }

  return(pfm_fg)
}




#### scoring kmerVect with pfm or pwm matrix ##########

  kmer_pfm_score <- function(kmerVect,pfmMat,normalize=TRUE)
    # calc freq for each kmer acc to pfm
    # if normalized each pfm col sum to 1
  {
    kmerLen= nchar(kmerVect[1])
    if (kmerLen!=ncol(pfmMat)) stop("length different for kmer and pfm")
    rownames(pfmMat)=c("A","C","G","T")
    if(normalize) pfmMat= sweep(pfmMat,2,colSums(pfmMat),"/")
    kmerFrag= fjComm::seqFregments(kmerVect,1L)
    prob_k= rep(1,length(kmerVect))
    for (i in 1:kmerLen)
    {
      prob_k=prob_k*pfmMat[kmerFrag[,i],i]
    }
    return(prob_k)
  }

  kmer_pwm_score <- function(kmerVect,pwmMat,normalize=TRUE)
    # calc energy for each kmer acc to pfm
    # if normalized each pwm col mean to 0
  {
    kmerLen= nchar(kmerVect[1])
    if (kmerLen!=ncol(pwmMat)) stop("length different for kmer and pfm")
    rownames(pwmMat)=c("A","C","G","T")
    if(normalize) pwmMat %<>% sweep(2,colMeans(.),FUN = "-")
    kmerFrag= fjComm::seqFregments(kmerVect,1L)
    energy_k= rep(0,length(kmerVect))
    for (i in 1:kmerLen)
    {
      energy_k=energy_k+pwmMat[kmerFrag[,i],i]
    }
    return(energy_k)
  }

#### inter-conver of pfm and pwm ##########

  pfm_to_pwm <-function(pfmMat, normalize=TRUE, pseudoCnt=0)
    # pwm energy in kT unit
    # if normalized each pwm col mean to 0
  {
    pwmMat= (pfmMat+pseudoCnt) %>% log() %>% `-`
    if (normalize) pwmMat %<>% sweep(2,colMeans(.),"-")
    pwmMat
  }

  pwm_to_pfm <-function(pwmMat, normalize=TRUE)
    # pwm energy in kT unit
    # if normalized each pfm col sum to 1
  {
    pfmMat= (-pwmMat) %>% base::exp()
    if (normalize) pfmMat %<>% sweep(2,colSums(.),"/")
    pfmMat
  }

  pfm_normalize <-function(pfmMat) pfmMat %>% sweep(2,colSums(.),"/")
  pwm_normalize <-function(pwmMat) pwmMat %>% sweep(2,colMeans(.),"-")


### EM-related ###

              nextPfm_calc <- function(kmerV, kcountV, pfmScoreV)
              {
                kmerLen=nchar(kmerV[1])
                kmerFrag= seqFregments(kmerV,1L) %>% as.data.frame() %>% mutate(pfmScore=pfmScoreV, kcounts=kcountV)
                cnt_all_pos= vector("list",kmerLen)
                for (i in 1:kmerLen) cnt_all_pos[[i]]=kmerFrag %>% group_by(.[[i]]) %>% dplyr::count(wt = pfmScore*kcounts)
                nextPfm=lapply(cnt_all_pos,function(tt){tt[[2]][match(tt[[1]],c("A","C","G","T"))]}) %>% {do.call(cbind,.)} %>% set_rownames(c("A","C","G","T"))
                nextPfm %>% pfm_normalize()
              }

      EM_pfm <- function(kmerV, kcountV, pfm=NULL, rounds=1L, bkcountV=NULL)
        # optimize current pfm by EM
        # if no pfm use random one
        # divid bk counts when provided
      {
        if (is.null(pfm)) pfm=runif(4*nchar(kmerV[1])) %>% matrix(nrow=4)
        currpfm=pfm
        for (round in 1:rounds)
        {
          pfmScoreV= kmer_pfm_score(kmerV,currpfm)
          currpfm= nextPfm_calc(kmerV,kcountV,pfmScoreV)
          # divide bk counts
          if(!is.null(bkcountV))
          {
            currBkpfm= nextPfm_calc(kmerV,bkcountV,pfmScoreV)
            currpfm= currpfm/currBkpfm
            # currpfm= ((currpfm %>% pfm_to_pwm())- (currBkpfm %>% pfm_to_pwm())) %>% pwm_to_pfm()
          }
        }
        currpfm %>% pfm_normalize()
      }


  EM_motif_from_seq <- function(seq_charV, motifs=5L, motifLen=10, EM_rounds=10L, is.file=FALSE)
    # de novo motif discovery
    # no rmdup yet
    # redundant segmentation to be optimized
  {
    if(is.file) {kcnt=fjComm::kmerCntBit(seq_charV %>% read_csv(col_names = FALSE) %>% .[[1]],k = motifLen,collapse = T)}
      else {kcnt=fjComm::kmerCntBit(seq_charV,k = motifLen,collapse = T)}
    result_motifs=vector("list",motifs)
    for (i in 1:motifs)
    {
      pfm=runif(4*motifLen) %>% matrix(nrow=4)
      pfm=EM_pfm(kcnt$kmer,kcnt$counts,pfm, rounds = EM_rounds)
      result_motifs[[i]]=pfm
    }
    result_motifs
  }



  moodsMap<-function(reads,pfm,batch=T,pval=0.001,returnRaw=TRUE,bindingPerLig=FALSE)
  {
    # returnRaw: true to return raw result, false to return counts by position
    # if hope to distinguish monomer form dimer, set bindingPerLig=TRUE
    selex=SELEXFile("tt");
    selex$seq=reads %>% Biostrings::DNAStringSet() %>% set_names(1:length(reads))
    selex$moodsMap(pwmFile =pfm, batch = batch, p = pval )

    moods=selex$moodsResult
    if(returnRaw) return(moods)
    if(bindingPerLig)
    {
      p_load(Biobase)
      moods %>%mutate(monomer_lig=isUnique(rangeNo)) %>%  group_by(strand,monomer_lig) %>% nest() %>% mutate(data=map(data,~table(.$pos) %>% as.data.frame() %>% set_colnames(c("pos","cnt"))  )) %>% unnest()
    }else{
      moods %>% group_by(strand) %>% nest() %>% mutate(data=map(data,~table(.$pos) %>% as.data.frame() %>% set_colnames(c("pos","cnt"))  )) %>% unnest()
    }
  }

