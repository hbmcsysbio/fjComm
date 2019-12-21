# highest kmer pair MI against distance (for TF motif length and cross-strand binding)
# 2016-10-24
library(fjComm)
kmerLen=3

test=F


maxGap=ifelse(test,10,0)
minGap=0; maxGap=130
testFileAdd= ifelse(test,"_test","")
verbose=T

calcOverWrite=T
# maxGap=85; minGap=67

# parameters to use when running on Mac
args_mac=c(
         reads_file="./Rtest/Trulig200v1IIIFJ4-4-CAP2-P10-TF377-TIIIc4_S754_L002_R1_001.peared_trimmed.fq.gz",
         # reads_file="/Volumes/Nutcase_DataZhu/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_147/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-4-CAP2-A06-TF99-RFX5IIIc4_S390_L002_R1_001.peared_trimmed.fq.gz",

         c0File="./Rtest/subset_200000_Trulig200v1IIIPEx200IIIc0_S6_R1_001.peared_trimmed.fq.gz",
         outDir="Rtest"
         )

# head-----------
scriptName<-"kmer_pair_MI_dist.R"
    isMac=T
    RlibDir= ifelse(isMac, "/Volumes/Nutcase_DataZhu/lib/Rlib/", "/wrk/data/fangjie/lib/Rlib/")

#-----------
usage=paste(' ~/R/R320/bin/Rscript ', scriptName," '<readsfile>' '<c0file>' <outputDir> \n",sep = "")
    args=args_mac; if(!isMac) {args=commandArgs(trailingOnly=TRUE)}
    if (!isMac && length(args)==0) {cat("usage:  ",usage); stop() }
filePattern=args[1]
c0File= args[2]
outDir=ifelse(is.na(args[3]), paste("./",scriptName,"/",sep = ""), args[3])

#---------------------------------------------------
    system(paste("md ", outDir,sep = ""));
    system(paste("md ", outDir,"\\plot\\",sep = ""))
    system(paste("md ", outDir,"\\one_off_calc\\",sep = ""))
    system(paste("md ", outDir,"\\data\\",sep = ""))
    cmd= paste("echo 'Rscript ",scriptName," ", paste(args,collapse = " "),"' >",outDir,"\\Rcmd_",scriptName,".txt", sep = ""); system(cmd) # memo cmd
#---------------------end head----------------------#

targetDir=dirname(filePattern)
allfiles<- list.files(path=targetDir, pattern = basename(filePattern),full.names = F) #read in all txt and fastq files
tmpDir=paste(outDir,"\\tmp\\",sep = ""); system(paste("md ",tmpDir,sep = ""));
#--------------------end fixed part------------------


# source("http://bioconductor.org/biocLite.R")  # for nutcase server
if (!require(pacman,quietly = T)){ install.packages("pacman") }
pacman::p_load(readr,stringr,ggplot2,dplyr,reshape2,Rcpp)
# source(paste0(RlibDir,"/Mylib/ggplot2_common.R"))
# source(paste0(RlibDir,"/Mylib/fileIO.R")) # for oneOffCalc
# source(paste0(RlibDir,"/Mylib/seqFile.R")) # for dedup
# sourceCpp(paste0(RlibDir,"/common/kcnt.cpp"))


################
# functions ---------
gapkCount_c0 <-function(c0File)
{
  allreads_c0= read_table(c0File,col_names = F) %>% .[[1]]
  # if (test && (allreads_c0 %>% length() > 10000)) {allreads_c0=allreads_c0[1:10000]}

  pseudoCnt= min((length(allreads_c0)/4^(kmerLen*2)/20) %>% as.integer(),10)
  gapped_kCnt_c0=gkmerCntBit(allreads_c0,1L,k=kmerLen,gapMins = minGap,gapMaxs = maxGap,pseudo = pseudoCnt) %>% melt
  colnames(gapped_kCnt_c0)=c("kmer","gapLen","counts")
  gapped_kCnt_c0= gapped_kCnt_c0 %>% mutate(gapLen= sub("n","",gapLen),kmer1= str_sub(kmer,1,kmerLen),kmer2=str_sub(kmer,kmerLen+1,2*kmerLen)) %>% mutate(gapLen=as.integer(gapLen))

  gapped_kCnt_df_c0= gapped_kCnt_c0 %>%
    dplyr::group_by(kmer1,gapLen) %>% mutate(prop_act_c0= counts/sum(counts)) %>%
    ungroup
  return(gapped_kCnt_df_c0)
}

gapkCount <-function(file_full)
{

  allreads= read_table(file_full,col_names = F)
  allreads= rmdup(allreads,1) # dedup if running on server
  allreads= allreads[[1]]

  # if (test && (allreads %>% length() > 10000)) {allreads=allreads[1:10000]}
  readsLen= allreads[1] %>% nchar()

  pseudoCnt= min((length(allreads)/4^(kmerLen*2)/20) %>% as.integer(),10)


  allkCnt= kmerCntBit(allreads,kmerLen,collapse = T,pseudo = pseudoCnt,asDf = F)  #counts for half of the gapped kmer
  prop_allkCnt=allkCnt/sum(allkCnt)

  # bk from current cycle
  gapped_kCnt=gkmerCntBit(allreads,1L,k=kmerLen,gapMins = minGap,gapMaxs = maxGap,pseudo = pseudoCnt) %>% melt
  colnames(gapped_kCnt)=c("kmer","gapLen","counts")
  gapped_kCnt= gapped_kCnt %>% mutate(gapLen= sub("n","",gapLen),kmer1= str_sub(kmer,1,kmerLen),kmer2=str_sub(kmer,kmerLen+1,2*kmerLen)) %>% mutate(gapLen=as.integer(gapLen))

  # calc foldchn_VS_mindiff
  gapped_kCnt_df= gapped_kCnt %>%
    dplyr::group_by(kmer1, gapLen) %>%
    mutate(prop_act= counts/sum(counts)) %>%
    ungroup %>%
    mutate(prop_exp_curr= prop_allkCnt[kmer2], prop_exp_c0=gapped_kCnt_df_c0[match(.$kmer,gapped_kCnt_df_c0$kmer), ]$prop_act_c0) %>%
    mutate(min_diff_exp= ifelse(abs(prop_exp_curr-prop_act)> abs(prop_exp_c0-prop_act), prop_exp_c0, prop_exp_curr), foldchn_VS_mindiff= log2(prop_act/min_diff_exp)) %>%
    arrange(foldchn_VS_mindiff)
  return(gapped_kCnt_df)
}

# body ========

# bk from c0 (kmer proportions)
# saved_oneOff_filename= paste0("../one_off_calc/",basename(c0File),testFileAdd,"_k",kmerLen,"_maxGap",maxGap,".csv")
# gapped_kCnt_df_c0= oneOffCalc(saved_oneOff_filename, calcFun =gapkCount_c0 , calcFunParamList =list(c0File), overWrite = calcOverWrite)
gapped_kCnt_df_c0=gapkCount_c0(c0File)

# for (file_name in allfiles)
# {
file_name=allfiles[1];
  cat(paste("processing ",file_name,"\n",sep = ""));
  file_full= paste(targetDir, file_name,sep = "/") #full path

  # saved_oneOff_filename= paste0("/data/",basename(file_full),testFileAdd,"_k",kmerLen,"_maxGap",maxGap,".csv")
  # gapped_kCnt_df= oneOffCalc(saved_oneOff_filename, calcFun =gapkCount, calcFunParamList =list(file_full) ,overWrite = calcOverWrite)
  gapped_kCnt_df=gapkCount(file_full)
  # gapped_kCnt_df$foldchn_VS_mindiff=log2(gapped_kCnt_df$prop_act/gapped_kCnt_df$prop_exp_c0)

  topEachGap= gapped_kCnt_df %>% group_by(gapLen) %>% arrange(foldchn_VS_mindiff) %>% top_n(5,foldchn_VS_mindiff)
  bottomEachGap= gapped_kCnt_df %>% group_by(gapLen) %>% arrange(foldchn_VS_mindiff) %>% top_n(-5,foldchn_VS_mindiff)

  topPlot= topEachGap %>% select(gapLen,foldchn_VS_mindiff) %>% mutate(end=factor("top", levels=c("top", "bottom")))
  bottomPlot= bottomEachGap %>% select(gapLen,foldchn_VS_mindiff) %>% mutate(end=factor("bottom", levels=c("top", "bottom")))
  topBottomPlot= rbind(topPlot,bottomPlot) #%>% mutate(motifLen=gapLen+ kmerLen)
  # optional
  longest=(topBottomPlot$gapLen %>% max)
  topBottomPlot= topBottomPlot %>% dplyr::filter(gapLen<= longest)

  p=ggplot() + trans_bk_xyaxis + geom_point(data=topBottomPlot, aes(gapLen, foldchn_VS_mindiff, color=end), size=0.2, alpha=1/2) +
    scale_y_continuous(expand = c(0,0))+ scale_x_continuous(expand = c(0,0))+
    geom_vline(xintercept = c(10,20,30,40,50,60,70,80,90),size=0.4,color="grey")+ theme(legend.position = c(.90, .50))

  #plot mean line
  meanDf= topBottomPlot %>% group_by(gapLen, end) %>% summarise(meanfdchn=mean(foldchn_VS_mindiff)) #%>% mutate(gapLen=as.integer(gapLen))
  p=p + geom_line(data = meanDf, aes(gapLen, meanfdchn, color=end))
  print(p)

  img=paste0(outDir,"/plot/", extract_well(file_name), "_", file_name,".pdf")
  ggsave(filename =img ,plot = p,width = 4,height =3 )
  # system( paste0("inkscape -f ",img," -b white -d 200 --export-png=",img,".png") );

  print(paste0("finish curr file",file_name))
# }

