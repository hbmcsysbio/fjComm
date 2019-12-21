# seq= read_csv("~/Nut_zhuData/Analysis/Analysis2/kcnt_test/testseq.txt",col_names = F)
# #
# for (i in 7:7){
# kcnt1=kmerCnt(seq$X1,k = i,diffLen = F,collapse = T,asDf = F)
# kcnt2=kmerCntBit(seq$X1,k = i,diffLen = F,collapse = T,asDf = F,all_possible_k = T)
# print(i)
# # print(all_equal(kcnt1,kcnt2)) #for non-collapse
# print(identical(kcnt1,kcnt2))
# }

# c0=SELEXFile("~/Nut_zhuData/seqFiles2/c0/subset_100000_Trulig147v1IIIPEx200IIIc0_S3_R1_001.peared_trimmed.fq.gz")
# # c0$getSeq()
# # c0$count_gpk()
# test=SELEXFile("~/Nut_zhuData/Analysis/Analysis2/kcnt_test/testseq.txt",c0_SELEXFile_obj = c0)
# test$getSeq()
# test$count_gpk(count_c0 = T,cmp_c0 = T)




# library(fjComm)
# bk=SELEXFile("~/Nut_zhuData/seqFiles/FAC_PE_EMSA_FJ003_0,25/FAC_Cyc4/lane8_NoIndex_L008_Assembly.fastq.assembled.fastq_barcodes/JustReads/FAC_14706CATG94N_4.txt"); bk$getSeq();
# sig=SELEXFile("~/Nut_zhuData/seqFiles/FAC_PE_EMSA_FJ003_0,25/FAC_Cyc4/lane8_NoIndex_L008_Assembly.fastq.assembled.fastq_barcodes/JustReads/FAC_14704GGTC94N_4.txt"); sig$getSeq()
# targetSeqsB=read_csv("~/Nut_zhuData/seqFiles/FAC_PE_EMSA_FJ003_0,25/FAC_Cyc4/lane8_NoIndex_L008_Assembly.fastq.assembled.fastq_barcodes/JustReads/FAC_14705AGGA94N_4.txt",col_names = F)$X1;
# targetSeqsUB=read_csv("~/Nut_zhuData/seqFiles/FAC_PE_EMSA_FJ003_0,25/FAC_Cyc4/lane8_NoIndex_L008_Assembly.fastq.assembled.fastq_barcodes/JustReads/FAC_14707GTAG94N_4.txt",col_names = F)$X1;
#
#
#
#   bk$count_gpk(gapMaxs = 20,raw=T,posInfo = T)
#   sig$count_gpk(gapMaxs = 20,raw=T,posInfo = T)
#   scoreCube= log2( (sig$kmerCntGap$result / sig$seqNum) / (bk$kmerCntGap$result / bk$seqNum) )
#   dimnames(scoreCube)=sig$kmerCntGap$dimnames
#
#
#   g_scoresB=scoring(targetSeqsB,scoreCube, k=3L, gapNo = 1,  gapMins = 0, gapMaxs = 20, diffLen = FALSE, posInfo = TRUE)
#   g_scoresUB=scoring(targetSeqsUB,scoreCube, k=3L, gapNo = 1,  gapMins = 0, gapMaxs = 20, diffLen = FALSE, posInfo = TRUE)
#
#   g_scoresB_noPos=scoring(targetSeqsB,scoreCube, k=3L, gapNo = 1,  gapMins = 0, gapMaxs = 20, diffLen = FALSE, posInfo = FALSE)
#   g_scoresUB_noPos=scoring(targetSeqsUB,scoreCube, k=3L, gapNo = 1,  gapMins = 0, gapMaxs = 20, diffLen = FALSE, posInfo = FALSE)
#
#
#
#   bk$count_gpk(gapMaxs = 0,raw=T,posInfo = T)
#   sig$count_gpk(gapMaxs = 0,raw=T,posInfo = T)
#   scoreCube= log2( (sig$kmerCntGap$result / sig$seqNum) / (bk$kmerCntGap$result / bk$seqNum) )
#   dimnames(scoreCube)=sig$kmerCntGap$dimnames
#
#   scoresB=scoring(targetSeqsB,scoreCube, k=3L, gapNo = 1,  gapMins = 0, gapMaxs = 0, diffLen = FALSE, posInfo = TRUE)
#   scoresUB=scoring(targetSeqsUB,scoreCube, k=3L, gapNo = 1,  gapMins = 0, gapMaxs = 0, diffLen = FALSE, posInfo = TRUE)
#
#   scoresB_noPos=scoring(targetSeqsB,scoreCube, k=3L, gapNo = 1,  gapMins = 0, gapMaxs = 0, diffLen = FALSE, posInfo = FALSE)
#   scoresUB_noPos=scoring(targetSeqsUB,scoreCube, k=3L, gapNo = 1,  gapMins = 0, gapMaxs = 0, diffLen = FALSE, posInfo = FALSE)
#
#   boxplot(g_scoresB,g_scoresUB,g_scoresB_noPos,g_scoresUB_noPos, scoresB,scoresUB,scoresB_noPos,scoresUB_noPos, names = c("g_b_pos","g_ub_pos","g_b","g_ub","b_pos","ub_pos","b","ub"), boxwex=0.7 )




