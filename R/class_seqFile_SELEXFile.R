# operations associated with SELEX reads

.SELEXFile<- R6Class("SELEXFile",
                     inherit = .seqFile,
                    # private ----------------
                     private = list(

                       # lookup_genome= NULL, #"BSgenome",
                       seqIni= function()
                       {
                         self$seqNum= nrow(self$seq); self$seqLen= nchar(self$seq[[1]][1])
                         self$seq= Biostrings::DNAStringSet(self$seq[[1]])
                         names(self$seq)=1:length(self$seq)
                       },

                       compare_c0 = function(target=NULL, pseudo=5L, ...)
                       { # target= 'kmerCnt'
                         if (is.null(target)) stop("specify target field to cmp for compare_c0")
                         c0= self$c0_SELEXFile_obj
                         if (is.null(c0)) stop("specify c0_SELEXFile_obj for current SELEXFile first")

                         # capture additional param
                         dots=list(...); for (i in seq_along(dots)) {assign(names(dots[i]),dots[[i]])}
                         if ((!count_c0) && is.null(get(target,envir = c0))) stop("count for c0 object first")

                         if (is.null(c0$seq)) c0$getSeq()
                         if (count_c0) {
                           if (target=="kmerCnt") {c0$kmerCnt= kmerCntBit(c0$seq %>% as.character(),k=k, collapse=collapse, diffLen=diffLen, asDf=asDf,pseudo = pseudo)}
                           else if (target=="kmerCntGap") {c0$kmerCntGap= gkmerCntBit(c0$seq %>% as.character(), gapNo = gapNo, k=k, gapMins = gapMins, gapMaxs = gapMaxs,diffLen = diffLen,posInfo = posInfo,all_possible_k = all_possible_k,pseudo = pseudo);  c0$kmerCntGap= c0$kmerCntGap %>% reshape2::melt() %>% dplyr::rename(kmer=Var1,gapLen=Var2,counts=value) %>% mutate(gapLen= sub("n","",gapLen),kmer1= stringr::str_sub(kmer,1,k),kmer2=stringr::str_sub(kmer,k+1,2*k)) %>% mutate(gapLen=as.integer(gapLen)) }
                           else {stop("not supported yet")}
                         }

                         pseudofreq= pseudo/self$seqNum

                         if (target=="kmerCnt"){
                           tmpCnt= get(target,envir = c0) %>% select(kmer,counts)
                           colnames(tmpCnt)[2]= "c0counts"
                           assign(target,
                                  merge(get(target,envir = self), tmpCnt,by = "kmer") %>%
                                    mutate(freq= counts/sum(counts), c0freq= c0counts/sum(c0counts), foldchn= log2((freq - c0freq) / (c0freq + pseudofreq) + 1)) %>%
                                    select(-freq,-c0freq),
                                  envir = self)
                         }else if (target=="kmerCntGap"){
                           tmpCnt= get(target,envir = c0) %>% select(kmer,counts,gapLen)
                           colnames(tmpCnt)[2]= "c0counts"
                           assign(target,
                                  merge(get(target,envir = self), tmpCnt,by = c("kmer","gapLen")) %>%
                                    mutate(freq= counts/sum(counts), c0freq= c0counts/sum(c0counts), foldchn= log2((freq - c0freq) / (c0freq + pseudofreq) + 1)) %>%
                                    select(-freq,-c0freq),
                                  envir = self)
                           }
                         else {stop("not supported yet")}
                         print("foldchn is log2(foldchn)")
                       }


                     ),

                     # public ----------------
                     public = list(
                       file = NULL,
                       seqNum = NULL,
                       seqLen = NULL,
                       seq =  NULL, #Biostring
                       moodsResult=NULL,
                       c0_SELEXFile_obj=NULL,
                       pfmFiles= NULL,
                       well=NULL,
                       # inherited: kmerCnt, kmerCntGap
                       # inherited method: moodsMap

                       initialize = function(file, c0_SELEXFile_obj) {
                         self$file <- file
                         self$c0_SELEXFile_obj <- c0_SELEXFile_obj
                       },

                       # get seq with or without rm duplicates
                       getSeq = function(dup_rm = FALSE)
                       {
                         # if(dup_rm)
                         # {
                         #   message(paste0("reading SELEX reads file (dedup): ", self$file))
                         #   self$seq= dedup(self$file, tmpDir= "~/_tmp/",keepFile=FALSE)
                         #   private$seqIni()
                         # }
                         # else
                         # {
                           message(paste0("reading SELEX reads file: ", self$file))
                           self$seq= readr::read_tsv(self$file,col_names = F)
                           if(dup_rm){self$seq= rmdup(seq_df= self$seq, 1)}
                           private$seqIni()
                         # }
                       },

                       # count kmer
                       count_k = function(k=4L,collapse=TRUE,diffLen=FALSE, asDf= TRUE, pseudo=5L, all_possible_k=T, cmp_c0=FALSE, count_c0=TRUE, rc=FALSE, raw=FALSE) {
                         if (is.null(self$seq)) stop("getSeq first before count k")
                         self$kmerCnt= kmerCntBit(self$seq %>% as.character(),k=k, collapse=collapse, diffLen=diffLen, asDf=asDf, pseudo = pseudo, all_possible_k = all_possible_k )
                         if(raw){return()}
                         if (rc)
                         {
                           if(collapse)
                           {
                             rckmer= (Biostrings::DNAStringSet(self$kmerCnt$kmer) %>% Biostrings::reverseComplement() %>% as.character)
                             # browser()
                             rcCnt= self$kmerCnt %>% mutate(kmer= rckmer)
                             self$kmerCnt$counts = rcCnt[match(self$kmerCnt$kmer,rcCnt$kmer), ]$counts + self$kmerCnt$counts
                           }
                           else # with pos info
                           {
                              rckmer=revComp(colnames(self$kmerCnt))
                              rcCnt= self$kmerCnt %>% set_colnames(rckmer) %>% .[,colnames(self$kmerCnt)]
                              self$kmerCnt=self$kmerCnt+rcCnt
                           }

                         }

                         if (cmp_c0) private$compare_c0(target='kmerCnt',pseudo=pseudo,k=k,collapse=collapse,diffLen=diffLen, asDf= asDf, count_c0=count_c0)
                       },


                       # count gapped kmer
                       count_gpk = function(k=3L, gapNo = 1,  gapMins = 0, gapMaxs = 100, diffLen = FALSE, posInfo = FALSE, all_possible_k = TRUE, pseudo = 5L, cmp_c0=FALSE, count_c0=TRUE, raw=FALSE) {
                         if (is.null(self$seq)) stop("getSeq first before count k")
                         self$kmerCntGap= gkmerCntBit(self$seq %>% as.character(),gapNo = gapNo, k=k, gapMins = gapMins, gapMaxs = gapMaxs,diffLen = diffLen,posInfo = posInfo,all_possible_k = all_possible_k,pseudo = pseudo)
                         if(raw){return()}
                            if (gapNo>1) {self$kmerCntGap= self$kmerCntGap %>% reshape2::melt() %>% dplyr::rename(kmer=Var1,gapLen=Var2,counts=value) }# does not work for 2 or more gap
                              else self$kmerCntGap= {self$kmerCntGap %>% reshape2::melt() %>% dplyr::rename(kmer=Var1,gapLen=Var2,counts=value) %>% mutate(gapLen= sub("n","",gapLen),kmer1= stringr::str_sub(kmer,1,k),kmer2=stringr::str_sub(kmer,k+1,2*k)) %>% mutate(gapLen=as.integer(gapLen))}
                         if (cmp_c0) private$compare_c0(target='kmerCntGap',gapNo = gapNo, k=k, gapMins = gapMins, gapMaxs = gapMaxs,diffLen = diffLen,posInfo = posInfo,all_possible_k = all_possible_k,pseudo = pseudo,  count_c0=count_c0)
                       },

                       getTFandWell =function()
                       {
                         self$TF= extract_TF(self$file)
                         self$well= extract_well(self$file)
                       }

                     )

)
# constructor fun #------------
SELEXFile <- function(file, c0_SELEXFile_obj=NULL) .SELEXFile$new(file=file, c0_SELEXFile_obj=c0_SELEXFile_obj)

# reads=SELEXFile("~/Nut_zhuData/seqFiles2/FJ4.4CAP2_PE/allreads/FJ4.4_PE_147/dual_trim/2_Trim_adaptor_kmer/Trulig147v1IIIFJ4-4-CAP2-A06-TF99-RFX5IIIc4_S390_L002_R1_001.peared_trimmed.fq.gz")
# reads$getSeq(dup_rm = F)
