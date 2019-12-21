# operations associated with sequences

.seqFile<- R6Class("seqFile",
                     # private ----------------
                     private = list(

                       # lookup_genome= NULL, #"BSgenome",
                       # seqIni= function()
                       # {
                       #   self$seqNum= nrow(self$seq); self$seqLen= nchar(self$seq[[1]][1])
                       #   self$seq= Biostrings::DNAStringSet(self$seq[[1]])
                       # }
                       #
                     ),

                     # public ----------------
                     public = list(
                         file = NULL,
                         seq =  NULL, #Biostring
                         moodsResult=NULL,
                         kmerCnt= NULL,
                         kmerCntGap= NULL,
                         TF=NULL,
                         pfmFiles=NULL,

                         initialize = function(file) {
                           self$file <- file
                         },

                        getSeq= function(){},

                        # mapping using moods cmd
                        moodsMap = function(pwmFile=c("test/*.{pfm,adm}","*.pfm"), p="0.0001", batch= TRUE,isMatrix=F) {
                           # -m example-data/matrices/*.{pfm,adm} -s example-data/seq/chr1-5k-55k.fa -p 0.0001  # -m can follow multiple files

                           if (!is.character(p)) p=as.character(p)#stop("p= character, not number!")
                           if (is.null(self$seq)) stop("getSeq first")
                           system("mkdir -p ~/_tmp")
                           rndSuffix= paste(sample(letters,40,replace = T),collapse = "")
                           tmpFileName= paste0("~/_tmp/tmp_", ifelse(is.character(self$file),basename(self$file),tempfile() %>% basename()),"_",rndSuffix,".fa")
                           tmpOutFile=paste0(tmpFileName,"_o")
                              if(isMatrix){tmppwmFile=paste0("~/_tmp/pwm_",tempfile() %>% basename(),".pfm"); write_tsv(as.data.frame(pwmFile),tmppwmFile,col_names = F); pwmFile=tmppwmFile}

                           Biostrings::writeXStringSet(self$seq, filepath = tmpFileName, format="fasta") #mapping
                           cmd=paste(sep =" ", "moods_dna.py -m ",paste(pwmFile,collapse = " ")," -s ",tmpFileName, " -p ",p, ifelse(batch," --batch "," "), "-o", tmpOutFile)
                           system(cmd)
                           self$moodsResult=readr::read_csv(tmpOutFile, col_names = F) %>% dplyr::rename(rangeNo=X1,pfmFile=X2,pos=X3,strand=X4,score=X5,match=X6) #%>% select(-X7)

                           system(paste(sep =" ", "rm",tmpFileName,tmpOutFile))
                              if(isMatrix){system(paste(sep =" ", "rm",tmppwmFile))}
                        },

                        getPfmFiles = function(guidefiles, TFname_full=self$TF, monomericFirst=TRUE, exactMatchFirst=TRUE, head=20)
                        {
                          # if (is.null(TFname_full)) TFname_full=self$TF
                          if(is.null(TFname_full)) stop("assign the TF name first")
                          TFname_noNum= str_match(TFname_full,"[A-Z]*")[1]
                          allpfms=character(0)
                          for(guidefile in guidefiles)
                          {
                            cmd=  paste0("cat ",guidefile," | grep -i '^",TFname_noNum,"[0-9A-Z-]*\t' | ") #some TF name like HOXB2L2 NKD2-4
                            # if (monomeric) cmd= paste0(cmd,"grep 'mono' | ")
                            cmd= paste0(cmd,'perl -lne \'@all=split /\\t/; $pfm_key="$all[2]_$all[3]_$all[4]"; @allpfm= glob("~/Nut_zhuData/SELEXphp/data/pfm/shortspacek/$pfm_key*.pfm"); chomp @allpfm; print "$all[0]_$all[3]_$all[4]_$all[7]_$all[8]"."\\t".$allpfm[0] if @allpfm;\'') # some like TAGTTA40NGTG_XEPAQ_NNNNTGCTGAC
                            allpfms= c(allpfms,system(cmd,intern = T))
                          }
                          allpfms= strsplit(allpfms,"\t") %>% as.data.frame() %>% t %>% as.data.frame() %>% dplyr::rename(pfmName=V1,file=V2)
                          if (monomericFirst){ allpfms= rbind((allpfms %>% dplyr::filter(grepl("mono",pfmName))),(allpfms %>% dplyr::filter(!grepl("mono",pfmName)))) }
                          if (exactMatchFirst){ allpfms= rbind((allpfms %>% dplyr::filter(grepl(TFname_full,pfmName))), (allpfms %>% dplyr::filter(!grepl(TFname_full,pfmName)))) }
                          allpfms= head(allpfms,head)
                          self$pfmFiles=allpfms %>% mutate(file=as.character(file),pfmName=as.character(pfmName))
                        }


                        # # count kmer
                        # count_k = function(k=4L,collapse=TRUE,diffLen=TRUE, asDf= TRUE, cmp_bk=FALSE) {
                        #   if (is.null(self$seq)) stop("getSeq first before count k")
                        #   self$kmerCnt= kmerCnt(self$seq %>% as.character(),k=k, collapse=collapse, diffLen=diffLen, asDf=asDf)
                        # },
                        #
                        # # count gapped kmer
                        # count_gpk = function(k=3L,diffLen=TRUE,maxGap=0, minGap=0, asDf= TRUE, cmp_bk=FALSE) {
                        #   if (is.null(self$seq)) stop("getSeq first before count k")
                        #   self$kmerCntGap= kmerCnt_allgap(self$seq %>% as.character(), k=k,maxGap = maxGap,diffLen = diffLen,minGap = minGap, asDf=asDf)
                        # }

                     )

)
# constructor fun #------------
seqFile <- function(file) .seqFile$new(file=file)
