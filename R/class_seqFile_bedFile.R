# operations associated with bed file

 .bedFile<- R6Class("bedFile",
                    inherit = .seqFile,
                    # private ----------------
                    private = list(
                      lookup_genome= NULL, #"BSgenome",
                      set_lookup_genome= function(){
                        private$lookup_genome= switch(self$genome,  # specify genome to extract sequence
                                                   hg19 = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                                                   hg19.masked = BSgenome.Hsapiens.UCSC.hg19.masked::BSgenome.Hsapiens.UCSC.hg19.masked,
                                                   mm9 = BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9,
                                                   mm9.masked = BSgenome.Mmusculus.UCSC.mm9.masked::BSgenome.Mmusculus.UCSC.mm9.masked,
                                                   stop("genome not found for the specified --- ",self$genome)
                        )
                      }
                    ),
                  # public ----------------
                    public = list(
                      file =  NULL, #"character",
                      bed =  NULL, #"data.frame",
                      seq =  NULL, #"Biostrings",
                      PWM =  NULL, #"matrix",
                      moodsResult =  NULL, #"ANY",
                      genome=  NULL, #"character",
                      from_Grange= NULL,

                      # inherited method: moodsMap

                      initialize = function(file, genome,from_Grange) {
                        self$file <- file
                        self$genome <- genome; private$set_lookup_genome() #; self$get_lookup_genome()
                        if (from_Grange)
                        {
                          self$bed = as.data.frame(file);colnames(self$bed)[1:3]= c("chr","start","end")
                        }else
                        {
                          self$bed <- readr::read_delim(file,"\t",col_names = FALSE,skip = 1)
                          colnames(self$bed)[1:3]= c("chr","start","end")
                          message(paste0("reading bed file", self$file))
                        }
                      },

                      # getSeq again if changed ranges
                      getSeq = function() {
                        self$seq= BSgenome::getSeq(private$lookup_genome, self$bed$chr, start=self$bed$start, end=self$bed$end)
                        names(self$seq)=paste0(1:length(self$seq),",", names(self$seq) )
                      },

                      get_lookup_genome=function() show(private$lookup_genome),

                      resizeRng = function(newWinSize=400, outOfRng.rm=TRUE){
                        self$bed = self$bed %>% dplyr::mutate(width=end-start+1, add=newWinSize-width, left= ceiling(add/2), right=floor(add/2), start=start-left,
                                                       end=end+right, width=end-start+1) %>% dplyr::select(-add,-left,-right)

                        if (outOfRng.rm){
                          chrLength= GenomeInfoDb::seqlengths(private$lookup_genome)-1
                          self$bed= self$bed %>% dplyr::filter(start>0 & end<=chrLength[chr])
                        }
                      }

                    )

                )
# constructor fun #------------
  bedFile <- function(file=character(),genome="hg19.masked",from_Grange=FALSE) .bedFile$new(file=file,genome=genome,from_Grange=from_Grange)



# test ----------------
# testBed=bedFile(file="/Volumes/Nutcase_DataZhu/Analysis/Nusel/ChIPseq/bedFile/Nature_2015/bed/GSM1505782_T_021113_h64.bed.peak.txt_new.txt")
# testBed$resizeRng(100000,outOfRng.rm = T); testBed$bed %>% show()
# testBed$getSeq()
# testBed$moodsMap(pwmFile="/Volumes/Nutcase_DataZhu/DefaultFilesforProc/PWM/PWM_for_Mapping/*.pfm", p="0.0001")

