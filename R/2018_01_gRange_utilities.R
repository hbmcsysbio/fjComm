# vplot

# bamfile="../../../allreads/bam/FJ49MChIP-ChIP-*D04*.bam" %>% Sys.glob()
# motifHits=readRDS("../0_pos_of_motifhits/D04_RFX5_hg19_motifHits.RDS"); motifLen=ncol(motifHits)
# motifGR=GRanges(seqnames = motifHits$chr,ranges = IRanges(start = motifHits$hitstart,end=motifHits$hitend),strand = motifHits$hitstrand)
# GrngPE <- readGAlignmentPairs(bamfile,use.names=TRUE, strandMode=0,param = ScanBamParam(which = motifGR_wide )) %>% granges(on.discordant.seqnames="drop")


# calc_Vplot <- function(PEReadsBam="*.bam",sitesGR="GR with str info", regionWidth=2001, for_insert_size=FALSE)
#   # sites has strand info
# {
#   pacman::p_load(GenomicRanges,GenomicAlignments,Rsamtools)
#   sitesGR_1bp=sitesGR %>%GenomicRanges::resize(width = 1,fix = "center")
#   sitesGR_wide= sitesGR %>% GenomicRanges::resize(width = regionWidth,fix = "center")
#   nearbyReadsGR <- readGAlignmentPairs(PEReadsBam,use.names=F, strandMode=0,param = ScanBamParam(which = sitesGR_wide,what = "isize"))
#     if(for_insert_size) return(mcols(nearbyReadsGR@first)[[1]] %>% abs)
#   nearbyReadsGR %<>% granges(on.discordant.seqnames="drop")
#   nearbyReadsGR_1bp= nearbyReadsGR %>%GenomicRanges::resize(width = 1,fix = "center")
#
#   nearest_=GenomicRanges::nearest(nearbyReadsGR_1bp,subject = sitesGR_1bp)
#   nearest_site=sitesGR_1bp[nearest_]
#   dist_= start(nearbyReadsGR_1bp)-start(nearest_site)
#   direc_=ifelse(strand(nearest_site)=="+",1,-1)
#   dist_=dist_*direc_
#   data.frame(center=dist_,lens=width(nearbyReadsGR))
#   # data.frame(center=dist_,lens=insert_sizes)
# }


grange_PE2SR_withsize<-function(PErange)
  # PE to SR grange, with isize retained
  # PE read with readGAlignmentPairs(PEReadsBam,use.names=F, strandMode=0,param = ScanBamParam(which = sitesGR_wide,what = "isize"))
{
  PErange=PErange[PErange@first@seqnames==PErange@last@seqnames]
  lens=PErange@first %>% mcols() %>% .[["isize"]] %>% abs()
  PErange %<>% granges() #granges(on.discordant.seqnames="drop",use.mcols = T)
  PErange$lens=lens
  PErange
}

calc_Vplot<-function(PEReadsBam="*.bam",sitesGR="GR with str info", regionWidth=2001, for_insert_size=FALSE)
  # sites has strand info
{
  pacman::p_load(GenomicRanges,GenomicAlignments,Rsamtools)
  sitesGR_1bp=sitesGR %>%GenomicRanges::resize(width = 1,fix = "center")
  sitesGR_wide= sitesGR %>% GenomicRanges::resize(width = regionWidth,fix = "center")
  nearbyReadsGR <- readGAlignmentPairs(PEReadsBam,use.names=F, strandMode=0,param = ScanBamParam(which = sitesGR_wide,what = "isize"))
  if(for_insert_size) return(mcols(nearbyReadsGR@first)[[1]] %>% abs)


  nearbyReadsGR %<>% grange_PE2SR_withsize()

  nearbyReadsGR_1bp= nearbyReadsGR %>%GenomicRanges::resize(width = 1,fix = "center")

  nearest_=GenomicRanges::nearest(nearbyReadsGR_1bp,subject = sitesGR_1bp)
  nearest_site=sitesGR_1bp[nearest_]
  dist_= start(nearbyReadsGR_1bp)-start(nearest_site)
  direc_=ifelse(strand(nearest_site)=="+",1,-1)
  dist_=dist_*direc_
  data.frame(center=dist_,lens=nearbyReadsGR_1bp %>% mcols() %>% .[["lens"]])
  # data.frame(center=dist_,lens=insert_sizes)
}


