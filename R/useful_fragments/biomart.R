# get TSS sites
pacman::p_load(biomaRt,fjComm)
ensembl=useMart(host='may2012.archive.ensembl.org',
                biomart='ENSEMBL_MART_ENSEMBL',
                dataset='mmusculus_gene_ensembl') # NCBIM37
# listDatasets(ensembl)
# filters = listFilters(ensembl)
# attributes = listAttributes(ensembl)

chrFilter= c(1:19,"X","Y") # only normal chrs
enQuery= fjComm::oneOffCalc("enQuery",biomaRt::getBM, list(attributes=c('chromosome_name', 'transcript_start'), #'name_1006'
                                                           filters = 'chromosome_name',
                                                           values = chrFilter,
                                                           mart = ensembl) ,asObj = T, overWrite = ifelse(test,F,T))
enQuery= enQuery %>% dplyr::rename(chr=chromosome_name,tss=transcript_start) %>% mutate(chr=paste0("chr",chr)) # GO=name_1006,
enQuery=split(enQuery,enQuery$chr)

# calculate tritile gcProb distClosestTSS
addFeatures <- function(result_obj)
{
  ## process result scored by DHS kmer/ ATI kmer
  pred= readRDS(result_obj)
  seg_pred= pred$GR %>% as_data_frame()
  # rm the last record which is out of bound
  seg_pred= seg_pred %>% group_by(seqnames) %>% dplyr::filter(end!=max(end)) %>% ungroup()
  # rm non DHS
  DHSseg_pred= seg_pred %>% dplyr::filter(DHSLabel==1) %>% arrange(desc(score)) %>% mutate(rank=rank(score))
  # devide into tritile
  num_per_tritile=(nrow(DHSseg_pred)/3) %>% floor()
  DHSseg_pred= DHSseg_pred %>% mutate(nrow1=1:nrow(DHSseg_pred)) %>% mutate(tritile= case_when(
    .$nrow1<= num_per_tritile ~ "tritile_1",
    .$nrow1>num_per_tritile & .$nrow1 <=num_per_tritile*2 ~ "tritile_2",
    TRUE ~ "tritile_3"
  ) )


  # GC content
  pacman::p_load(Biostrings,BSgenome.Mmusculus.UCSC.mm9.masked)

  # letterProb=with(DHSseg_pred, BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm9.masked, names=seqnames, start=start, end=end)) %>% letterFrequency("ACGT",OR = 0 ,as.prob = T)
  # gcProb=letterProb[,"G"]+letterProb[,"C"]

  letterProb= with(DHSseg_pred, BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm9.masked, names=seqnames, start=start, end=end)) %>% dinucleotideFrequency(as.prob = T)
  gcProb= letterProb[,"CG"]
  totalProb= rowSums(letterProb) # "NNNNNNNNNNN" will give total prob "0", should be filtered out
  DHSseg_pred = DHSseg_pred %>% mutate(gcProb=gcProb,totalProb=totalProb) %>% dplyr::filter(totalProb!=0)

  pacman::p_load(matrixStats)
  # as.character is important, otherwise not able to take the correct item in enQuery list
  # take outer substration, each row correspond to each row of DHSseq_pred
  DHSseg_pred=DHSseg_pred %>% group_by(seqnames) %>% mutate(distClosestTSS=outer((start+end+1)/2, enQuery[[first(seqnames) %>% as.character()]]$tss, "-") %>% abs %>% rowMins()) %>% ungroup()
  return(DHSseg_pred %>% dplyr::select(-nrow1,-totalProb,-strand,-width,-DHSLabel))
}
