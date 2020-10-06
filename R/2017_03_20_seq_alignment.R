

align_AA_DBD <- function(TFnames=c("DLX1","DLX2","DLX3"),outDir="~/Nut_zhuData/Analysis/misc/protein_alignment/test/",phylo=FALSE)
{
  pacman::p_load(msa,seqinr,ape,readxl)
  seqTableFile="~/Desktop/Experiment/Resourse/!!TF clones/!!ORF_arr_all/Fangjie_all_ORF_info_3ver.xlsx"
  seqTable= read_excel(seqTableFile, 5)
  seqTable= seqTable[!is.na(seqTable[[4]]),]

  table_subset=seqTable[match(TFnames, seqTable[[4]]),]

  mySequences <- AAStringSet(table_subset$`AA-sequence`)
  names(mySequences)=TFnames

  seqAlign=msa(mySequences)

  system(paste0("mkdir -p ",outDir))
  setwd(outDir)

  msaPrettyPrint(seqAlign, output="pdf",
                 logoColors="rasmol", shadingMode="similar",showLogo = "none",
                 showLegend=FALSE, askForOverwrite=FALSE)

  if (phylo)
  {
    seqinr= msaConvert(seqAlign, type="seqinr::alignment")
    d <- dist.alignment(seqinr, "identity")
    hemoTree <- nj(d)
    pdf(paste0(get_scriptpath(),"/tree.pdf"))
    plot(hemoTree, main="Phylogenetic Tree")
    dev.off()
  }

  pacman::p_load(Biobase)
  Biobase::openPDF(paste0(outDir,"/seqAlign.pdf"), bg=TRUE)
}



align_AA<-function(string_vector,outDir=".",seqNames=NULL,filename=NULL,phylo=FALSE)
{
  # align amino acid sequences, input is a character vector containing aa sequences
  pacman::p_load(msa,seqinr,ape)
  mySequences <- string_vector %>% AAStringSet()
  names(mySequences)=seqNames

  seqAlign=msa(mySequences,order = "input")

  system(paste0("mkdir -p ",outDir))
  setwd(outDir)

  msaPrettyPrint(seqAlign, output="pdf",file = filename %>% paste0(".pdf"),
                 logoColors="rasmol", shadingMode="similar",showLogo = "none",
                 showLegend=FALSE, askForOverwrite=FALSE)

  if (phylo)
  {
    seqinr= msaConvert(seqAlign, type="seqinr::alignment")
    d <- dist.alignment(seqinr, "identity")
    hemoTree <- nj(d)
    pdf(paste0(outDir,"/",filename,"_tree.pdf"))
    plot(hemoTree, main="Phylogenetic Tree")
    dev.off()
  }

}


align_AA_UniprotID<-function(UniprotID_vector,outDir=".",filename=NULL,seqNames=NULL,phylo=FALSE)
{
  pacman::p_load(Rcpi)
  if(!is.null(seqNames))seqNames=seqNames[UniprotID_vector!="NA"]
  UniprotID_vector=UniprotID_vector[UniprotID_vector!="NA"]
  seqs=getFASTAFromUniProt(UniprotID_vector,2)
  string_vector=seqs %>% stringr::str_extract("(?<=\n)(.*\n*){1,10}$") %>% {gsub("\n","",.)}
  align_AA(string_vector,outDir = outDir,seqNames = seqNames, filename = filename,phylo = phylo)
}



species_names=qw("human xenopus tomato arabidopsis rice maize")
hist_H2A=qw("Q6FI13 P06897 P25469 O23628 Q6ZL43 P40280")
hist_H2B=qw("P62807 P02281 O65821 Q9LQQ4 A3AGM4 P30755")
hist_H3=qw("P68431 P84233 NA P59226 Q0JCT1 P69246")
hist_H4=qw("P62805 P62799 P35057 P59259 NA P62787")


align_AA_UniprotID(UniprotID_vector = hist_H2A,filename = "hist_H2A",seqNames=species_names)
align_AA_UniprotID(UniprotID_vector = hist_H2B,filename = "hist_H2B",seqNames=species_names)
align_AA_UniprotID(UniprotID_vector = hist_H3,filename = "hist_H3",seqNames=species_names)
align_AA_UniprotID(UniprotID_vector = hist_H4,filename = "hist_H4",seqNames=species_names)


