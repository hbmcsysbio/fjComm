

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



align_AA<-function(string_vector,outDir=".",filename=NULL,phylo=FALSE)
{
  # align amino acid sequences, input is a character vector containing aa sequences
  pacman::p_load(msa,seqinr,ape)
  mySequences <- string_vector %>% AAStringSet()
  names(mySequences)=species_names

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


align_AA_UniprotID<-function(UniprotID_vector,outDir=".",filename=NULL,phylo=FALSE)
{
  seqs=getFASTAFromUniProt(UniprotID_vector,2)
  string_vector=seqs %>% stringr::str_extract("(?<=\n)(.*\n*){1,10}$") %>% {gsub("\n","",.)}
  align_AA(string_vector,outDir = outDir, filename = filename,phylo = phylo)
}

