biomart_info_from_id<-function(filter_val=c("ENSG00000101187", "ENSG00000167772", "ENSG00000105173"),filter='ensembl_gene_id',retrive=c('ensembl_gene_id','hgnc_symbol','wikigene_name','uniprot_gn_symbol','chromosome_name'))
{
  pacman::p_load(biomaRt)
  ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  biomaRt::getBM(attributes=retrive,filters = filter,values = filter_val,mart = ensembl)
}
