TF_classify <- function(TFnames)
{
  classify_table= read_csv("~/Nut_zhuData/DefaultFilesforProc/TF_classification/TFclassify_table.csv")
  classify_table[match(TFnames,classify_table$TF),]$family
  # classify_table %>% dplyr::slice(match(TFnames,classify_table$TF)) %>% .$family
}

TF_color_acc_family <- function(TFfamily)
{
  classify_table= read_csv("~/Nut_zhuData/DefaultFilesforProc/TF_classification/TFclassify_table.csv")
  # family=classify_table[match(TFnames,classify_table$TF),]$family
  # classify_table %>% dplyr::slice(match(TFnames,classify_table$TF)) %>% .$family
  colors=table(classify_table$family); colors=rainbow(length(colors)) %>% set_names(names(colors))
  colors["CUT"]="#CCCC00FF"
  colors[TFfamily]
}
