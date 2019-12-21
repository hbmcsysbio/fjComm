add_file_col<-function(df,files,pattern_b="deg-",pattern_a="-III",new_col_name="c4f_0deg")
{
  TFs=files %>% str_extract(paste0("(?<=",pattern_b,").*(?=",pattern_a,")"))
  files=files %>% str_replace("/Users/zhu","~")
  orders=match(df$TF,TFs)
  df[[new_col_name]] =files[orders]
  if (!identical(df[[new_col_name]] %>% str_extract(paste0("(?<=",pattern_b,").*(?=",pattern_a,")")),df$TF )) stop("TF match err!!")
  df
}
