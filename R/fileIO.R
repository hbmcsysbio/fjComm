# read and write files

# write files for PHP, under outDir
outputPHP <- function ()
{

}

extract_well <- function(file_nameOnly,keep0=FALSE)
{
  library(pacman);  pacman::p_load(stringr)
  if (keep0) {file_nameOnly=sub("(-[A-Q]0{0,1})([[:digit:]]{1,2}-)","\\1\\2",file_nameOnly)} else {file_nameOnly=sub("(-[A-Q])0{0,1}([[:digit:]]{1,2}-)","\\1\\2",file_nameOnly)}
  well=stringr::str_match(file_nameOnly,"-([A-Q][[:digit:]]{1,2})-")[2]
  well
}
extract_batch <- function(file) (file %>% str_match("III(FJ\\d.\\d)-"))[,2] %>% str_replace("-",".")


extract_TF<-function(file_nameOnly,noNum=FALSE,sufix_rm=F)
{
  library(pacman);  pacman::p_load(stringr)

  # FJ11
  TF=str_extract(file_nameOnly,"(?<=-s-TF-).*(?=III)")


  if (is.na(TF)) {TF=str_match(file_nameOnly,"-([A-Z]*[-[0-9]]{0,2}[A-Z]{0,2}[0-9]{0,2})-?IIIc[[:digit:]]")[2]}
  # adj for E2F4, T
  if (is.na(TF)){TF=str_match(file_nameOnly,"-([A-Z]*[[0-9]]*)-?IIIc[[:digit:]]")[2] }

  # FJ4.5
  if (is.na(TF)) {TF=str_match(file_nameOnly,"-[A-P]{1}[0-9]{1,2}-(.*?)-?IIIc[[:digit:]]")[2]}

  if(noNum) {TF= str_match(TF,"[A-Z]*")[1]}
  if (sufix_rm) TF= str_replace(TF,"(.*?)ooo(.*)","\\1")
  TF
}



extract_cyc <- function(file) file %>% str_match("IIIc(\\d)(_|-)") %>% .[1,2] %>% as.integer()

rm_suffix <- function(filename) {stringr::str_match(basename(filename),"(.*)\\.[^\\.]*$")[2]}
add_well <- function(filename)
{
  currBase=basename(filename)
  well=extract_well(currBase)
  paste0(well,"_",currBase)
}

# obtain C12_dslfjsldjf
add_well_rm_suff <- function(filename)
{
  filename=rm_suffix(filename)
  add_well(filename)
}

save_gg <- function( ggImage, subfolder="/plot/", currfileName, add= "p1", suffix=".png", width=NA, height=NA, scale=1, outDir=NA)
{
  if (is.na(outDir)) {outDir=parent.frame()$outDir} # use the outDir var of the calling Env
  system(paste0("mkdir -p ", outDir, subfolder))
  outFileName= paste0(add_well_rm_suff(currfileName),"_",add,suffix)
  ggsave(plot=ggImage, filename =outFileName, path = paste0(outDir,subfolder), width = width,height = height,scale = scale)
}




# asObj=T if save obj directly
# do oneOff calc and save, or read from file if exist
# "*" OK in filename, default Fun is OK for dataframe
oneOffCalc <- function(saved_oneOff_filename, calcFun, calcFunParamList=list(), asObj=TRUE, overWrite=FALSE, saveFun=dfSaveFun, readFileFun=dfReadFileFun, useScriptPath=FALSE, outDir="")
{
  targetName= if (outDir!="") paste0(outDir,"/", saved_oneOff_filename) else saved_oneOff_filename
     script_dir=ifelse(exists("script_path_from_fun"),script_path_from_fun,fjComm::get_scriptpath())

  if ((dirname(targetName)=="." || useScriptPath) && (outDir==""))  {targetName=paste0(script_dir,"/",saved_oneOff_filename)} # specify if not specified dir name
  if (!dir.exists(dirname(targetName))) dir.create(dirname(targetName),recursive = TRUE) # create target dir if not exist

  saved_file=Sys.glob(targetName)

  if (asObj) {saveFun=saveRDS; readFileFun=readRDS}

  if (length(saved_file)==0 || overWrite)
  {
    print("new oneOffCalc >>>>> to "); print(targetName)
    calcResult= do.call(calcFun, calcFunParamList)
    system(paste0("mkdir -p ",dirname(targetName)))
    saveFun(calcResult, targetName)
  }else
  {
    print(paste0("oneOffCalc read from: ",targetName))
    calcResult= readFileFun(saved_file)
  }
  return(calcResult)
}

        dfSaveFun<-function(df, filename)
        {
          write.csv(df, filename, quote = F, row.names = F)
        }

        dfReadFileFun<-function(file)
        {
          df=readr::read_csv(file, col_names = T)
            #fread(file,header = T)
          return(df)
        }

readClipBoard <- function(sep="\t",colNames = F,rowNames = F) {
  rowNames= if (rowNames) {rowNames=1} else {rowNames=NULL}
  read.table(pipe("pbpaste"), sep="\t", header=colNames,row.names = rowNames)
  }

writeClipboard = function(x,sep="\t",col.names=T,...) {
  write.table(x
              ,file = pipe("pbcopy")
              ,sep=sep
              ,col.names = col.names
              ,row.names = F
              ,quote = F,...)
}


save_df_to_file <- function(df,fileName)
{
  script_dir=ifelse(exists("script_path_from_fun"),script_path_from_fun,fjComm::get_scriptpath())

  if (dirname(fileName)=="."){fileName=paste0(script_dir,"/",fileName)} # specify if not specified dir name
  if (!dir.exists(dirname(fileName))) dir.create(dirname(fileName),recursive = TRUE) # create target dir if not exist

  write.csv(df, fileName, quote = F, row.names = F)
}
