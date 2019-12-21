# macro to delete all var and reload fjComm
clear_ <- gtools::defmacro(var,expr={rm(list=ls(all=TRUE));
  if(!is.na(fjComm::get_scriptpath())) setwd(fjComm::get_scriptpath());
  if("package:fjComm" %in% search()) detach("package:fjComm", character.only = TRUE, unload = TRUE);
  library(fjComm)
  })


detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

get_scriptpath <- function() {
  # location of script can depend on how it was invoked:
  # source() and knit() put it in sys.calls()
  path <- NULL

  if(!is.null(sys.calls())) {
    # get name of script - hope this is consisitent!
    path <- as.character(sys.call(1))[2]
    # make sure we got a file that ends in .R, .Rmd or .Rnw
    if (grepl("..+\\.[R|Rmd|Rnw]", path, perl=TRUE, ignore.case = TRUE) )  {
      script_path_from_fun <<- dirname(path) #!!!!!!!!!
      return(dirname(path))
    } else {
      message("Obtained value for path does not end with .R, .Rmd or .Rnw: ", path)
    }
  } else{
    # Rscript and R -f put it in commandArgs
    args <- commandArgs(trailingOnly = FALSE)
  }

  script_path_from_fun <<- dirname(path) #!!!!!!!
  return(dirname(path))
}

# script.dir=getSrcDirectory(function(x) {x})

stopQuietly <- function(...) {
  print(">>>>>> no error stop ©©")
  blankMsg <- sprintf("\r%s", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} # stopQuietly()

qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))



namedVect_to_df <- function(namedVect, name_label="name", value_label="value" ) namedVect %>% as.data.frame() %>% tibble::rownames_to_column("temp_") %>% set_colnames(c(name_label, value_label))



