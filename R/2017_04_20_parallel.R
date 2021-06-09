parallel <- function(list_or_vect= 1:10, fun= function(x){x}, workers=15, ...)
{
  pacman::p_load(BiocParallel)
  curr_workers= workers
    max_workers=parallel::detectCores()-2
    if (curr_workers>max_workers) curr_workers=max_workers

    # register(SerialParam())
    # register(MulticoreParam())
  parallelParam= MulticoreParam(workers = curr_workers)
  bplapply(list_or_vect, FUN = fun, BPPARAM = parallelParam,...)
}


parallel_cmd_gen<-function(cmd="samtools view -c {}",params=c("file1","file2","file3"),threads=8)
  glue::glue("parallel -j {threads} -k {cmd} ::: ") %>%paste0(paste0(params,collapse = " "))
