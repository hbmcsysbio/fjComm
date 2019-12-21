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


