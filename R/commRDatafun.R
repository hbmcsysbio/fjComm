mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  # tf147 %<>%  mutate_cond(TF==TFname & use!=0, upload=1)  //var should be created beforehand
  condition <- eval(substitute(condition), .data, envir)
  condition[is.na(condition)] = FALSE
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
