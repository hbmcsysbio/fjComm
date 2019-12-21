largeListToDf <- function(list)
{
  n <- length(list[[1]])
  structure(list, row.names = c(NA, -n), class = "data.frame")
}

letter2num <- function(x) { sapply(x, function(x) {utf8ToInt(x) - utf8ToInt("A") + 1L}) }
num2letter <- function(x) {  sapply(x, function(x) {intToUtf8 (x + utf8ToInt("A") - 1L)})  }
