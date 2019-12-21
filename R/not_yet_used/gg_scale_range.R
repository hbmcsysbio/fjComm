library(scales)
library(ggthemes)

range_breaks <- function(expand=c(0.04, 0), ..., threshold=0.5){
  function(x) {
    spread <- range(x)
    m <- matrix(c(1+expand[1], expand[1], expand[1], 1+expand[1]),
                ncol=2, byrow=TRUE)
    limits <- m %*% (spread + c(1,-1)*expand[2]) / (1+2*expand[1])

    pretty_breaks <-  scales::extended_breaks(...)(x)
    #     pretty_breaks <- pretty(x, ...)

    spacing <- min(diff(pretty_breaks))

    ## don't include pretty breaks if outside range
    keep <- pretty_breaks > limits[1] & pretty_breaks < limits[2]
    pretty_breaks <- pretty_breaks[keep]

    ## prevent potential label overlap at the edges
    n <- length(pretty_breaks)
    clash <- c(abs(limits[1] - pretty_breaks[1]) < threshold * spacing,
               abs(limits[2] - pretty_breaks[n]) < threshold * spacing)
    remove <- -c(1,n)[clash]

    all_breaks <- if(any(clash))
      c(limits, pretty_breaks[remove]) else
        c(limits, pretty_breaks)

    sort(all_breaks)
  }
}


scale_x_tufte <-  function(breaks = range_breaks(expand), ...,
                           labels = as.character,
                           expand=c(0.04, 0))
  continuous_scale(c("x", "xmin", "xmax", "xend", "xintercept"),
                   "position_c", identity,
                   breaks = breaks, ..., labels=labels,
                   expand = expand, guide = "none")

scale_y_tufte <-  function(breaks = range_breaks(expand), ...,
                           labels = as.character,
                           expand=c(0.04, 0))
  continuous_scale(c("y", "ymin", "ymax", "yend", "yintercept"),
                   "position_c", identity,
                   breaks = breaks, ..., labels=labels,
                   expand = expand, guide = "none")

tufte_x_axis <- function(...)
  geom_segment(x=-Inf, xend=-Inf, y=-Inf, yend=Inf, ..., .inherit.aes=FALSE)
tufte_y_axis <- function(...)
  geom_segment(x=-Inf, xend=-Inf, y=-Inf, yend=Inf, ..., .inherit.aes=FALSE)


range_breaks <- function(expand=c(0.04, 0), ..., threshold=0.5, digits = 1, FUN = scales::extended_breaks){
  function(x) {
    spread <- range(x)
    m <- matrix(c(1+expand[1], expand[1], expand[1], 1+expand[1]),
                ncol=2, byrow=TRUE)
    limits <- m %*% (spread + c(1,-1)*expand[2]) / (1+2*expand[1])
    limits <- round(limits, digits)

    pretty_breaks <-  FUN(...)(x)

    spacing <- min(diff(pretty_breaks))

    ## don't include pretty breaks if outside range
    keep <- pretty_breaks > limits[1] & pretty_breaks < limits[2]
    pretty_breaks <- pretty_breaks[keep]

    ## prevent potential label overlap at the edges
    n <- length(pretty_breaks)
    clash <- c(abs(limits[1] - pretty_breaks[1]) < threshold * spacing,
               abs(limits[2] - pretty_breaks[n]) < threshold * spacing)
    remove <- -c(1,n)[clash]

    all_breaks <- if(any(clash))
      c(limits, pretty_breaks[remove]) else
        c(limits, pretty_breaks)


    sort(all_breaks)
  }
}
