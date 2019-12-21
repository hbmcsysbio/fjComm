# common settings for ggplot2
# source(paste0(RlibDir,""))
# 1point = 0.35mm

# p_table <- ggplot_gtable(ggplot_build(p))

gg_theme_Publication <- function(base_size=9, base_family="Helvetica") {
  pacman::p_load(grid)
  pacman::p_load(ggthemes)

  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = rel(11/9), hjust = 0.5),
            text = element_text(),
            # panel.background = element_rect(colour = NA),
            # plot.background = element_rect(colour = NA),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = rel(9/9)),
            axis.title.y = element_text( margin=unit(c(0.5,0.5,0.5,0.5), "mm") ),
            axis.title.x = element_text( margin=unit(c(0.5,0.5,0.5,0.5), "mm")),
            axis.text = element_text(size=rel(9/9)),
            axis.line = element_line(colour="black",size = 0.4),
            axis.ticks = element_line(size = 0.4),
            panel.grid.major = element_blank(), #element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            # legend.position = "bottom",
            # legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.background = element_rect(fill = "transparent",colour = NA),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic",size=rel(9/9)),
            legend.text = element_text(size=rel(8/9)),
            plot.margin = unit(c(3,3,1.5,1.5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(size=rel(8/9)),
            # strip.background = element_blank(),
            legend.key.height=unit(8,"pt")
    ))

}

gg_theme_Publication_diag= gg_theme_Publication()+ theme(axis.line = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())
gg_theme_bordered_diag= theme(panel.border = element_rect(size = 0.5,colour = "black"))

gg_scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

gg_scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

gg_axis_x_noExp <- function(continuous=T,...) { if (continuous) scale_x_continuous(expand =c(0,0),...) else scale_x_discrete(expand =c(0,0))}
gg_axis_y_noExp <- function(continuous=T,...) { if (continuous) scale_y_continuous(expand =c(0,0),...) else scale_y_discrete(expand =c(0,0)) }











p_basic=ggplot()+ scale_x_continuous(expand = c(0,0)) + scale_y_discrete(expand = c(0,0))+ ggplot2::theme_bw()

gg_theme_transparent <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))

trans_bk <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  # axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  # axis.line = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA)
)

# most useful with + scale_x_continuous(expand = c(0,0))+scale_y_discrete(expand = c(0,0))
trans_bk_xyaxis <- theme(
  # axis.title.x = element_blank(),
  # axis.title.y = element_blank(),
  # axis.text.x = element_blank(),
  # axis.text.y = element_blank(),
  # axis.ticks = element_blank(),
  panel.grid = element_blank(),
  # axis.line = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA),
  legend.background = element_blank()
)

# dummy plot for positioning
createDummy<- function(x=100,y=100)
{
  df <- data.frame(a = c(1,x), b = c(1,y));
  base <- ggplot(df, aes(a, b)) +  geom_blank() +  gg_theme_transparent+ scale_x_continuous(expand = c(0,0))+ scale_y_continuous(expand = c(0,0))
}

no_expand_lab <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
)

ggbasic= ggplot()+ theme_bw() + scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0))
ggbasic_dis= ggplot()+ theme_bw() + scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c(0,0))


