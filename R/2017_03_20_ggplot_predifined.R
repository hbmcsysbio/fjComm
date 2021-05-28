


# colors
  gg_heat_rainbow= c("black", "purple4", "#3d33ff", "#03efd7", "green", "yellow", "red")
  gg_heat_rainbow_plus= c(gg_heat_rainbow,"deeppink","plum1","thistle1")
  gg_steelblue_red= c("white","steelblue","red")
  gg_green_red=c("#ceffce","red")
  gg_bkcolor <- function(bkcolor="white") {theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_rect(fill = bkcolor))}


  gg_to_polar <- function(gg=last_plot(),limits=c(0,10),breaks = c(0,0.1,limits[2]/2))
  {
    gg+theme(axis.line = element_blank(), panel.grid.major.y =element_line(size = 0.4), panel.grid.major.x = element_blank(),legend.position = "bottom") +
      coord_polar() + scale_y_continuous( limits=limits,breaks = breaks)
  }

# predifined plots
  gg_heat2D_MI <- function(dataFrame="col1,2 pos, col3 value", grad_colors=gg_heat_rainbow, bkcolor="navy",limits=NULL) # 1st and 2st col are pos info, 3rd col is the value
  {
    names_col=colnames(dataFrame)

    legend_col="white"; if (bkcolor=="white") legend_col="black"
    if(is.null(limits)) {limits= dataFrame[[3]] %>% quantile(c(0.005,0.995))} # c(0.2,2) for maxBias

    meanDiag15= dataFrame %>% dplyr::filter(abs(dataFrame[[2]]-dataFrame[[1]])<=15) %>% .[[3]] %>% mean %>% prettyNum(digits=4)
    anno_grob= grid::grobTree(grid::textGrob(paste0("mean_diag (<=15)\n",meanDiag15),x=0.25, y=0.1, hjust=0, gp= grid::gpar(col=legend_col,fontsize=15,fontface="italic")))
    disp_text=with(dataFrame, ifelse(is.na(get(names_col[5])),character(0),paste(get(names_col[1]),get(names_col[2]),get(names_col[5]),get(names_col[6]),sep = "_") ))

    p=ggplot()+ theme_bw() + geom_raster(data = dataFrame, aes(get(names_col[1]),get(names_col[2]),fill = get(names_col[3]),text=disp_text  )) +labs(x=names_col[1],y=names_col[2])+
      scale_fill_gradientn(colours= grad_colors , na.value="grey",limits=limits, name = names_col[3], oob=scales::squish )+
      theme(legend.text = element_text(colour=legend_col), legend.title = element_text(colour=legend_col))+
      # guides(fill = guide_legend(title.theme = element_text(size=15, face="italic", colour = "red", angle = 45)))+
      scale_x_continuous(expand = c(0, 0),position = "top", breaks =  gg_breaks(dataFrame[[1]],end_plus_middle=T))+scale_y_continuous(expand = c(0, 0), breaks =  gg_breaks(dataFrame[[2]],end_plus_middle=T)) +
      theme(axis.text=element_text(size=6),axis.title=element_text(size=13,face="bold"), axis.text.x = element_text(angle = 90, hjust = 1))+ xlab("(bp)")+ylab("(bp)")+
      annotation_custom(anno_grob)+
      theme(legend.position = c(.85, .30),legend.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_rect(fill = bkcolor) )
    return(p+gg_theme_Publication())
  }

  gg_heat2D_MI_rotate<- function(dataFrame="col1,2 pos, col3 value", grad_colors=gg_heat_rainbow, bkcolor="navy",limits=NULL) # 1st and 2st col are pos info, 3rd col is the value
  {
    # browser()
    names_col=colnames(dataFrame)
    legend_col="white"; if (bkcolor=="white") legend_col="black"
    if(is.null(limits)) {limits= dataFrame[[3]] %>% quantile(c(0.005,0.995))} # c(0.2,2) for maxBias
    dataFrame %<>% mutate(pos1_=pos1,pos2_=pos2) %>% mutate(pos1=(pos1_+2+pos2_)/2, pos2=(pos2_-pos1_-3)/2) %>% mutate(id=paste0(pos1_,pos2_))
    dataFrame %<>% group_by(id) %>% nest() %>% mutate(data=map(data,function(x){ tmp=x[c(1,1,1,1),]; tmp$pos1=x$pos1 %>% {c(.,.+.5,.,.-.5)}; tmp$pos2= x$pos2 %>% {c(.+.5,.,.-.5,.)}; tmp })) %>% unnest(cols = c(data))#unnest(.drop = F)
    p=ggplot()+ theme_bw() + geom_polygon(data = dataFrame, aes(get(names_col[1]),get(names_col[2]),fill = get(names_col[3]),group=id )) +
      labs(x=names_col[1],y=names_col[2])+
      scale_fill_gradientn(colours= grad_colors , na.value="grey",limits=limits, name = names_col[3], oob=scales::squish )+
      scale_x_continuous(expand = c(0, 0))+scale_y_continuous(expand = c(0, 0)) +
      xlab("(bp)")+ylab("(bp)")+gg_theme_Publication()+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.line.y = element_blank())

    return(p)
  }

  gg_heat2D_diag <- function(dataFrame="col1,2 pos, col3 value", grad_colors=gg_heat_rainbow, bkcolor="navy",limits=NULL) # 1st and 2st col are pos info, 3rd col is the value
  {
    names_col=colnames(dataFrame)

    diagSpacing= with(dataFrame, abs(pos2-pos1)) %>% min
    dataFrame= dataFrame %>% dplyr::filter(abs(pos2-pos1)==diagSpacing)
    if(is.null(limits)) {limits= dataFrame[[3]] %>% quantile(c(0.005,0.995))} # c(0.2,2) for maxBias
    p=ggplot()+ theme_bw() + geom_raster(data = dataFrame, aes(get(names_col[1]),1,fill = get(names_col[3]))) +labs(x=names_col[1])+
      scale_fill_gradientn(colours= grad_colors , na.value="grey",limits=limits, name = names_col[3], oob=scales::squish )+
      scale_x_continuous(expand = c(0, 0),position = "bottom", breaks =  gg_breaks(dataFrame[[1]],end_plus_middle=T))+ xlab("(bp)")+ guides(fill=FALSE)+

      scale_y_continuous(expand = c(0, 0)) +
      theme(axis.text=element_text(size=6),axis.title=element_text(size=13,face="bold"), axis.text.x = element_text(angle = 0, hjust = 0.5), axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.title.y = element_blank())+
      theme(legend.position = c(.85, .30),legend.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_rect(fill = bkcolor) )

    return(p+gg_theme_Publication()+theme(axis.line = element_blank(),axis.ticks.y = element_blank(), axis.text.y = element_blank(),axis.title.y = element_blank()))
  }

  # plot heat diag with vector
  gg_heat2D_diag_vect <- function(Vect,grad_colors=gg_heat_rainbow, bkcolor="navy",limits=NULL)
    gg_heat2D_diag(data_frame(pos1=1:length(Vect),pos2=1:length(Vect),signal=Vect),grad_colors=grad_colors, bkcolor=bkcolor,limits=limits)

  # violin plot with vectors
  gg_vioplot <- function(...)
  {
    # parse dot and get the var name
    dots=list(...)
    dotNames= stringr::str_match(paste0(deparse(substitute(list(...))),collapse = ""),"list\\((.*)\\)")[2] %>% stringr::str_split(",\\s*")
    dotNames= dotNames[[1]]

    # construct df for plot
    repEach= lapply(dots,length)
    totalLen= repEach%>% unlist %>% sum
    cnt=0
    elementNames=lapply(repEach,function(x){cnt<<- cnt+1;rep(dotNames[cnt],x)})
    df= data.frame(unlist(elementNames),unlist(dots))
    colnames(df)=c("V1","V2")
    df= df %>% mutate(V1=as.factor(V1), V2=as.numeric(V2))

    # # plot
    pacman::p_load(ggpubr)
    ggviolin(df, x = "V1", y = "V2", fill = "V1",
             # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"))+ xlab("")+ ylab("value") + scale_fill_discrete(guide=FALSE)
  }

# misc: additions use "+" in gg
  #>> cust
  gg_axis_x_labels<- function(x_labels=c("tick1","tick2","tick3"),...){scale_x_discrete(labels=x_labels,...)}
  gg_axis_x_label_angle<- function(angle=45,...){theme(axis.text.x = element_text(angle = angle, hjust = 1),...)}

  gg_anno_grob <- function(text,col="black",x=0.25, y=0.9, fontsize=9,fontface="italic")
  {
    pacman::p_load(grid)
    anno_grob= grid::grobTree(grid::textGrob(text,x=x, y=y, hjust=0, gp= grid::gpar(col=col,fontsize=fontsize,fontface=fontface)))
    return(annotation_custom(anno_grob))
  }
  gg_add_corr <- function(var1,var2,method="pearson",x=0.25, y=0.9,fontsize=9, fontface= "italic", col="black")
  {
    corr=cor(var1,var2,method =method )
    gg_anno_grob(text = paste0(method," corr: ",corr %>% prettyNum(digit=2)),col = col, x = x,y = y,fontsize = fontsize,fontface = fontface)
  }



# gen breaks
  axis_breaks<- function(numVec, fromMin=T, interval=5, alt_from=0 )
  {
    minVal= alt_from
    if(fromMin) minVal=min(numVec) %>% ceiling
    maxVal= max(numVec) %>% floor
    ticks=(minVal:max(numVec)) [(minVal:max(numVec) %% interval)==0]
    if(length(ticks)<3) ticks=axis_breaks(numVec*10,interval = interval, alt_from = alt_from, fromMin = fromMin)/10
    return(ticks)
  } # only display ticks end with 5 or 0

  adjustFactor <-function(num)
  {
    adjFactor=10^(-floor(log10(num-0.00000000000001)));# avoid error when exactly 10
    # adj=num*adjFactor
    return(adjFactor)
  }

  # if limits is provided, numVect is not useful
  gg_breaks <-function(numVect,end_plus_middle=T,fromMin=T, alt_from=0, interval=10, limits=NULL,two_ticks=T)
  {
    if(!is.null(limits))numVect=limits
    tick_max=max(numVect)
    tick_min=min(numVect)
    adjFactor= adjustFactor(tick_max-tick_min)
    if (end_plus_middle)
    {
      # if((tick_max-tick_max %>% floor)<=0.1*(tick_max-tick_min))tick_max=
      decimalplaces <- function(x) {
        if ((x %% 1) != 0) {
          nchar(strsplit(sub('0+$', '', format(x,scientific = F)), ".", fixed=TRUE)[[1]][[2]])
        } else {
          return(0)
        }
      }
      tick_max=tick_max
      tick_min=tick_min
      tick_middle=(tick_max+tick_min)/2
      max_digits=max(decimalplaces(tick_min),decimalplaces(tick_max))
      if (decimalplaces(tick_middle)>max_digits) tick_middle=floor(tick_middle*10^max_digits)/10^max_digits
      ticks_3= c(tick_min,tick_middle,tick_max)
      if(two_ticks)ticks_3=c(tick_min %>% base::ceiling(),tick_max %>% base::floor())
    }
    else{
      numVect=numVect*adjFactor
      ticks=axis_breaks(numVect,interval = interval, alt_from = alt_from, fromMin = fromMin)/adjFactor
      tick_min=min(ticks)
      tick_max=ifelse(((length(ticks)-3) %% 2==0), max(ticks), ticks[length(ticks)-1])
      ticks_3= c(tick_min,mean(c(tick_min,tick_max)),tick_max)
      if(two_ticks)ticks_3=ticks_3[-2]
    }
    return(ticks_3)
  }
  scale_x_continuous1=purrr::partial(scale_x_continuous,expand=c(0,0),breaks=gg_breaks)
  scale_y_continuous1=purrr::partial(scale_y_continuous,expand=c(0,0),breaks=gg_breaks)


  # generate scale with 0 in the middle
  gg_gen_scale_with_quan <- function(vect, quan1, quan2)
  {
    quans=quantile(vect,c(quan1, quan2))
    (c(quans[1],0,quans[2])-quans[1])/(quans[2]-quans[1])
  }

  gg_multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
      print(plots[[1]])

    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }


  # path= ifelse(exists("script_path_from_fun"),script_path_from_fun,fjComm::get_scriptpath())
# save related
  gg_save_all<- function(plot, width=4,height=4, path= ".", newName=NULL, save_plotly=FALSE, ... )
  {
    width=width/2.54; height=height/2.54 # inch to cm
    if (is.null(newName)) filename= deparse(substitute(plot)) else filename= newName
    # if (path!=fjComm::get_scriptpath()){path=paste0(fjComm::get_scriptpath(),"/",path)}
    fullpath= paste0(path,"/",filename)
    if (!dir.exists(dirname(fullpath))) dir.create(dirname(fullpath),recursive = T)
    gg_save_pdf(plot = plot,width=width, height=height, path = path, filename=filename, pt=TRUE, ... );
# browser()
    if(system("uname",intern = T)=="Darwin") {gg_save_png(plot = plot, width=width,height=height, path = path, filename=filename, pt=TRUE, ... )}
    else {system( paste0("inkscape -f ",fullpath,".pdf -b white -d 300 --export-png=",fullpath,".png") )}
    # print("path: " %>% paste0(path))
    # print("filename: " %>% paste0(filename))
    # print("filename: " %>% paste0(filename))
    if(save_plotly) {gg_save_plotly(plot = plot, path=path, filename=filename)}
  }

  gg_save_pdf<- function(plot, width=4,height=4, path= ".", filename=NULL, pt=FALSE, ... )
  {
    if(!pt){width=width/2.54; height=height/2.54} # inch to cm
    filename = ifelse(is.null(filename),paste0(deparse(substitute(plot)),".pdf"),paste0(filename,".pdf") )
    # if (path!=fjComm::get_scriptpath()){path=paste0(fjComm::get_scriptpath(),"/",path)}
    fullpath= paste0(path,"/",filename)
    if (!dir.exists(dirname(fullpath))) dir.create(dirname(fullpath),recursive = T)
    ggsave(plot = plot,width=width, height=height, path = path,filename = filename, ... )
  }

  gg_save_png<- function(plot, width=4,height=4, path= ".", filename=NULL, pt=FALSE, ... )
  {
    if(!pt){width=width/2.54; height=height/2.54} # inch to cm
    filename = ifelse(is.null(filename),paste0(deparse(substitute(plot)),".png"),paste0(filename,".png"))
    # if (path!=fjComm::get_scriptpath()){path=paste0(fjComm::get_scriptpath(),"/",path)}
    fullpath= paste0(path,"/",filename)
    if (!dir.exists(dirname(fullpath))) dir.create(dirname(fullpath),recursive = T)
    # if(system("uname",intern = T)!="Darwin") {gg_save_pdf(plot = plot,width=width, height=height, path = path, filename=filename, pt=TRUE, ... )};
    ggsave(plot = plot,width=width, height=height, path = path,filename = filename, ... )
  }

  gg_save_plotly <- function(plot, path=".", filename=NULL)
  {
    # filename= deparse(substitute(plot))
    filename= ifelse(is.null(filename),paste0(deparse(substitute(plot)),".html"),paste0(filename,".html") )
    # if (path!=fjComm::get_scriptpath()){path=paste0(fjComm::get_scriptpath(),"/",path)}
    filename= paste0(path,"/",filename)
    if (!dir.exists(dirname(filename))) dir.create(dirname(filename),recursive = T)
    plot=plotly::ggplotly(plot)
    htmlwidgets::saveWidget(plot, file = paste0(getwd(),"/",filename))
  }
  gg_save_diag<- function(plot, width = 4.25,height = 1.8, path= ".", newName=NULL, ... )
  {
    width=width/2.54; height=height/2.54 # inch to cm
    if (is.null(newName)) filename= deparse(substitute(plot)) else filename= newName
    gg_save_pdf(plot = plot,width=width, height=height, path = path, filename=filename, pt=TRUE, ... );
    gg_save_png(plot = plot, width=width,height=height, path = path, filename=filename, pt=TRUE, ... )
  }





#
  ggheat=function(m, rescaling='none', clustering='none', labCol=T, labRow=T, border=FALSE,
                  heatscale= c(low='blue',high='red'))
  {
    ## m=matrix(data=sample(rnorm(100,mean=0,sd=2)), ncol=10)
    ## this function makes a graphically appealing heatmap (no dendrogram) using ggplot
    ## whilst it contains fewer options than gplots::heatmap.2 I prefer its style and flexibility

    ## the function can be be viewed as a two step process
    ## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
    ## using simple options or by a user supplied function
    ## 2. with the now resahped data the plot, the chosen labels and plot style are built

    require(reshape2)
    require(ggplot2)

    ## you can either scale by row or column not both!
    ## if you wish to scale by both or use a differen scale method then simply supply a scale
    ## function instead NB scale is a base funct

    if(is.function(rescaling))
    {
      m=rescaling(m)
    }
    else
    {
      if(rescaling=='column')
        m=scale(m, center=T)
      if(rescaling=='row')
        m=t(scale(t(m),center=T))
    }

    ## I have supplied the default cluster and euclidean distance- and chose to cluster after scaling
    ## if you want a different distance/cluster method-- or to cluster and then scale
    ## then you can supply a custom function

    if(is.function(clustering))
    {
      m=clustering(m)
    }else
    {
      if(clustering=='row')
        m=m[hclust(dist(m))$order, ]
      if(clustering=='column')
        m=m[,hclust(dist(t(m)))$order]
      if(clustering=='both')
        m=m[hclust(dist(m))$order ,hclust(dist(t(m)))$order]
    }
    ## this is just reshaping into a ggplot format matrix and making a ggplot layer

    rows=dim(m)[1]
    cols=dim(m)[2]
    melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows) ,melt(m))
    melt.m %<>% mutate(Var2=colInd)
    g=ggplot(data=melt.m)

    ## add the heat tiles with or without a white border for clarity
    if(border==TRUE)
      g2=g+geom_rect(aes(xmin=Var2-0.5,xmax=Var2+0.5,ymin=rowInd-1,ymax=rowInd, fill=value),colour='white')
    if(border==FALSE)
      g2=g+geom_rect(aes(xmin=Var2-0.5,xmax=Var2+0.5,ymin=rowInd-1,ymax=rowInd, fill=value))

    ## add axis labels either supplied or from the colnames rownames of the matrix

    # if(labCol==T)
    #   g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
    # if(labCol==F)
    #   g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))

    if(labRow==T)
      g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m),expand = c(0,0))
    if(labRow==F)
      g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows),expand = c(0,0))

    ## get rid of grey panel background and gridlines

    g2=g2+theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA),
                panel.background=element_rect(fill=NA, colour=NA))

    ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
    return(g2+scale_fill_continuous("", heatscale[1], heatscale[2]))

  }


  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }


