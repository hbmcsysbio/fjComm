gg_distMDS_2Dplot <- function(distMat, labels=NULL, color="black", fig_anno="", textsize=3, textalpha=0.5, pointsize=2,pointalpha=1)
{
  # visualize distance matrix with 2D xyplot
  # if(is.null(labels))labels=rownames(distMat)
  p_load(ggrepel)
  fit=cmdscale(distMat,eig=TRUE, k=2)
  plotdf=data.frame(x=fit$points[,1], y=fit$points[,2])#,color=color,labels=labels)
  p=ggplot(plotdf,aes(x,y))+geom_point(aes(color=I(color)),size=pointsize,alpha=pointalpha) +xlab("")+ylab("")+gg_anno_grob(x = 0.1,fig_anno)+ coord_fixed(ratio = 1)
  if(!is.null(labels)) p=p+geom_text_repel(aes(label=labels,color=I(color)),size=textsize,alpha=textalpha)
  p
}


gg_PCA_2Dplot <- function(Mat, labels=NULL, color="black", fig_anno="", textsize=3, textalpha=0.5, pointsize=2,pointalpha=1)
{
  # visualize PCA dim1,2 with 2D xyplot
  # if(is.null(labels))labels=rownames(distMat)
  p_load(ggrepel,FactoMineR)
  pca = PCA(Mat, graph = F,ncp = 5)
  pcaplot2=pca$ind$coord %>% as.data.frame() #%>%  mutate(text=labels,color=color)
  p=ggplot(pcaplot2)+geom_point(aes(x=Dim.1, y=Dim.2,color=I(color)),size=pointsize,alpha=pointalpha) +xlab("")+ylab("")+gg_anno_grob(x = 0.1,fig_anno)+ coord_fixed(ratio = 1)
  if(!is.null(labels)) p=p+geom_text_repel(aes(x=Dim.1,y=Dim.2,label=labels,color=I(color)),size=textsize,alpha=textalpha)
  p
}

gg_hist <- function(Arr,binwidth=1,fill="white",color="grey30",...)
{
  force(Arr)
  ggplot()+geom_histogram(aes(Arr),binwidth = binwidth,fill=fill,color=color,...)
}

gg_histly <- function(Arr,binwidth=1,fill="white",color="grey30",...)
{
  pacman::p_load(plotly)
  (ggplot()+geom_histogram(aes(Arr),binwidth = binwidth,fill=fill,color=color,...)) %>% ggplotly()
}

gg_2Dheat_dotly <-function(df,scale_dot_size=TRUE)
  # visualize 2D MI heatmap top kmer combs
{
  dot_sizes= if(scale_dot_size) {df$topMIsum %>% {(.-min(.))/max(.)*1.2} } else{1}
  (ggplot(df)+geom_point(aes(pos1,pos2,color=topk1,text=paste0(topk1,pos2-pos1-3,topk2),size=I(dot_sizes) ))) %>% ggplotly()
}


gg_multi_ggoutput<-function(plotlist,ncol=1,col_to_row=FALSE,magnify=100)
{
  plotNum=length(plotlist)
  cols=ncol
  rows=(plotNum/cols) %>% base::ceiling()
  plotNum_=cols*rows

  options(warn=-1)
  layoutMat=matrix(1:plotNum_,ncol=ncol,byrow = TRUE);
  if(col_to_row) {layoutMat %<>% t; cols=rows; rows=ncol}
  plotlist %<>% purrr::map(~ggplotGrob(.))

  placeholder= createDummy(5,5) %>% ggplotGrob()
  plotlist= c(plotlist, vector("list", plotNum_-plotNum) %>% purrr::map(~placeholder))
  curr_p= createDummy((cols)*magnify,(rows)*magnify)+theme(plot.margin = margin())

  for(row in 1:nrow(layoutMat) )
  {
    for (col in 1:ncol(layoutMat))
    {
      curr_p= curr_p+annotation_custom(grob = plotlist[[layoutMat[row,col]]], xmin = (col-1)*magnify, xmax = (col)*magnify, ymin =(rows-row)*magnify, ymax = (rows-row+1)*magnify)
      # browser()
    }
  }
  options(warn=0)
  curr_p
}
