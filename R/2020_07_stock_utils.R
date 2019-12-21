# fjComm::clear_()


stock_indices_plot <-function()
{
  # pacman::p_load(tidyquant)
  # pacman::p_load(smooth)
  pacman::p_load(quantmod)
  startDate = '2018-01-01';
  names=qw("Shanghai  Shengzhen  Hengsheng  Nikkei225  S&P500 Nasdaq DAX STOXX50  GOLD")
  symbol=qw("000001.SS  399001.SZ  ^HSI  ^N225  ^GSPC  ^IXIC  ^GDAXI  ^STOXX50E GC=F")

  allData=vector("list",length(names))

  for (i in 1:length(allData)) {
    allData[[i]]=getSymbols(symbol[i],src="yahoo",warnings = F,auto.assign = F,from=startDate) %>%
      as.data.frame() %>% set_colnames(qw("open high low close volume adjusted")) %>%
      mutate(date=rownames(.) %>% as.Date(),middle=high/2+low/2,label=names[i],volume=ifelse(volume==0,NA,volume)) %>%
      mutate(normVol=scales::rescale(volume)*sd(middle,na.rm = T)+mean(middle,na.rm = T))#, ma30=middle %>% rollmean(30,fill = NA))
  }

  allData=do.call(rbind,allData)


  ggplot(allData)+ geom_line(aes(date,normVol),color="grey")+#geom_line(aes(date,ma30),color="purple")+
    geom_boxplot(aes(lower=pmin(open,close), middle=middle, upper=pmax(open,close), ymin=low,
                    ymax=high, x=date, group=date, color=I(ifelse(open>close,"green","red"))),size=0.5,stat='identity',alpha=.8,)+
              theme(axis.text.x = element_text(angle = 90))+scale_x_date(expand = c(0,0),date_breaks = "2 month", date_labels =  "%b %Y")+
              ylab(NULL)+xlab(NULL)+facet_wrap(.~label,ncol = 2,scales = "free_y")

  }



# getSymbols("GOOG",src="yahoo",from = '2018-01-01',to = '2019-01-01')
# barChart(Shanghai,subset='last 6 months')
# candleChart(GOOG,col = TRUE,theme ="white",subset='last 2 months')
# addADX()
# available source: yahoo  google  FRED  oanda

