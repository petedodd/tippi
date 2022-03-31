## this file is to generate composite outputs
library(here)
library(data.table)
library(ggplot2)
library(scales)
library(ggpubr)
library(glue)
gh <- function(x) glue(here(x))
dg <- 1
ft <- function(x)format(x,digits=dg,nsmall=dg)
source(here('../dataprep/tippifunctions.R'))

## make single graphs
grphs <- list()
for(page in 1:3){
  for(qty in c("tx","px")){
    cat(qty,"...\n")
    ## --- read data
    ## load data
    if(qty=='tx'){ #select qty
      D <- fread(here('../dataprep/outdata/T.csv'))
    } else{
      D <- fread(here('../dataprep/outdata/P.csv'))
    }
    D <- D[age==shhs[page,aged]] #select age
    D[,Facility:=site]
    D[,index:=1:nrow(D)] #for start/ends
    setkey(D,Country,Facility)

    ## empirical data
    ## sites
    siteeffects <- D[,.(site.effect = (Intervention.Num/
                                       Intervention.FT)/
                          (Baseline.Num/Baseline.FT),
                        int.number = Intervention.Num,
                        country=Country,
                        index)]
    fwrite(siteeffects,file=gh('outdata/siteeffects_{qty}_{shhs[page,age]}.csv'))
    ## country
    countryeffects <- D[,.(country.effect =
                             (sum(Intervention.Num)/
                              sum(Intervention.FT))/
                             (sum(Baseline.Num)/sum(Baseline.FT)),
                           int.number=sum(Intervention.Num) ),
                        by=.(country=Country)]
    fwrite(countryeffects,file=gh('outdata/countryeffects_{qty}_{shhs[page,age]}.csv'))

    ## MA results
    MC <- fread(gh('outdata/{qty}MC{shhs[page,age]}.csv'))
    MC$country <- factor(MC$country,
                         levels=rev(MC$country),ordered = TRUE)
    MC[,txt:=paste0(ft(`50%`)," (",ft(`2.5%`)," - ",ft(`97.5%`),")")]

    ## --- make graph
    ## title
    ttl <- glue(ifelse(qty=='tx',
                       'anti-tuberculosis treatment initiation',
                       'tuberculosis preventive therapy initiation')) +
      ", " + shhs[page,aged] + ' years'
    ## graph
    MAP <-
      ggplot(MC,aes(country,`50%`)) +
      geom_point(size=2) +
      geom_point(data=siteeffects[is.finite(site.effect)],
                 aes(country,site.effect,
                     size=int.number,col=country),
                 shape=1) +
      geom_point(data=countryeffects,aes(country,country.effect,
                                         size=int.number,col=country),
                 shape=4) +
      geom_hline(yintercept = 1,lty=2,col='darkgrey')+
      geom_point(size=2) +
      geom_errorbar(aes(ymin=`2.5%`,ymax=`97.5%`),width=0) +
      geom_text(aes(x=country,y=90,label=txt),size=3)+
      scale_y_log10(limits=c(5e-2,1.4e2),
                    label=comma) + #NOTE one DRC site omitted
      ylab('Rate ratio (log scale)')+
      xlab('Country')+
      coord_flip() +
      theme_classic() + ggpubr::grids()+ ggtitle(ttl)+
      labs(size='Number (intervention)')+
      guides(colour="none") + 
      theme(legend.position = 'bottom')

    ## save
    grphs[[glue("{qty}{shhs[page,age]}")]] <- MAP
  }
}
## make joined graph
GA <- ggarrange(plotlist = grphs,
                ncol=2,nrow=3,
                labels = paste0(letters[1:6],")"),
                common.legend = TRUE,legend='top')
ggsave(filename=here("graphs/MAll.eps"),w=13,h=12)
ggsave(GA,filename=here("graphs/MAll.png"),w=13,h=12)


## gather summaries for reporting
## tx
T1 <- fread(here('outdata/txMC04.csv'))
T2 <- fread(here('outdata/txMC514.csv'))
T3 <- fread(here('outdata/txMC014.csv'))

TALL <- rbindlist(list(
    T1[country=='SUMMARY',.(`50%`,`2.5%`,`97.5%`)],
    T2[country=='SUMMARY',.(`50%`,`2.5%`,`97.5%`)],
    T3[country=='SUMMARY',.(`50%`,`2.5%`,`97.5%`)]
))

fwrite(format(TALL,digits=3,nsmall=3),file=here('outdata/txALL.csv'))

## px
P1 <- fread(here('outdata/pxMC04.csv'))
P2 <- fread(here('outdata/pxMC514.csv'))
P3 <- fread(here('outdata/pxMC014.csv'))

PALL <- rbindlist(list(
    P1[country=='SUMMARY',.(`50%`,`2.5%`,`97.5%`)],
    P2[country=='SUMMARY',.(`50%`,`2.5%`,`97.5%`)],
    P3[country=='SUMMARY',.(`50%`,`2.5%`,`97.5%`)]
))

fwrite(format(PALL,digits=3,nsmall=3),file=here('outdata/pxALL.csv'))
