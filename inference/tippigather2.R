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

## qty <- 'tx'; page <- 2
## shhs[page,age]

## make single graphs
grphsE <- grphs <- list()
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
    fwrite(siteeffects,
           file=gh('outdata/siteeffects_{qty}_{shhs[page,age]}.csv'))
    ## country
    countryeffects <- D[,.(country.effect =
                             (sum(Intervention.Num)/
                              sum(Intervention.FT))/
                             (sum(Baseline.Num)/sum(Baseline.FT)),
                           int.number=sum(Intervention.Num) ),
                        by=.(country=Country)]
    fwrite(countryeffects,
           file=gh('outdata/countryeffects_{qty}_{shhs[page,age]}.csv'))

    ## MA results
    load(gh('outdata/bsmy_{qty}_{shhs[page,age]}.Rdata'))
    load(gh('outdata/fsmy_{qty}_{shhs[page,age]}.Rdata'))

    ## reorder
    lvl <- unique(as.character(D$Country))
    lvl <- rev(lvl)
    siteeffects$country <- factor(siteeffects$country,levels=lvl)
    countryeffects$country <- factor(countryeffects$country,levels=lvl)
    fsmy$country <- factor(fsmy$country,levels=lvl)
    bsmy$country <- factor(bsmy$country,levels=lvl)
    clz <- RColorBrewer::brewer.pal(length(lvl),'Paired')

    ## harmonization
    siteeffects[,c('RR.mid','RR.lo','RR.hi'):=NA]
    countryeffects[,c('RR.mid','RR.lo','RR.hi'):=NA]
    fsmy[,c('RR.lo','RR.hi'):=NA]

    ## --- make graph
    ## title
    ttl <- glue(ifelse(qty=='tx',
                       'anti-tuberculosis treatment initiation',
                       'tuberculosis preventive therapy initiation')) +
      ", " + shhs[page,aged] + ' years'

    MAP <- 
      ggplot(data=bsmy,aes(x=country,
                           y=RR.mid,ymin=RR.lo,ymax=RR.hi,
                           col=country)) +
      geom_point(size=2) +
      geom_point(data=fsmy,size=2,col=2,shape=15) +
      geom_point(data=siteeffects[is.finite(site.effect)],
                 aes(x=country,y=site.effect,
                     size=int.number),
                 shape=1) +
      geom_point(data=countryeffects,
                 aes(x=country,y=country.effect,
                     size=int.number),
                 shape=4) +
      ## scale_color_manual(values = clz)+
      geom_hline(yintercept = 1,lty=2,col='darkgrey')+
      geom_point(size=2,col='black') +
      geom_errorbar(width=0,col='black') +
      scale_y_sqrt(limits=c(0,150)) + #NOTE one DRC site omitted
      ylab('Rate ratio (square root scale)')+
      xlab('Country')+
      coord_flip() +
      theme_classic() + ggpubr::grids()+ ggtitle(ttl)+
      labs(size='Number (intervention)')+
      guides(colour="none") + 
      theme(legend.position = 'bottom')
    ## MAP

    ## save
    grphs[[glue("{qty}{shhs[page,age]}")]] <- MAP

    txtd <- copy(bsmy)
    txtd[,txt:=paste0(ft(RR.mid)," (",ft(RR.lo)," - ",ft(RR.hi),")")]

    MAPE <- MAP + geom_text(data=txtd,
                            aes(x=country,y=0.1,label=txt), #90
                            size=3,col='black')+
      scale_y_log10(limits=c(5e-2,1.4e2),
                    label=comma) +
      ylab('Rate ratio (log scale)')

    ## save
    grphsE[[glue("{qty}{shhs[page,age]}")]] <- MAPE
  }
}

## make joined graphs
GA <- ggarrange(plotlist = grphs,
                ncol=2,nrow=3,
                labels = paste0(letters[1:6],")"),
                common.legend = TRUE,legend='top')
ggsave(filename=here("graphs/MAll2.eps"),w=13,h=12)
ggsave(GA,filename=here("graphs/MAll2.png"),w=13,h=12)

GA <- ggarrange(plotlist = grphsE,
                ncol=2,nrow=3,
                labels = paste0(letters[1:6],")"),
                common.legend = TRUE,legend='top')
ggsave(filename=here("graphs/MAllE2.eps"),w=13,h=12)
ggsave(GA,filename=here("graphs/MAllE2.png"),w=13,h=12)
