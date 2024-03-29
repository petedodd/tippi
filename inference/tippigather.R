## this file is to generate composite outputs
library(here)
library(data.table)
library(ggplot2)
library(scales)
library(ggpubr)
library(glue)
dg <- 1
ft <- function(x)format(x,digits=dg,nsmall=dg)
source(here('dataprep/tippifunctions.R'))
infpath <- here('inference')
gh <- function(x) glue(here(infpath,x))

## qty <- 'tx'; page <- 2
## shhs[page,age]

## make single graphs
grphsE <- grphs <- MCF <- CDX <- list()
for(page in 1:3){
  for(qty in c("tx","px")){
    cat(qty,"...\n")
    ## --- read data
    ## load data
    if(qty=='tx'){ #select qty
      D <- fread(gh('../dataprep/outdata/T.csv'))
    } else{
      D <- fread(gh('../dataprep/outdata/P.csv'))
    }
    D <- D[age==shhs[page,aged]] #select age
    D[,Facility:=site]
    D[,index:=1:nrow(D)] #for start/ends
    setkey(D,Country,Facility)

    ## convergence dx
    tmp <- scan(file=gh('outdata/cdx_{qty}_{shhs[page,age]}.txt'))
    tmp <- as.data.table(tmp)
    tmp[,variable:=c('median Rhat','min Rhat','max Rhat','median ESS','min ESS','max ESS')]
    tmp[,qty:=qty]; tmp[,age:=shhs[page,aged]]
    CDX[[glue('{qty}_{shhs[page,age]}')]] <- tmp

    ## empirical data
    ## sites
    siteeffects <- D[,.(site.effect = (Intervention.Num/
                                       Intervention.FT)/
                          (Baseline.Num/Baseline.FT),
                        int.number = Intervention.Num,
                        country=Country)]
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

    ## reorder & add circonflex
    lvl <- unique(as.character(gsub("ote","ôte",D$Country)))
    lvl <- rev(lvl)
    siteeffects$country <- factor(gsub("ote","ôte",siteeffects$country),levels=lvl)
    countryeffects$country <- factor(gsub("ote","ôte",countryeffects$country),levels=lvl)
    bsmy$country <- factor(gsub("ote","ôte",bsmy$country),levels=lvl)
    clz <- RColorBrewer::brewer.pal(length(lvl),'Paired')

    ## harmonization
    siteeffects[,c('RR.mid','RR.lo','RR.hi'):=NA]
    countryeffects[,c('RR.mid','RR.lo','RR.hi'):=NA]

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

    ## assemble compare methods data
    bsmy[,method:='Bayesian site<country model']
    countryeffects[,method:='empirical']
    tmp <- rbindlist(list(
      countryeffects[,.(country,RR=country.effect,
                        RR.lo,RR.hi,method)],
      bsmy[,.(country,RR=RR.mid,RR.lo,RR.hi,method)]
    ))
    tmp[,age:=shhs[page,aged]]
    tmp[,qty:=ifelse(qty=='px','TPT','ATT')]
    MCF[[glue("{qty}{shhs[page,age]}")]] <- tmp
  }
}

## convergence dx
CDX <- rbindlist(CDX)
CDX <- dcast(CDX,qty+age~variable,value.var = 'tmp')
nmz <- names(CDX)
nmz <- nmz[c(1,2,6,8,4,5,7,3)]
setcolorder(CDX,nmz)
fwrite(CDX,file=gh('outdata/CDX.csv'))


## make joined graphs
GA <- ggarrange(plotlist = grphs,
                ncol=2,nrow=3,
                labels = paste0(letters[1:6],")"),
                common.legend = TRUE,legend='top')
ggsave(filename=gh("graphs/MAll2.eps"),w=13,h=12)
ggsave(GA,filename=gh("graphs/MAll2.png"),w=13,h=12)

GA <- ggarrange(plotlist = grphsE,
                ncol=2,nrow=3,
                labels = paste0(letters[1:6],")"),
                common.legend = TRUE,legend='top')
ggsave(filename=gh("graphs/MAllE2.eps"),w=13,h=12)
ggsave(GA,filename=gh("graphs/MAllE2.png"),w=13,h=12)


## model comparison
MCF <- rbindlist(MCF)
MCF <- merge(MCF,
             MCF[method=='Bayesian site<country model',
                 .(ref=RR,qty,age,country)],
             by=c('country','age','qty'),
             all.x=TRUE)
MCF <- MCF[method %in%c('Bayesian site<country model','empirical')]


GP <- ggplot(MCF,aes(RR,RR/ref,
                     ymin=RR.lo/ref,ymax=RR.hi/ref,
                     col=country,shape=method))+
  geom_hline(yintercept = 1,col=2,alpha=0.5,lty=2)+
  geom_point()+
  geom_errorbar(width=0)+
  scale_x_log10()+  scale_y_log10()+
  facet_grid(age~qty,scales='free')+
  theme_bw()+
  xlab('Incidence rate ratio')+
  ylab('Ratio compared to reference')

ggsave(GP,filename=gh("graphs/MCF.png"),w=13,h=12)
ggsave(GP,filename=gh("graphs/MCF.pdf"),w=13,h=12)


MCF <- MCF[!grepl('meta',method)]
MCF2 <- dcast(MCF,country+age+qty~method,value.var = c('RR','RR.lo','RR.hi'))

GP <- ggplot(MCF2[age!='0-14'],
             aes(`RR_Bayesian site<country model`,
                 RR_empirical,
                 xmin=`RR.lo_Bayesian site<country model`,
                 xmax=`RR.hi_Bayesian site<country model`,
                 col=country))+
  geom_abline(intercept = 0,slope=1,col=2,alpha=0.5,lty=2)+
  geom_point()+
  geom_errorbarh(height=0)+
  scale_x_log10()+  scale_y_log10()+
  facet_wrap(age~qty,scales='free')+
  theme_bw()+
  xlab('Bayesian mixed model IRR')+
  ylab('Empirical IRR')
## GP

ggsave(GP,filename=gh("graphs/MCF2.png"),w=13,h=12)
ggsave(GP,filename=gh("graphs/MCF2.pdf"),w=13,h=12)
