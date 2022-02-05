## this file is to generate composite outputs
library(here)
library(data.table)
library(ggplot2)
library(ggpubr)

## load graph data:
grphs <- list()
load(here('data/txMA304.Rdata'))
grphs[['tx04']] <- MAP
load(here('data/pxMA304.Rdata'))
grphs[['px04']] <- MAP
load(here('data/txMA3514.Rdata'))
grphs[['tx514']] <- MAP
load(here('data/pxMA3514.Rdata'))
grphs[['px514']] <- MAP
load(here('data/txMA3014.Rdata'))
grphs[['tx014']] <- MAP
load(here('data/pxMA3014.Rdata'))
grphs[['px014']] <- MAP


## make joined graph
ggarrange(plotlist = grphs,
          ncol=2,nrow=3,
          labels = paste0(letters[1:6],")"),
          common.legend = TRUE,legend='top')

ggsave(filename=here("graphs/MAll.eps"),w=10,h=12)
ggsave(filename=here("graphs/MAll.png"),w=10,h=12)

## TODO
## additional titles on graphs
## log scale
## text on graphs


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
