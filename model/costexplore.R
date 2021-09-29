## trying to understand incremental costs etc
library(ggplot2)
library(data.table)

library(here)
load(here("data/T.Rdata"))
load(here("data/PT.Rdata"))

## ATT


## PT
## costACF.soc
## cost.soc, cost.int


## PT[,costPT.soc:= (hhc/ptentry)*traceperhhcpt*socu_a.thct +
##         (1-hhc/ptentry)*socu_a.hct + socu_a.tpt]
## PT[,costPT.int:= (hhc/ptentry)*traceperhhcpt*intu_a.thct +
##         (1-hhc/ptentry)*intu_a.hct + intu_a.tpt]


## ## cost (not including thct)
## PT[,costACF.soc.hh:= ## coprev HHCM'd
##         totindexscreen*(Screened*socu_a.hct +
##                         Presumptive*socu_a.tbe +
##                         Diagnosed*socu_a.att)]
## PT[,costACF.soc.nhh:=
##         ## coprev not HHCM'd - not including PHC screening cost
##         (totindexscreen*Diagnosed*cdr*(socu_a.att+socu_a.tbe))]
## PT[,costACF.int.hh:= ## coprev HHCM'd
##         totindexscreen*(Screened*intu_a.hct +
##                         Presumptive*intu_a.tbe +
##                         Diagnosed*intu_a.att)]
## PT[,c('costACF.soc.hh','costACF.soc.nhh','costACF.int'):=
##         .(sum(costACF.soc.hh),
##           sum(costACF.soc.nhh),
##           sum(costACF.int.hh)),
##    by=.(id,country)] #sum over ages - per PT



PTS <- PT[,.(RR=mean(RR), #intervention effect
             ## total costs
             costPT.soc=sum(costPT.soc),costPT.int=sum(costPT.int),
             costACF.soc.hh=sum(costACF.soc.hh),
             costACF.soc.nhh=sum(costACF.soc.nhh),
             costACF.int=sum(costACF.int.hh),
             ## activities
             prophh=mean(hhc/ptentry),
             traceperpt=mean(traceperhhcpt),
             totindexscreen=mean(totindexscreen),
             Screened=mean(Screened),
             Presumptive=mean(Presumptive),
             Diagnosed=mean(Diagnosed),
             ## unit costs
             socu_a.thct=mean(socu_a.thct),
             socu_a.hct=mean(socu_a.hct),
             socu_a.tpt=mean(socu_a.tpt),
             socu_a.tbe=mean(socu_a.tbe),
             socu_a.att=mean(socu_a.att),
             intu_a.thct=mean(intu_a.thct),
             intu_a.hct=mean(intu_a.hct),
             intu_a.tpt=mean(intu_a.tpt),
             intu_a.tbe=mean(intu_a.tbe),
             intu_a.att=mean(intu_a.att)
             ),
          by=.(country)]


ptsm <- melt(PTS,id='country')

looka <- unique(ptsm$variable)

look <- grep('cost',looka,value=TRUE)
## overall costs
ggplot(ptsm[variable %in% look],aes(country,value,
                                    col=variable,group=variable))+
    geom_point()+
    geom_line()+
    ggtitle('Total costs, SoC & INT (relative)')+
    scale_y_sqrt() + ylab('Relative value (sqrt scale)')

ggsave(here('graphs/drivers_pt1.pdf'),w=6,h=5)

look <- grep('u_',looka,value=TRUE)
looks <- grep('soc',look,value=TRUE)
looki <- grep('int',look,value=TRUE)

## unit costs (SOC)
ggplot(ptsm[variable %in% looks],aes(country,value,
                                    col=variable,group=variable))+
    geom_point()+
    geom_line()+
    ggtitle('Unit costs, SoC')+
    scale_y_sqrt() + ylab('Relative value (sqrt scale)')

ggsave(here('graphs/drivers_pt2us.pdf'),w=6,h=5)

## unit costs (INT)
ggplot(ptsm[variable %in% looki],aes(country,value,
                                     col=variable,group=variable))+
    geom_point()+
    geom_line()+
    ggtitle('Unit costs, INT')+
    scale_y_sqrt() + ylab('Relative value (sqrt scale)')

ggsave(here('graphs/drivers_pt2ui.pdf'),w=6,h=5)




## activities etc
look <- looka[!grepl('cost|u_',looka)]
ggplot(ptsm[variable %in% look],aes(country,value,
                                    col=variable,group=variable))+
    geom_point()+
    geom_line()+
    ggtitle('Cascade & effect variables')+
    scale_y_sqrt() + ylab('Relative value (sqrt scale)')

ggsave(here('graphs/drivers_pt3a.pdf'),w=6,h=5)


CD[grep('hct',act)]
corfac #largest

fn <- here('indata/resource.int.csv')
RI <- fread(fn)

names(RI)[names(RI)=='CDI'] <- "Cote d'Ivoire"
RIM <- melt(RI,id='metric')
names(RIM)[2] <- 'country'
RIM <- RIM[metric!='Number of sites reporting']

ptv <- c('Number of Index cases with contact tracing done',
         'PT initiation among HIV entry point',
         'PT initiation among contacts',
         'PT initiation all')
ptvn <- c('PThhcu5','PThhcu5pc','PTHIVentryu5',
          'PTHIVentryu5pc','Ptcompletepc')
attv <- c("Screened for symptoms","Presumptive TB identified",
          "Presumptive TB tested on Xpert","Diagnosed with TB",
          "Treated for DS-TB")
attvn <- c("bacpos","bacposu5","bacposu5pc","TxDenom","TxSuccess",
           "TxSuccesspc")
otherv <- c("CXRamongclindx")

## NOTE correction factor for
corfac <- RIM[metric %in% c(ptv[1],attv[1])]




corfac <- dcast(corfac,country ~ metric,value.var = 'value')
corfac[,scale:=`Screened for symptoms` / `Number of Index cases with contact tracing done`]
corfac <- corfac[country %in% cns,.(country,scale)]-


## things related to correction factor
ggplot(corfac,aes(metric,value,
                  col=country,group=country))+
    geom_point()+
    geom_line()+
    ## facet_wrap(~country,ncol=1)+
    ggtitle('HH & clinic screening')+
    coord_flip()+
    scale_y_sqrt() + ylab('Value (sqrt scale)')


ggsave(here('graphs/drivers_pt4cor.pdf'),w=12,h=5)
