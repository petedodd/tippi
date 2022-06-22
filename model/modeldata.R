## data preparation and cascade modelling
rm(list=ls())
library(here)
library(glue)
library(ggplot2)
library(data.table)
library(ggthemes)
library(scales)
library(readxl)
library(discly)
rexel <- function(x,...) as.data.table(read_excel(x,...))
ssum <- function(x) sqrt(sum(x^2))
gh <- function(x) glue(here(x))

set.seed(1234)
cns <- c("Cameroon","Cote d'Ivoire","DRC","Kenya",
         "Lesotho","Malawi","Tanzania","Uganda","Zimbabwe")
cnisos <- c('CMR','CIV','COD','KEN','LSO','MWI','TZA','UGA','ZWE')


## ===== COUNTRY KEY
## NOTE if inference hasn't been run, skip to SKIP below and load files before continuing
countrykey <- data.table(iso3=cnisos,country=cns)
setkey(countrykey,iso3)
countrykey

save(countrykey,file=here('data/countrykey.Rdata'))

## ====== EMPIRICAL EFFECT DATA
## px & tx data
edat <- list()
for(qty in c('px','tx')){
  for(age in c('04','514')){
    fn1 <- gh('../inference/outdata/siteeffects_{qty}_{age}.csv')
    fn2 <- gh('../inference/outdata/countryeffects_{qty}_{age}.csv')
    SE <- fread(fn1)
    CE <- fread(fn2)
    SE <- SE[int.number>0 & is.finite(site.effect),.(sdl=sd(log(site.effect))),by=country]
    CE <- merge(CE,SE,by='country')
    CE[,le:=log(country.effect)]

    ## use log normal sample with country mean and site variance
    SS <- list()
    for(cn in CE[,country]){
      SS[[cn]] <- data.table(country=cn,
                             id=1:1e4,
                             RR=exp(rnorm(1e4,
                                          CE[country==cn,le],
                                          sd = CE[country==cn,sdl])))
    }
    SS <- rbindlist(SS)
    SS <- dcast(SS,id~country,value.var = 'RR')
    SS[,quant:=qty]
    if(age=='04')
      SS[,age:='0-4']
    else
      SS[,age:='5-14']
    ind <- glue(qty) + '.' + age
    edat[[ind]] <- SS
  }
}
edat <- rbindlist(edat)
edatm <- melt(edat,id.vars = c('id','quant','age'))
edatm <- as.data.table(edatm)
edat <- edatm[,.(id,quant,age,country=variable,RR=value)]

save(edat,file=here('data/edat2.Rdata'))

## ====== EFFECT DATA FROM INFERENCE (new)
## px & tx data
edat <- list()
for(qty in c('px','tx')){
  for(age in c('04','514')){
    fn1 <- glue(here('../inference/outdata/MT_'))
    fn1 <- fn1 + qty + '_' + age + '.Rdata'
    load(fn1)
    colnames(MT) <- countrykey$country
    MT <- as.data.table(MT)
    MT <- MT[sample(nrow(MT),1e4,replace=TRUE)] #resample to 1e4
    MT[,id:=1:nrow(MT)]
    SS <- melt(MT,id='id')
    SS[,quant:=qty]
    if(age=='04')
      SS[,age:='0-4']
    else
      SS[,age:='5-14']
    ind <- glue(qty) + '.' + age
    edat[[ind]] <- SS
  }
}
edat <- rbindlist(edat)
edat <- edat[,.(id,quant,age,country=variable,RR=value)]

save(edat,file=here('data/edat3.Rdata'))


## start here (SKIP)
load(file=here('data/edat3.Rdata'))
load(file=here('data/countrykey.Rdata'))

## ## ====== EFFECT DATA FROM INFERENCE
## ## px & tx data
## edat <- list()
## for(qty in c('px','tx')){
##     for(age in c('04','514')){
##         fn1 <- glue(here('../inference/outdata/'))
##         fn1 <- fn1 + qty + 'SS' + age + '.Rdata'
##         fn2 <- glue(here('../inference/data/'))
##         fn2 <- fn2 + qty + 'CK.' + age + '.Rdata'
##         load(fn1)
##         load(fn2)
##         CK <- CK[order(cno)]
##         names(SS) <- CK$country
##         SS[,quant:=qty]
##         SS[,id:=1:nrow(SS)]
##         if(age=='04')
##             SS[,age:='0-4']
##         else
##             SS[,age:='5-14']
##         ind <- glue(qty) + '.' + age
##         edat[[ind]] <- SS
##     }
## }
## edat <- rbindlist(edat)
## edatm <- melt(edat,id.vars = c('id','quant','age'))
## edatm[,RR:=exp(value)]

## edat <- edatm[,.(id,quant,age,country=variable,RR)]

## save(edat,file=here('data/edat.Rdata'))
## load(file=here('data/edat.Rdata'))

## ===== DISCOUNTED LIFE-YEARS TABLES
## NOTE discount rate set here
discount.rate <- 0.03
## make life-years
fn <- here('data/LYK.Rdata')

if(TRUE){## if(!file.exists(fn)){
    ## calculate discounted life-years
    ## template:
    LYT <- data.table(age=0:14,
                      age_group=c(rep('0-4',5),rep('5-14',10)),
                      LYS=0.0,LYS0=0.0)
    ## make country/age key
    LYK <- list()
    for(iso in cnisos){
        ## iso <- cn
        tmp <- copy(LYT)
        tmp[,iso3:=iso]
        for(ag in tmp$age) tmp[age==ag,LYS:=discly(iso,ag,2020,
                                                   dr=discount.rate)]
        for(ag in tmp$age) tmp[age==ag,LYS0:=discly(iso,ag,2020,
                                                   dr=0)]
        LYK[[iso]] <- tmp
    }
    LYK <- rbindlist(LYK)
    ## assume unweighted & collapse
    LYK <- LYK[,.(LYS=mean(LYS),LYS0=mean(LYS0)),
               by=.(iso3,age=age_group)]
    setkey(LYK,age)
    save(LYK,file=fn)
} else {

    load(file=fn)
}


## ===== CASCADE DATA

## screening by entry-point
SBEP <- fread(here('indata/screened_by_entrypoint.csv'))
SBEP[,iso3:=c("CMR","CIV","COD","KEN","LSO","MWI","UGA","TZA","ZWE")] #NOTE different order
SBEP <- SBEP[order(iso3)]
save(SBEP,file=here('data/SBEP.Rdata'))

## blextract1.csv  blextract2.csv  resoure.int.csv
## baseline cascade data TODO from report
## bl extract 1:
## pt, pthiv from table O3.3 pg 15
fn <- here('indata/blextract1.csv')
B1 <- fread(fn)
## bl extract 2:
## see figure 10 pg 19 of report for PT
## see figure 10 pg 18 of report for ATT
fn <- here('indata/blextract2.csv')
B2 <- fread(fn)
## cascades & resources during intervention
## - from 'cascade data for Nyasha' spreadsheet
fn <- here('indata/resource.int.csv')
RI <- fread(fn)
str(RI)
## correctype <- function(x) if(is.character(x[1])){ as.numeric(gsub(",","",x))} else {x}
## rnmz <- names(RI)[c(-1)]
## RI[,(rnmz) := lapply(.SD,correctype),.SDcols=rnmz] #remove "," and convert char to num


## melted inverention cascade data
names(RI)[names(RI)=='CDI'] <- "Cote d'Ivoire"
RIM <- melt(RI,id='metric')
names(RIM)[2] <- 'country'
RIM <- RIM[metric!='Number of sites reporting']

RIM[,unique(metric)]

## NOTE make sure pc are correct for way used
## item row names
## PT breakdown by entry point
ptv <- c('Number of Index cases with contact tracing done',
         'PT initiation among HIV entry point',
         'PT initiation among contacts',
         'PT initiation all')
## PT breakdown by age
ptvn <- c('PThhcu5','PThhcu5pc','PTHIVentryu5',
          'PTHIVentryu5pc','Ptcompletepc')
## ATT cascade
attv <- c("Screened for symptoms","Presumptive TB identified",
          "Presumptive TB tested on Xpert","Diagnosed with TB",
          "Treated for DS-TB")
## bac+ & tx success
attvn <- c("bacpos","bacposu5","bacposu5pc","TxDenom","TxSuccess",
           "TxSuccesspc")
otherv <- c("CXRamongclindx")

## restrict to particular intervention (ie ATT or TPT focussed)
PTT <- RIM[metric %in% ptv]
ATT <- RIM[metric %in% attv]
atv <- RIM[1:5,metric]
ATT$metric <- factor(ATT$metric,levels=atv,ordered = TRUE)

## ATT <- ATT[country %in% cns] #restrict
## PTT <- PTT[country %in% cns] #restrict

## NOTE setting CIV to no HIV-entry point PT? TODO check
PTT[is.na(value),value:=0]

## absolute
ggplot(ATT,aes(metric,value,fill=metric)) +
    geom_bar(stat = 'identity',position='dodge') +
    facet_wrap(~country,scales='free')+
    scale_fill_colorblind() +
    xlab('') + scale_y_log10(label=comma)+ ylab('Nuber (log scale)')+
    theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))

ggsave(filename=here('graphs/cascade_intervention_abs.pdf'),
       w=11,h=10)
ggsave(filename=here('graphs/cascade_intervention_abs.png'),
       w=11,h=10)

## absolute
ggplot(PTT,aes(metric,value,fill=metric)) +
    geom_bar(stat = 'identity',position='dodge') +
    facet_wrap(~country,scales='free')+
    scale_fill_colorblind() +
    xlab('') + scale_y_log10(label=comma)+ ylab('Nuber (log scale)')+
    theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))

ggsave(filename=here('graphs/cascade_PTintervention_abs.pdf'),
       w=11,h=10)
ggsave(filename=here('graphs/cascade_PTintervention_abs.png'),
       w=11,h=10)

## amount of contact tracing per contact started PT
PTC <- dcast(PTT,country~metric)
PTC[,traceperhhcpt:=`Number of Index cases with contact tracing done`/`PT initiation among contacts`]
PTC <- PTC[,.(country,traceperhhcpt)]

save(PTC,file=here('data/PTC.Rdata'))

## what proportion of PT is via HIV?
PTFH <- PTT[metric %in% ptv[c(2,4)]]
PTFH[,tot:=sum(value),by=country]
PTFH <- PTFH[metric==ptv[2]]
PTFH <- PTFH[,.(country,ptinhiv=value/tot,hiv=value,allpt=tot)]

save(PTFH,file=here('data/PTFH.Rdata'))

## PT age splits
PTTn <- RIM[metric %in% ptvn]
PTTn <- PTTn[country %in% cns]
PTTn[is.na(value),value:=0] #CID HIV=0 NOTE 

save(PTTn,file=here('data/PTTn.Rdata'))

## ATT success and tx chars
ATTn <- RIM[metric %in% attvn]
ATTn <- ATTn[country %in% cns]

save(ATTn,file=here('data/ATTn.Rdata'))

## intervention cascade
## relative
ATR <- copy(ATT)
ATR[metric=='Treated for DS-TB',tx:=value,by=country]
ATR[,frac:=value/tx[5],by=country]
ATR <- ATR[metric!='Treated for DS-TB']

ggplot(ATR,aes(metric,frac,color=country,group=country)) +
    geom_point() + geom_line()+
    xlab('') +
    scale_y_log10()+
    ylab('Number per child treated for TB (log scale)')+
    theme_classic() + ggpubr::grids()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))

ggsave(filename=here('graphs/intervention_cascade.pdf'),
       w=6,h=5)
ggsave(filename=here('graphs/intervention_cascade.png'),
       w=6,h=5)


ATR <- ATR[,.(metric,country,frac)]

save(ATR,file=here('data/ATR.Rdata')) #cascade per ATT initiation

load(file=here('data/ATR.Rdata')) #cascade


## output cascade data: costs
CD <- fread(here('indata/unit_costs_updated_2022-05-25.csv'))
CD[,iso3:=gsub('DRC','COD',Country)]
CD[,iso3:=gsub('CDI','CIV',iso3)]
CD[,V1:=NULL]
save(CD,file=here('data/CD.Rdata')) #cascade

load(file=here('data/CD.Rdata')) #cascade
## mapping costs to activities
CD[,metric:=Activity]
CD[,unique(Activity)]
ATR[,unique(metric)]
CD[,metric:=NA]
CD[,.(Country,Activity,uc.soc,uc.soc.sd,uc.int)]
ATR

## NOTE adding sample collection costs on to Xpert costs
tmp <- CD[Activity=="Sample collection",
          .(Country,uc.soc1=uc.soc,uc.soc.sd1=uc.soc.sd,uc.int1=uc.int)]
CD <- merge(CD,tmp,by='Country')
CD[Activity=="Xpert testing",.(Country,uc.soc,uc.soc.sd,uc.int)] #check
CD[Activity=="Xpert testing", #add on
   c('uc.soc','uc.soc.sd','uc.int'):=
       .(uc.soc+uc.soc1,
         sqrt(uc.soc.sd^2+uc.soc.sd1^2),
         uc.int+uc.int1)]
CD[,c('uc.soc1','uc.soc.sd1','uc.int1'):=NULL] #remove temporary data
CD[Activity=="Xpert testing",.(Country,uc.soc,uc.soc.sd,uc.int)] #check

## --- splitting screening types
unique(CD[,.(Activity,metric)]) #metric is from ATR, Activity is from CD
## NOTE in ATR, only have "Screened for symptoms" -- needs splitting using SBEP

SBEPm <- melt(SBEP[,.(iso3,CBhhcm,FBhhcm,HIVicf,nonHIVicf)],id='iso3')
SBEPm[,total:=sum(value),by=iso3]
SBEPm[,epfrac:=value/total] #entry-point fraction for screening
SBEPm[,sum(epfrac),by=iso3] #check

## merge cateogries
SBEPm[variable=='CBhhcm',Activity:='Community hhci']
SBEPm[variable=='FBhhcm',Activity:='Facility hhci']
SBEPm[variable=='HIVicf',Activity:='Screening in HIV clinic']
SBEPm[variable=='nonHIVicf',Activity:='Screening in non-HIV clinic']

## do merge, add epfrac =1 for rest
CD <- merge(CD,SBEPm[,.(iso3,Activity,epfrac)],
            by=c('iso3','Activity'),all=TRUE)
CD[is.na(epfrac),epfrac:=1]

## --- mapping between metrics and activities
screenacts <- CD[,grep("creen|hhci",unique(Activity),value=TRUE)]

CD[Activity %in% screenacts,
   metric:="Screened for symptoms"] #NM/PJD -- adding all screenings

CD[Activity=='Presumptive TB evaluation',
   metric:="Presumptive TB identified"] #NM

CD[Activity=='Xpert testing',
   metric:="Presumptive TB tested on Xpert"]

CD[Activity=='TB treatment',
   metric:="TB treatment"]

CD[Activity=='TPT treatment',
   metric:="TPT treatment"]

## NOTE these are not in below and dropped from cascade merge
## OK for sample collection which is already buily into Xpert testing

CD[Activity=='Sample collection',
   metric:="Sample collection"]

CD[Activity=='Chest x-ray',
   metric:="Chest x-ray"]

unique(CD[,.(Activity,metric)]) #check

## NOTE need to change this column type to NA real from logical
CD[,uc.int.sd:=NULL]
CD[,uc.int.sd:=NA_real_]

## --- build X-ray costs in (see row in resource.int.csv) ---
RIM[metric=='CXRamongclindx'] #NOTE this is among presumptives, but the denominator is really those bac-ve
## this means that as a pc applied to /all/ TB it is likely an overestimate (and conservative for CEA)

tmp <- CD[metric=='Chest x-ray']
tmp <- merge(tmp,countrykey,by='iso3')
tmp <- merge(tmp,RIM[metric=='CXRamongclindx',.(country,value)],by='country')
## multiply X-ray cost by proportion of time used
tmp[,c('E.uc.soc','E.uc.int','E.uc.soc.sd','E.uc.int.sd'):=.(uc.soc*value,uc.int*value,uc.soc.sd*value,0)]
tmp <- tmp[,.(iso3,E.uc.soc,E.uc.int,E.uc.soc.sd,E.uc.int.sd)]
tmp[,metric:='Presumptive TB identified'] #what this will be added to

CD <- merge(CD,tmp,by=c('iso3','metric'),all.x=TRUE)
CD[metric=='Presumptive TB identified',
   c('uc.soc','uc.soc.sd','uc.int','uc.int.sd'):=
        .(uc.soc+E.uc.soc,
          sqrt(uc.soc.sd^2+E.uc.soc.sd^2),
          uc.int+E.uc.int,
          sqrt(0^2+E.uc.int.sd^2))]
CD[metric=='Presumptive TB identified'] #check
CD[,c('E.uc.soc','E.uc.soc.sd','E.uc.int','E.uc.int.sd'):=NULL] #remove dummy data
## --- X-ray end ---

## aggregate over metrics before merge
CD[is.na(uc.soc.sd),uc.soc.sd:=0]
CD[is.na(uc.int.sd),uc.int.sd:=0]

## summing means and sqrtsumsqr for sd
CD <- CD[,.(uc.soc=sum(uc.soc*epfrac),uc.soc.sd=ssum(uc.soc.sd*epfrac),
            uc.int=sum(uc.int*epfrac),uc.int.sd=ssum(uc.int.sd*epfrac)
            ), by=.(iso3,metric)]


## merge into cascade
ATR <- merge(ATR,countrykey,by='country')

ART2 <- rbind(ATR[,.(iso3,metric,frac)],
              data.table(iso3=cnisos,metric='TB treatment',frac=1.0)) #TB treatment added to cascade

## adding dx, which isn't in itself associated with resources
CD2 <- rbindlist(list(
    CD[,.(iso3,metric,uc.soc,uc.soc.sd,uc.int)],
    data.table(iso3=cnisos,metric='Diagnosed with TB',
               uc.soc=0.0,uc.soc.sd=0.0,uc.int=0.0)))

## NOTE this drops unused init costs (X-ray, sample, TPT)
ART2 <- merge(ART2,CD2,by=c('iso3','metric'))
ART2[iso3=='MWI']

save(ART2,file=here('data/ART2.Rdata')) #cascade + costs

load(file=here('data/ART2.Rdata')) #cascade + costs

## ==================== baseline data ===========
B1
B2 #NOTE this is aggregated over country

## --- compare cascade
## baseline
D <- B1[!is.na(prtb)]
D <- D[,.(country,prtb,xpert,Tbdx)]
names(D)[2:4] <- atv[2:4]
D <- melt(D,id='country')
D[,period:='Baseline']
D[variable=='Diagnosed with TB',dx:=value,by=country]
D <- D[order(country,variable)]
D[,frac:=value/dx[3],by=country]
## D <- D[variable!='Diagnosed with TB']


## intervention
D2 <- ATT[metric!='Treated for DS-TB']
D2[metric=='Diagnosed with TB',dx:=value,by=country]
D2[,frac:=value/dx[4],by=country]
## D2 <- D2[metric!='Diagnosed with TB']
D2[,period:='Intervention']


## join
DB <- rbind(D2[,.(country,metric,frac,period)],
            D[,.(country,metric=variable,frac,period)])

vrs <- DB[,unique(metric)]
DB[,metric:=factor(metric,levels=vrs,ordered = TRUE)]

DBR <- DB[country %in% cns]

ggplot(DBR,aes(metric,frac,col=country,group=paste(period,country),
              lty=period)) +
    geom_point() + geom_line() +
    scale_y_log10(label=comma)+ xlab('')+
    ylab('Number per children diagnosed for TB (log scale)')+
    theme_classic() + ggpubr::grids()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))

## ggsave(filename=here('graphs/cascade_compare.pdf'),w=7,h=8)
ggsave(filename=here('graphs/cascade_compare.png'),w=7,h=8)


DBR <- DBR[metric %in% vrs[2:3]]

ggplot(DBR,aes(period,frac,col=country,lty=metric,
               group=paste(country,metric))) +
    geom_point() + geom_line()+
    facet_wrap(~metric) +
    ylab('Number per child diagnosed with TB') + xlab('Period')## +
    ## theme_classic() + ggpubr::grids()

## ggsave(filename=here('graphs/cascade_compare2.pdf'),w=8,h=5)
ggsave(filename=here('graphs/cascade_compare2.png'),w=8,h=5)

DBR

## comparison of increases in resources vs effect
## edatm <- edat[,.(RR=mean(RR)),by=.(quant,age,country)]

DBC <- dcast(DBR,country+metric ~ period,value.var = 'frac')
DBC[,ratio:=Intervention/Baseline]
## compute average for each metric to fill in gaps
tmp <- DBC[,.(rat=mean(ratio,na.rm=TRUE)),by=metric]
DBC <- merge(DBC,tmp,by='metric')
DBC[is.na(Baseline),ratio:=rat] #replace missing with average
DBC[,rat:=NULL]                 #remove additional data
tmp <- DBC[metric=='Presumptive TB identified']
tmp[,metric:="Screened for symptoms"]
DBC <- rbind(DBC,tmp) #NOTE apply presumptive to screened also


## this gives the int/soc ratio of presumptive and testing per treatment
save(DBC,file=here('data/DBC.Rdata')) #ratios int v bl
load(file=here('data/DBC.Rdata')) #ratios int v bl

## === cost tables
DBC <- merge(DBC,countrykey,by='country')
ART2 <- merge(ART2,DBC,by=c('iso3','metric'),all.x=TRUE)
ART2[,c('country','Baseline','Intervention'):=NULL]
ART2[is.na(ratio),ratio:=1.0]
ART2[,.(iso3,metric,ratio)]

## merge and compute both sets of costs
xtra <- ART2[,.(cost.soc=sum(uc.soc*frac/ratio),
                cost.soc.sd=ssum(frac*uc.soc.sd/ratio),
                cost.int=sum((uc.int+uc.soc)*frac),
                cost.int.sd=ssum(frac*uc.soc.sd)
                ),
             by=iso3]

xtra[,frac.soc:=paste0(round(cost.soc),' (',
                       round(cost.soc.sd),')')]
xtra[,frac.int:=paste0(round(cost.int),' (',
                       round(cost.int.sd),')')]
xtra <- melt(xtra[,.(iso3,frac.soc,frac.int)],id='iso3')
xtra[,metric:='cost']

cascadetab.soc <- rbind(ART2[,.(iso3,metric,
                                frac=paste0(
                                    round(frac/ratio,2)))],
                        xtra[variable=='frac.soc',
                             .(iso3,metric,frac=value)])
cascadetab.int <- rbind(ART2[,.(iso3,metric,
                                frac=paste0(round(frac,2)))],
                        xtra[variable=='frac.int',
                             .(iso3,metric,frac=value)])


cascadetab.soc <- dcast(cascadetab.soc,iso3~metric,
                        value.var = 'frac')
cascadetab.int <- dcast(cascadetab.int,iso3~metric,
                        value.var = 'frac')

## this block is for table 1 below
ctsm <- melt(cascadetab.soc,id='iso3')
ctsm[,condition:='SOC']
ctim <- melt(cascadetab.int,id='iso3')
ctim[,condition:='INT']
ctm <- rbind(ctsm,ctim)

ccs <- c("iso3","Screened for symptoms",
         "Presumptive TB identified",
         "Presumptive TB tested on Xpert",
         "Diagnosed with TB",
         "TB treatment","cost")
setcolorder(cascadetab.soc,ccs)
setcolorder(cascadetab.int,ccs)
cascadetab.int[,iso3:=NULL]
cascadetab <- cbind(cascadetab.soc,cascadetab.int)
cascadetab

fwrite(cascadetab,file=here('outdata/cascadetab.csv'))

## ------ reformatting for upper part of Table 1 -----
ctm <- merge(ctm,countrykey,by='iso3')
ctm[,iso3:=NULL]
ctm <- ctm[order(country,condition,variable)]
ctmc <- dcast(ctm,condition + variable ~ country,value.var = 'value')
ctmc$condition <- factor(ctmc$condition,levels=c('SOC','INT'),ordered=TRUE)
ctmc$variable <- factor(ctmc$variable,levels=ccs[-1],ordered=TRUE)
setkey(ctmc,condition,variable)
Table1ATT <- copy(ctmc)

save(Table1ATT,file=here('data/Table1ATT.Rdata'))



## NOTE below here is WHO data and doesn't relate to costing
## === making CDRs and getting some relevant HIV data ===
whodir <- glue('~/Dropbox/Documents/WHO_TBreports/data2020/')

## === load WHO notifications
## http://www.who.int/tb/country/data/download/en/
fn <- here('indata/TB_notifications_2020-10-15.csv')
N <- fread(fn)

nmz <- paste0('newrel_',c('m04','f04',
                          'm59','f59',
                          'm1014','f1014',
                          'm014','f014'))
nmz <- c('iso3','year',nmz)

## reduce to relevant data: 2018 has no NA
NP <- N[year==2018,..nmz]
NP <- melt(NP,id.vars = c('iso3','year'))
NP[,sex:=ifelse(grepl("f",variable),'F','M')]
NP[,age:=gsub("[a-z]|_","",variable)]
NP <- NP[iso3 %in% cnisos]
NP
NP <- NP[age %in% c('04','014')]
NP


## === load WHO age-specific incidence estimates
fn <- here('indata/TB_burden_age_sex_2020-10-15.csv')
A <- fread(fn)
## keep only relevant categories
A <- A[year==2019]
A <- A[sex!='a']
A <- A[age_group %in% c('0-4','0-14','5-14')]
A <- A[risk_factor=='all']
A[,age:=gsub('-','',age_group)]
## harmonize namings
A[sex=='f',sex:='F']
A[sex=='m',sex:='M']
unique(A[,.(sex,age)])                  #check
A[,best.sd:=(hi-lo)/3.92]
A <- A[iso3 %in% cnisos]
A

## HIV
fn <- here('indata/TB_burden_countries_2020-10-15.csv')
H <- fread(fn)
H <- H[year==2019,.(iso3,e_tbhiv_prct,e_tbhiv_prct_lo,e_tbhiv_prct_hi)]
H <- H[iso3 %in% cnisos]
H[,hiv:=e_tbhiv_prct/100]
H[,hiv.sd:=(e_tbhiv_prct_hi-e_tbhiv_prct_lo)/392] #adults of course
H <- H[,.(iso3,hiv,hiv.sd)]
hfn <- here('data/H.Rdata')
save(H,file=hfn)


## === merge data
AN <- merge(NP[,.(iso3,sex,age,notes=value)],
            A[,.(iso3,sex,age,inc=best,lo,hi)],
            by=c('iso3','sex','age'),all.x=TRUE,all.y=FALSE)


ANO <- AN[age=='014']
ANY <- AN[age=='04']
ANB <- merge(ANY[,.(iso3,sex,notes.04=notes,inc.04=inc,
                    inc.sd.04=(hi-lo)/3.92)],
             ANO[,.(iso3,sex,notes.514=notes,inc.514=inc,
                    inc.sd.514=(hi-lo)/3.92)],
             by=c('iso3','sex')
             )
CDRu <- ANB[,.(rel.sd=mean(inc.sd.514/inc.514)),by=iso3]
ANB[,c('notes.514','inc.514'):=.(notes.514-notes.04,inc.514-notes.04)]
CDR <- ANB[,.(notes.04=sum(notes.04),inc.04=sum(inc.04),
              notes.514=sum(notes.514),inc.514=sum(inc.514)),by=iso3]
CDR[,c('cdr04','cdr514'):=.(notes.04/inc.04,notes.514/inc.514)]
CDR[,totnotes:=notes.04+notes.514]
CDR[,c('frac04','frac514'):=.(notes.04/totnotes,notes.514/totnotes)]
CDR <- melt(CDR[,.(iso3,frac04,frac514,cdr04,cdr514)],id='iso3')
CDR[,qty:='cdr']
CDR[grepl('frac',variable),qty:='frac']
CDR[,age:='0-4']
CDR[grepl('14',variable),age:='5-14']
CDR <- CDR[,.(iso3,qty,age,value)]
CDR <- merge(CDR,CDRu,by='iso3',all.x=TRUE)
CDR[,cdr.v:=NA_real_]
CDR[,cdr.v:=(rel.sd*value)^2]
CDR[,rel.sd:=NULL]
CDR <- CDR[iso3 %in% cnisos]
cdrfn <- here('data/CDR.Rdata')
save(CDR,file=cdrfn)


## Figure 1
## HHCM ACF
## 
## 0-4 screened at facility = 4661 (60.1%)
## 0-4 presumptive TB at facility =725 (48.2%)
## 0-4 diagnosed with TB at facility= 110 (49.3%)
## 
## 0-4 screened at community = 11285 (46.8%)
## 0-4 presumptive TB at community =725 (38.0%)
## 0-4 diagnosed with TB at community= 244 (44.2%)
## 
## 
## 0-14 screened at facility = 7422
## 0-14 presumptive TB at facility = 1503
## 0-14 diagnosed with TB at facility= 223
## 
## 0-14 screened at community = 24093
## 0-14 presumptive TB at community = 1546
## 0-14 diagnosed with TB at community= 558

## 47139
## 25944
## 17990
## 13533
## 31515

## HHCM cascade across countries
HHCMtop <- data.table(
    value=c(47139,
            25944,
            17990,
            13533,
            31515
            ),
    activity=c(
        "index cases",
        "estimated B+ PTB index cases",
        "estimated TPT-eligible contacts <5 years",
        "index cases with HHCM",
        "household contacts screened"
    )
)
HHCMtop[,c('age','mode'):=NA]

aty <- c('Screened','Presumptive','Diagnosed')
mds <- c('facility','household')
ags <- c('0-4','0-14','5-14')
yf <- c(4661,725,110)
yc <- c(11285,725,244)
df <- c(7422,1503,223)
dc <- c(24093,1546,558)
mdsl <- rep(mds,each=3)
mdsl <- rep(mdsl,3)

HHCM <- data.table(activity=rep(aty,2*3),
                   mode=mdsl,
                   age=rep(ags,each=2*3),
                   value=c(yf,yc,df,dc,df-yf,dc-yc)
                   )

HHCM <- rbind(HHCMtop,HHCM)

fwrite(HHCM,file=here('data/HHCM.csv'))
save(HHCM,file=here('data/HHCM.Rdata'))


## ===== SIMPLE CONVERSIONS FROM CSV TO RDATA FOR CONSISTENCY
## BL extract HIV
BL <- fread(here('indata/blextract1.csv'))
save(BL,file=here('data/BL.Rdata'))
## INT cascade data
INT <- fread(here('indata/resource.int.csv'))
names(INT)[names(INT)=='CDI'] <- "Cote d'Ivoire"
save(INT,file=here('data/INT.Rdata'))
## CE thresholds
CET <- fread(here('indata/TIPPIresults - CEthresholds.csv'))
CET[,`1x GDP`:=as.numeric(gsub(",","",`1x GDP`))] ## format threshold data
CET[,`3x GDP`:=as.numeric(gsub(",","",`3x GDP`))]
CETM <- CET[,.(iso3=Country,`1x GDP`,`3x GDP`,Y1a,Y1b,Y2a,Y2b)]
CETM <- melt(CETM,id='iso3')
CETM[,value:=as.numeric(value)]
CETM <- merge(CETM,countrykey,by='iso3')
names(CETM)[2] <- 'threshold'
tmp <- CETM[threshold=='1x GDP'] #add in 1/2 GDP as in Chi Vassall
tmp[,value:=value/2]
tmp[,threshold:='0.5x GDP']
CETM <- rbind(CETM,tmp)
save(CETM,file=here('data/CETM.Rdata'))
## modelling parmeters
PD <- read.csv(here('indata/TIPPIresults - TIPPIparms.csv'))
save(PD,file=here('data/PD.Rdata'))


## --- data for ASM age splits file & BC
D <- fread(here('../dataprep/outdata/D.csv'))
T <- fread(here('../dataprep/outdata/T.csv'))
P <- fread(here('../dataprep/outdata/P.csv'))


## BC rate ratios
TC <- T[,.(ATT.Baseline.Num=sum(Baseline.Num),
           ATT.Intervention.Num=sum(Intervention.Num),
           Baseline.FT=sum(Baseline.FT),
           Intervention.FT=sum(Intervention.FT)),by=.(Country)]
PC <- P[,.(PT.Baseline.Num=sum(Baseline.Num),
           PT.Intervention.Num=sum(Intervention.Num),
           Baseline.FT=sum(Baseline.FT),
           Intervention.FT=sum(Intervention.FT)),by=.(Country)]
TC[,ATT.Baseline.Rate:=ATT.Baseline.Num/Baseline.FT];
TC[,ATT.Intervention.Rate:=ATT.Intervention.Num/Intervention.FT];
PC[,PT.Baseline.Rate:=PT.Baseline.Num/Baseline.FT];
PC[,PT.Intervention.Rate:=PT.Intervention.Num/Intervention.FT];

BC <- merge(TC[,.(Country,ATT.Baseline.Num,ATT.Intervention.Num,
                  ATT.Baseline.Rate,ATT.Intervention.Rate)],
            PC[,.(Country,PT.Baseline.Num,PT.Intervention.Num,
                  PT.Baseline.Rate,PT.Intervention.Rate)],
            by='Country'
            )

BC <- BC[,.(Country,
            BR=ATT.Baseline.Rate/PT.Baseline.Rate,
            IR=ATT.Intervention.Rate/PT.Intervention.Rate
            )]

save(BC,file=here('data/BC.Rdata'))


## ASM age splits
(P <- P[,.(Baseline=sum(Baseline.Num),
           Intervention=sum(Intervention.Num)),by=.(Country,age)])
(D <- D[,.(Baseline=sum(Baseline.Num),
           Intervention=sum(Intervention.Num)),by=.(Country,age)])
(T <- T[,.(Baseline=sum(Baseline.Num),
           Intervention=sum(Intervention.Num)),by=.(Country,age)])

P[,qty:='pt']
D[,qty:='dx']
T[,qty:='tx']


ASM <- rbindlist(list(P,D,T))
ASM <- melt(ASM,id=c('Country','age','qty'))
ASM <- dcast(ASM,Country + qty + variable ~ age)
names(ASM)[c(1,3)] <- c('country','period')
ASM[,`0-14`:=NULL]

save(ASM,file=here('data/ASM.Rdata'))


## completing Table 1
BCI <- BC[,.(Country,BR=1/BR,IR=1/IR)] #report other way up
Table1PT <- transpose(BCI,make.names = 'Country')

Table1PT[,c('condition','variable'):=.(rep('PT',2),
                                       c('Baseline ratio PT:ATT intiations',
                                         'Intervention ratio PT:ATT intiations'))]

setcolorder(Table1PT,c('condition','variable',BC$Country))

save(Table1PT,file=here('data/Table1PT.Rdata'))

load(file=here('data/Table1PT.Rdata'))



## PT via HIV clinic -- hago?
tmp1 <- transpose(PTFH[,.(country, 1e2*ptinhiv)],make.names = TRUE)
tmp1[,c('condition','variable'):=.('PT','PT via HIV clinic')]
setcolorder(tmp1,names(Table1PT))


## HHCM community based*; *=calculated on screens
tmp2 <- SBEP[,.(country, 1e2*CBhhcm/(CBhhcm+FBhhcm))]

tmp2[country=='CDI',country:="Cote d'Ivoire"]
tmp2 <- transpose(tmp2,make.names = TRUE)
tmp2[,c('condition','variable'):=.('PT','HHCM community-based')]
setcolorder(tmp2,names(Table1PT))

rqty <- c("Number of Index cases with contact tracing done",
          "PT initiation among contacts")
## Started on PT
tmp3 <- RIM[metric %in% rqty[2],.(country,value)]
tmp3 <- transpose(tmp3,make.names = TRUE)
tmp3[,c('condition','variable'):=.('PT',rqty[2])]
setcolorder(tmp3,names(Table1PT))


## Households screened
tmp4 <- RIM[metric %in% rqty[1],.(country,value)]
tmp4 <- transpose(tmp4,make.names = TRUE)
tmp4[,c('condition','variable'):=.('PT',rqty[1])]
setcolorder(tmp4,names(Table1PT))


## Children screened
tmp5 <- SBEP[,.(country, (CBhhcm+FBhhcm))]

tmp5[country=='CDI',country:="Cote d'Ivoire"]
tmp5 <- transpose(tmp5,make.names = TRUE)
tmp5[,c('condition','variable'):=.('PT','Child HH contacts screened')]
setcolorder(tmp5,names(Table1PT))


## Presumptive TB
## NOTE this does this & next, but it is same across countries
## NOTE this is also based on BL data for presum/prevalence - used to model
## TODO query whether best data available by country
## TODO these are the relevant qties used in model
tmpr <- B2[metric %in%
           c('Number of Index cases with contact tracing done',
             'Presumptive TB identified','Diagnosed with TB')]
tmpr[,perHH:=total/
        tmpr[metric=='Number of Index cases with contact tracing done',
             total]]


tmps <- SBEP[,.(country, kids=(CBhhcm+FBhhcm))]


tmps[country=='CDI',country:="Cote d'Ivoire"]
kidsperhh <- merge(tmps,
                   RIM[metric %in% rqty[1],.(country,nhh=value)],
                   by='country')
kidsperhh[,kidsphh:=kids/nhh]
kidsperhh[,tb:=nhh*tmpr[metric=='Diagnosed with TB',perHH]]
kidsperhh[,pr:=nhh*tmpr[metric=='Presumptive TB identified',perHH]]
kidsperhh[,c('tbperkid','prperkid'):=.(tb/kids,pr/kids)]

tmp6 <- transpose(kidsperhh[,.(country,pr)],make.names = TRUE)
tmp6[,c('condition','variable'):=.('PT','Child HH contacts presumed TB')]
setcolorder(tmp6,names(Table1PT))

## Diagnosed TB
tmp7 <- transpose(kidsperhh[,.(country,tb)],make.names = TRUE)
tmp7[,c('condition','variable'):=.('PT','Child HH contacts diagnosed with TB')]
setcolorder(tmp7,names(Table1PT))

Table1PTxtraA <- rbindlist(list(tmp1,tmp2))
Table1PTxtraB <- rbindlist(list(tmp3,tmp4,tmp5,tmp6,tmp7))

## make B table relative to top row
Table1PTxtraB <- cbind(Table1PTxtraB[,.(condition,variable)],
                       Table1PTxtraB[,lapply(.SD,function(x) x/x[1]),
                                     .SDcols=names(Table1PTxtraB)[-c(1:2)]])


Table1PT <- rbind(Table1PT,Table1PTxtraA)
Table1PT <- rbind(Table1PT,Table1PTxtraB)

save(Table1PT,file=here('data/Table1PT.Rdata'))


## Cost per PT initiation, $ (SD)
## TODO -- see model script



## TODO questions
## HIV mix -- different under intervention?
