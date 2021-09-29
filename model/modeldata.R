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

## ====== EFFECT DATA FROM INFERENCE
## px & tx data
edat <- list()
for(qty in c('px','tx')){
    for(age in c('04','514')){
        fn1 <- glue(here('../inference/outdata/'))
        fn1 <- fn1 + qty + 'SS' + age + '.Rdata'
        fn2 <- glue(here('../inference/data/'))
        fn2 <- fn2 + qty + 'CK.' + age + '.Rdata'
        load(fn1)
        load(fn2)
        CK <- CK[order(cno)]
        names(SS) <- CK$country
        SS[,quant:=qty]
        SS[,id:=1:nrow(SS)]
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
edatm[,RR:=exp(value)]

edat <- edatm[,.(id,quant,age,country=variable,RR)]
save(edat,file=here('data/edat.Rdata'))

## load(file=here('data/edat.Rdata'))

## ===== COUNTRY KEY
(cns <- edat[,unique(country)])
## cnisos <- c('CMR','CIV','COD','KEN','LSO','MWI','TZA','UGA','ZWE')
cnisos <- c('CMR','CIV','COD','KEN','LSO','MWI','UGA','ZWE')

countrykey <- data.table(iso3=cnisos,country=cns)
setkey(countrykey,iso3)
countrykey

save(countrykey,file=here('data/countrykey.Rdata'))

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


## ===== TODO up to here
## BL data etc: blextract1.csv  blextract2.csv  resoure.int.csv
fn <- here('indata/blextract1.csv')
B1 <- fread(fn)
fn <- here('indata/blextract2.csv')
B2 <- fread(fn)
fn <- here('indata/resource.int.csv')
RI <- fread(fn)
str(RI)
correctype <- function(x) if(is.character(x[1])){ as.numeric(gsub(",","",x))} else {x}
rnmz <- names(RI)[c(-1)]
RI[,(rnmz) := lapply(.SD,correctype),.SDcols=rnmz] #remove "," and convert char to num

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
corfac <- corfac[country %in% cns,.(country,scale)]
save(corfac,file=here('data/corfac.Rdata'))

## restrict to particular intervention
PTT <- RIM[metric %in% ptv]
ATT <- RIM[metric %in% attv]
atv <- RIM[1:5,metric]
ATT$metric <- factor(ATT$metric,levels=atv,ordered = TRUE)

ATT <- ATT[country %in% cns] #restrict
PTT <- PTT[country %in% cns] #restrict


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

## what proportion of PT is via HIV
PTFH <- PTT[metric %in% ptv[c(2,4)]]
PTFH[,tot:=sum(value),by=country]
PTFH <- PTFH[metric==ptv[2]]
PTFH <- PTFH[,.(country,ptinhiv=value/tot,hiv=value,allpt=tot)]

save(PTFH,file=here('data/PTFH.Rdata'))

PTTn <- RIM[metric %in% ptvn]
PTTn <- PTTn[country %in% cns]
PTTn[is.na(value),value:=0]

save(PTTn,file=here('data/PTTn.Rdata'))

ATTn <- RIM[metric %in% attvn]
ATTn <- ATTn[country %in% cns]

save(ATTn,file=here('data/ATTn.Rdata'))

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
save(ATR,file=here('data/ATR.Rdata')) #cascade

## output cascade data 
CD <- fread(here('indata/TIPPIresults - costs.csv'))
CD[,iso3:=gsub('DRC','COD',Country)]
CD <- CD[iso3!='']
CD[,uc.soc.sd:=abs(uc.soc.sd)]
save(CD,file=here('data/CD.Rdata')) #cascade

load(file=here('data/ATR.Rdata')) #cascade

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

## map unit costs to activities
CD[Activity=='Presumptive TB evaluation',
   metric:="Presumptive TB identified"] #NM
CD[Activity=='TB contact investigation',
   metric:="Screened for symptoms"] #NM
CD[Activity=='Sample collection',
   metric:="Sample collection"] #NM
CD[Activity=='Xpert testing',
   metric:="Presumptive TB tested on Xpert"]
CD[Activity=='TB treatment',
   metric:="TB treatment"]

CD2 <- CD[!is.na(metric)]
ATR <- merge(ATR,countrykey,by='country')

ART2 <- rbind(ATR[,.(iso3,metric,frac)],
              data.table(iso3=cnisos,metric='TB treatment',frac=1.0)) #TB treatment added to cascade

CD2 <- rbindlist(list(
    CD2[,.(iso3,metric,uc.soc,uc.soc.sd,uc.int)],
    data.table(iso3=cnisos,metric='Diagnosed with TB',
               uc.soc=0.0,uc.soc.sd=0.0,uc.int=0.0)))## ,
    ## data.table(iso3=cnisos,metric='Presumptive TB identified',
    ##            uc.soc=0.0,uc.soc.sd=0.0,uc.int=0.0)))


ART2 <- merge(ART2,CD2,
              by=c('iso3','metric'))
ART2[is.na(uc.int),uc.int:=0.0]
ART2[is.na(uc.soc.sd),uc.soc.sd:=0]
ART2[iso3=='MWI']
save(ART2,file=here('data/ART2.Rdata')) #cascade + costs

load(file=here('data/ART2.Rdata')) #cascade + costs

## see after ratio computations
## xtra <- ART2[,.(cost=sum(uc.soc*frac),cost.sd=ssum(frac*uc.soc.sd)),
##              by=iso3]

## xtra[,frac:=paste0(round(cost),' (',round(cost.sd),')')]
## xtra[,metric:='cost']

## cascadetab <- rbind(ART2[,.(iso3,metric,frac=paste0(round(frac,2)))],
##                     xtra[,.(iso3,metric,frac)])


## cascadetab <- dcast(cascadetab,iso3~metric,value.var = 'frac')
## ccs <- c("iso3","Screened for symptoms","Presumptive TB identified",
##          "Presumptive TB tested on Xpert","Diagnosed with TB",
##          "TB treatment","cost")
## setcolorder(cascadetab,ccs)

## fwrite(cascadetab,file=here('outdata/cascadetab.csv'))

## CD <- dcast(CD,iso3~Activity,value.var = c('uc.soc','uc.soc.sd'))

## ==================== baseline data ===========
B1
B2

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

ggsave(filename=here('graphs/cascade_compare.pdf'),w=7,h=8)
ggsave(filename=here('graphs/cascade_compare.png'),w=7,h=8)


DBR <- DBR[metric %in% vrs[2:3]]

ggplot(DBR,aes(period,frac,col=country,lty=metric,
               group=paste(country,metric))) +
    geom_point() + geom_line()+
    facet_wrap(~metric) +
    ylab('Number per child diagnosed with TB') + xlab('Period')## +
    ## theme_classic() + ggpubr::grids()

ggsave(filename=here('graphs/cascade_compare2.pdf'),w=8,h=5)
ggsave(filename=here('graphs/cascade_compare2.png'),w=8,h=5)

DBR

## see SOS cost parms
## TB assesment costs ~ 10 USD
## TB ATT cost ~ 140 USD
## Xpert ~ 20 USD
## PT cost parms

## PT cascade -- currently assume zero baseline


## comparison of increases in resources vs effect
edatm <- edat[,.(RR=mean(RR)),by=.(quant,age,country)]

DBC <- dcast(DBR,country+metric ~ period,value.var = 'frac')
DBC[,ratio:=Intervention/Baseline]
DBC <- DBC[!is.na(ratio)]

edatm <- edatm[quant=='tx',.(RR=mean(RR)),country]
DBC <- merge(DBC,edatm,by='country',all.x = TRUE,all.y=FALSE)

ggplot(DBC,aes(RR,ratio,col=country)) +
    geom_point() +
    facet_wrap(~metric)


DBC <- DBC[,.(country,metric,ratio)]

exns <- ATR[,unique(country)]
## exclude countries that have a ratio
exns <- exns[!exns %in% c(#'Kenya','Lesotho',
                            'Malawi',
                            'Uganda',
                            'Zimbabwe')]
exns <- as.character(exns)
mrcs <- ATR[,unique(metric)]
mrcs <- mrcs[grepl('Presumptive',mrcs)]

DBCE <- data.table(expand.grid(country=c(exns),metric=mrcs))
DBCE[,ratio:=NA_real_]
DBCE <- rbind(DBC,DBCE)
DBCE[metric=='Presumptive TB identified',tbi:=mean(ratio,na.rm=TRUE)]
DBCE[,tbi:=mean(ratio,na.rm=TRUE),by=metric]
DBCE[is.na(ratio),ratio:=tbi] #means for NAs
DBCE[,tbi:=NULL]

save(DBCE,file=here('data/DBCE.Rdata')) #ratios int v bl
load(file=here('data/DBCE.Rdata')) #ratios int v bl


## === cost tables
DBCE <- merge(DBCE,countrykey,by='country')
ART2 <- merge(ART2,DBCE,by=c('iso3','metric'),all.x=TRUE)
ART2[,country:=NULL]
tmp <- DBCE[metric=='Presumptive TB identified']
ART2 <- merge(ART2,tmp[,.(iso3,tt=ratio)],all.x=TRUE,by='iso3')
ART2[metric=='Screened for symptoms',ratio:=tt]
ART2[is.na(ratio),ratio:=1.0]
ART2[,tt:=NULL]


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

## NOTE below here is WHO data and doesn't relate to costing
## === making CDRs and getting some relevant HIV data ===
whodir <- glue('~/Dropbox/Documents/WHO_TBreports/data2020/')

## === load WHO notifications
## http://www.who.int/tb/country/data/download/en/
fn <- whodir + 'TB_notifications_2020-10-15.csv'
N <- fread(fn)

nmz <- paste0('newrel_',c('m04','f04','m59','f59','m1014','f1014','m014','f014'))
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
fn <- whodir + 'TB_burden_age_sex_2020-10-15.csv'
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
fn <- whodir + 'TB_burden_countries_2020-10-15.csv'
H <- fread(fn)
H <- H[year==2019,.(iso3,e_tbhiv_prct,e_tbhiv_prct_lo,e_tbhiv_prct_hi)]
H <- H[iso3 %in% cnisos]
H[,hiv:=e_tbhiv_prct/100]
H[,hiv.sd:=(e_tbhiv_prct_hi-e_tbhiv_prct_lo)/392] #adults of course
H <- H[,.(iso3,hiv,hiv.sd)]
save(H,file=here('data/H.Rdata'))


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

save(CDR,file=here('data/CDR.Rdata'))

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

## TODO questions
## HIV mix -- different under intervention