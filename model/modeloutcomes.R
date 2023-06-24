## flags for sensitivity analyses
shell <- FALSE #whether running from shell script or not
if(shell){
  ## running from shell
  args <- commandArgs(trailingOnly=TRUE)
  print(args)
  SA <- args[1]                  #none,base/lo/hi,cdr,txd
  if(SA == 'none'){
    SA <- ''
  }
} else { #set by hand
  rm(list=ls()) #clear all
  shell <- FALSE #whether running from shell script or not
  ##sensitivity analyses (mostly for PT):
  ## '' = basecase
  ## 'discr'='base'/'lo'/'hi'
  ## 'cdr' = making cdr higher for incidence
  ## 'txd' = making the completion influence tx/pt outcome
  sacases <- c('','base','cdr','txd', 'dalys')
  SA <- sacases[1]
}
SAT <- ifelse(SA=='txd','txd','') #SA relevant to Tx
ACF <- 1

## libraries
install.packages(here)
library(here)
library(data.table)
library(HEdtree)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(scales)
library(glue)

## for CEAC plotting
source(here('../dataprep/tippifunctions.R')) #CEAC & plotting utils
## modpath <- here('model')
modpath <- '' #seem to get slightly different behaviour wrt sub .here's mac vs linux
gh <- function(x) glue(here(modpath,x))

## ===== INPUT DATA
## many of these are made by modeldata.R
load(file=gh('data/edat.Rdata')) #effect data from inference NOTE using empirical atm
load(file=gh('data/LYK.Rdata'))  #LYs discounted
if(SA %in% c('hi','lo')){
  sa.drn <- ifelse(SA=='lo',0,5)
  fn <- gh('data/LYK_{sa.drn}.Rdata')
  load(fn)
}
load(file=gh('data/DBC.Rdata')) #cascade ratios for int v bl
load(file=gh('data/ATR.Rdata')) #ATT cascade
load(file=gh('data/ART2.Rdata')) #ATT cascade & costs NOTE 2 update w/costs
load(file=gh('data/CDR.Rdata')) #CDR
load(file=gh('data/H.Rdata')) #HIV by country from baseline data
load(file=gh('data/PTFH.Rdata')) #PT from HIV split
load(file=gh('data/PTC.Rdata')) #PT cascade: HH screened per PT init
load(file=gh('data/CD.Rdata'))  #rawer cost data (made from csv if modeldata.R)
load(file=gh('data/ASM.Rdata')) #age splits from pre/post data
load(file=gh('data/BC.Rdata')) #PT v ATT split from pre/post data
load(file=gh('data/HHCM.Rdata'))#HHCM cascade screen:presume:dx x FB/CB
load(file=gh('data/BL.Rdata'))#BL data: tbdx, hiv, prtb, Xpert, pt, pthiv
load(file=gh('data/INT.Rdata')) #INT cascade data from spreadsheet
load(file=gh('data/CETM.Rdata'))         #CE thresholds
load(file=gh('data/SBEP.Rdata')) #screening by entry point (new) (FB/CB HHCM vs HIV+/- ICF by country)
load(file=gh('data/PD.Rdata'))           #modelling parmeters
PZ <- parse.parmtable(PD)              #make into parm object

## --- settings
set.seed(1234)
ceactop <- 3e3 #top to plot in CEAC curves

## country key
CK <- data.table(iso3=unique(LYK$iso3)) #NOTE this is where the countries involved are coded
CK[iso3=='CMR',country:='Cameroon']
CK[iso3=='CIV',country:="Cote d'Ivoire"]
CK[iso3=='COD',country:='DRC']
CK[iso3=='KEN',country:='Kenya']
CK[iso3=='LSO',country:='Lesotho']
CK[iso3=='MWI',country:='Malawi']
CK[iso3=='TZA',country:='Tanzania']
CK[iso3=='UGA',country:='Uganda']
CK[iso3=='ZWE',country:='Zimbabwe']
CK
LYK <- merge(LYK,CK,by='iso3')

## restrict intervention resource data to countries used
keep <- c('metric',CK$country)
INT <- INT[,..keep]
INT[is.na(INT)] <- 0

## part0
## ================= base PSA data for both analyses ===================
## make PSA df for outcome data
PSA <- makePSA(max(edat$id),PZ,
               dbls = list(c('hivartOR:mn','hivartOR:sg')))
PSA[,id:=1:nrow(PSA)]
## make ontx HIV+ outcomes
PSA[,ontxHAY:=ilogit(logit(ontxY) + `hivartOR:mn` + `hivartOR:sg`)] #on ART
PSA[,ontxHAO:=ilogit(logit(ontxO) + `hivartOR:mn` + `hivartOR:sg`)] #on ART
## PSA[,ontxHAY:=pmax(ontxHAY,ontxY)]; PSA[,ontxHAO:=pmax(ontxHAO,ontxO)]

CFRdata <- PSA[,.(id,notxY,ontxY,
                  notxO,ontxO,
                  notxHAY,notxHAO,
                  ontxHAY,ontxHAO)] #,ontxHAY,ontxHAO
CFRdatam <- melt(CFRdata,id='id')
CFRdatam[,age:='5-14']; CFRdatam[grepl('Y',variable),age:='0-4']
CFRdatam[,variable:=gsub("Y$|O$","",variable)]
CFRdatam <- dcast(CFRdatam,id+age~variable,value.var = 'value')
CFRdatam[,c('dN','dA'):=.(notx-ontx,notxHA-ontxHA)] #delta-CFR
mean(CFRdatam$dA<0) #about 1% of the time random draw means -ve effect
CFRdatam[dA<0,notxHA:=ontxHA * notx/ontx] #remove by adhoc assumption
CFRdatam[,dA:=notxHA-ontxHA] #delta-CFR

## changes in ATT or PT success
nmz <- names(INT)[-1]
txsuccess <- data.table(
    country=nmz,
    BL=unlist(INT[metric=='TxSuccesspcBL',..nmz]),
    INT=unlist(INT[metric=='TxSuccesspc',..nmz])
)
ptsuccess <- data.table(
    country=nmz,
    BL=unlist(INT[metric=='PtcompletepcBL',..nmz]),
    INT=unlist(INT[metric=='Ptcompletepc',..nmz])
)

## write out
fwrite(txsuccess,file=gh('outdata/txsuccess.csv'))
fwrite(ptsuccess,file=gh('outdata/ptsuccess.csv'))


tmp <- CFRdatam[,.(age,dN,dA)]
tmp <- melt(tmp,id='age')

## NOTE see also PSA manipulation specific to PT in second part

## HHCM cascade in aggregate: activity per index
## NOTE aggregated over mode
HHCM <- HHCM[,.(value=sum(value)),by=.(age,activity)]
## make relative to index cases with HHCM
HHCM[,value:=value/HHCM[activity=="index cases with HHCM",value]]

## part1
## ================= ATT component ============================

## int/soc ratios for cascades at different stages
E1 <- merge(DBC[,.(country,metric,ratio)],CK,by='country')

## merge against cascade/cost data
K <- merge(ART2,E1,by=c('iso3','metric'),all.x = TRUE)
K[is.na(ratio),ratio:=1] #TB tx or dx
K[,country:=NULL]
K[iso3=='ZWE'] #check
##NOTE this same structure (K) is developed in model data
## would be best to ensure consistency and remove from this file

## computing costs by activity: soc -> int cascade change + increment
KA <- copy(K)
KA[,cost.soc:=frac*uc.soc/ratio] #cost assuming resource use same as soc
KA[,cost.int:=frac*(uc.soc+uc.int)] #including scale-up
KA[iso3=='MWI']                     #inspect


## formatting and plotting
KAM <- melt(KA[metric!='Diagnosed with TB',
               .(metric,iso3,cost.soc,cost.int)],
            id.vars = c('metric','iso3'))
KAM <- merge(KAM,CK,by='iso3')
KAM$metric <- factor(KAM$metric,
                     levels=c('Screened for symptoms',
                              'Presumptive TB identified',
                              'Presumptive TB tested on Xpert',
                              'TB treatment'),
                     ordered = TRUE)
KAM[grepl('soc',variable),variable:='SOC']
KAM[grepl('int',variable),variable:='Intervention']


GP <- ggplot(KAM,aes(iso3,value,fill=metric)) +
  geom_bar(stat='identity') +
  facet_wrap(~variable) +
  scale_fill_colorblind() +
  scale_y_continuous(label=comma)+
  xlab('Country') + ylab('Cost per child treated (USD)') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
if(!shell) GP

ggsave(GP,file=gh('graphs/cost_cascade.png'),h=6,w=10)
## ggsave(GP,file=gh('graphs/cost_cascade.pdf'),h=6,w=10)
KA[iso3=='LSO']                     #inspect

GP <- ggplot(KAM[variable=='Intervention'],
       aes(country,value,fill=metric)) +
  geom_bar(stat='identity',position='fill') +
  facet_wrap(~variable) +
  scale_fill_colorblind() +
  scale_y_continuous(label=percent)+
  xlab('Country') + ylab('Fraction of cost per child treated by stage') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1),
          legend.position = 'top',legend.title = element_blank())
if(!shell) GP


ggsave(GP,file=gh('graphs/cost_cascade2.png'),h=6,w=9)
## ggsave(GP,file=gh('graphs/cost_cascade2.pdf'),h=6,w=9)


## compute total costs & SD
K2 <- K[,.(cost.soc=sum(frac*uc.soc/ratio),
           cost.soc.sd=ssum(frac*uc.soc.sd/ratio),
           cost.int=sum(frac*(uc.soc+uc.int)),
           cost.int.sd=ssum(frac*uc.soc.sd)),
        by=iso3]
K2 <- merge(K2,CK,by='iso3')

## merge cost data into effect data
T <- edat[quant!='px']
T <- merge(T,K2,by='country',all.x=TRUE)

## gamma samples for cost uncertainty (with safeties for 0)
T[,costv.soc:=rgamma(nrow(T),
                  shape=(cost.soc/(cost.soc.sd+1e-6))^2,
                  scale = cost.soc.sd^2/(cost.soc+1e-6))]
T[,costv.int:=rgamma(nrow(T),
                  shape=(cost.int/(cost.int.sd+1e-6))^2,
                  scale = cost.int.sd^2/(cost.int+1e-6))]
## multiply unit costs by volume (ie # treated)
T[,c('costt.soc','costt.int'):=.(costv.soc,costv.int*RR)]

## incremental cost
T[,Dcost:=costt.int-costt.soc]

## HIV prevalence from BL data
H <- BL[country %in% CK$country,.(country,hiva=hiv,Tbdx)]
H <- H[rep(1:nrow(H),each=max(PSA$id))]
H[,id:=rep(1:max(PSA$id),nrow(CK)-1)] #NOTE missing TZA
H <- rbind(H,H)
H[,age:=rep(c('0-4','5-14'),each=nrow(H)/2)]
H[,hiv:=rbeta(nrow(H),hiva,Tbdx)]
H <- H[,.(id,country,age,hiv)]
## fill mean for TZA
tmp <- H[,.(country='Tanzania',hiv=mean(hiv)),by=.(id,age)]
H <- rbind(H,tmp)
H <- merge(H,CK,by='country')

## merge in HIV prevalence
T <- merge(T,H[,.(id,iso3,age,hiv)],
           by=c('id','age','iso3'),all.x = TRUE)

## make data on average CFR change
## merge CFR data into cost/effect data
T <- merge(T,CFRdatam,by=c('id','age')) #merge in CFR
T <- merge(T,txsuccess,by='country')    #change in tx success kk


## calculate mean differential CFR: assume all HIV on ART
T[,CFRnotx:=notx*(1-hiv) + notxHA*hiv]
T[,CFRtx:=ontx*(1-hiv) + ontxHA*hiv]
T[,CFRtx.soc:=CFRtx]
if(SA=='txd'){                        #sensitivity analysis
    T[,CFRtx.soc:=CFRtx * (1-BL)/(1-INT)] #scale up tx CFR by non-success
}


## incremental lives saved - change in implied by RR
## before - 1 ontx:(RR-1) notx
## after - RR ontx: 0     notx
T[,deaths.int:=CFRtx*RR]         #deaths in intervention
T[,deaths.soc:=CFRtx.soc + (RR-1)*CFRnotx]         #deaths in soc
T[,LYD.int:=0.25*0.331*RR]         # life years lived with disease in intervention
T[,LYD.soc:=1*0.25*0.331 + (RR-1)*0.5*0.331]         # life years lived with disease in soc
T[,LS:=deaths.soc-deaths.int]
T[,dYLD:=LYD.soc-LYD.int]
if(SA=='txd'){                        #sensitivity analysis
    T[,LS.hiv0:=(ontx*(1-BL)/(1-INT) + (RR-1)*notx)-ontx*RR]
} else {
    T[,LS.hiv0:=(ontx + (RR-1)*notx)-ontx*RR]
}



## u5/o5 split 
## S <- data.table(age=c('0-4','5-14'),frac=c(0.6,0.4))
## S <- CDR[qty=='frac',.(iso3,age,frac=value)] #frac us WHO est case mix
S <- ASM[qty=='tx']
S <- S[rep(1:nrow(S),each=max(T$id))]
S <- S[,frac:=rbeta(nrow(S),`0-4`,`5-14`)] #approx flat conjugate prior
S[,id:=rep(1:max(T$id),nrow(S)/max(T$id))]
S <- dcast(S[,.(country,id,period,frac)],
           country+id~period,value.var = 'frac')
T <- merge(T,S,by=c('country','id'),all.x=TRUE)
T[age=='5-14',Baseline:=1-Baseline] #defined as prop u5
T[age=='5-14',Intervention:=1-Intervention] #defined as prop u5
T[iso3=='ZWE' & id==1] #check
names(T)[names(T)=='Baseline'] <- 'frac'
names(T)[names(T)=='Intervention'] <- 'fracI'

## merge in life-expectancy & calculate DALY changes
T <- merge(T,LYK[,.(iso3,age,LYS,LYS0)],by=c('iso3','age'),all.x=TRUE)
T[,c('dLYS','dLYS0','dLYS.nohiv'):=.(LYS*LS,LYS0*LS,LYS*LS.hiv0)]

if(SA=='dalys') {
  T[,c('dDALY','dDALY0','dDALY.nohiv'):=.(dLYS,dLYS0,dLYS.nohiv)]
} else 
{
  # T[,dYLD:=0.25*0.331]         # saved years lived with disease in intervention
  T[,c('dDALY','dDALY0','dDALY.nohiv'):=.(dLYS+dYLD,dLYS0+dYLD,dLYS.nohiv+dYLD)]
}

## calculate differential costs & DALYs over ages
## NOTE frac here is age split
T1 <- T[,.(cost.soc=sum(costt.soc*frac),
           cost.int=sum(costt.int*frac),
           Dcost=sum(Dcost*frac),
           LS=sum(LS*frac),
           ## deaths.soc=sum((deaths.int+LS)*frac),
           tx=sum(RR*frac),
           deaths.soc=sum(deaths.soc*frac),
           deaths.int=sum(deaths.int*frac),
           # dYLD=sum(dYLD*frac),
           dDALY0=sum(dDALY0*frac),
           dDALY=sum(dDALY*frac),
           dDALY.nohiv=sum(dDALY.nohiv*frac)),
        by=.(country,id)]
## T[,mean(RR-1),by=.(country,age)]
T1 <- merge(T1,CK,by='country') #country iso3 merged on

## same with age - basically just consistent renaming
T2 <- T[,.(cost.soc=(costt.soc),
           cost.int=(costt.int),
           Dcost=(Dcost),
           LS=(LS),
           tx = RR,
           ## deaths.soc=((deaths.int+LS)),
           deaths.soc=(deaths.soc),
           deaths.int=(deaths.int),
           # dYLD=(dYLD),
           dDALY0=(dDALY0),
           dDALY=(dDALY),
           dDALY.nohiv=(dDALY.nohiv)),
        by=.(country,id,age)]
T2 <- merge(T2,CK,by='country') #country iso3 merged on

## inspect
summary(T1)

T1$country <- gsub("ote","ôte",T1$country)
CETM$country <- gsub("ote","ôte",CETM$country)

## --- CEA and CEAC plots ---

## CEA plot
GP <- ggplot(T1,aes(dDALY,Dcost)) +
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    geom_point(alpha=0.1,shape=1) +
    geom_abline(data=CETM[threshold %in% c("0.5x GDP","1x GDP")],
                aes(intercept=0,slope=value,col=threshold))+
    facet_wrap(~country) +
    scale_y_continuous(label=comma) +
    xlab('Incremental discounted life-years saved')+
    ylab('Incremental cost (USD)')+
    theme(legend.position = "top" )
if(!shell) GP

fn1 <- gh('graphs/CEall') + SAT + '.png'
## fn2 <- gh('graphs/CEall') + SAT + '.pdf' 
ggsave(GP,file=fn1,w=10,h=10); ## ggsave(GP,file=fn2,w=10,h=10)
## PDF versions too big

## make CEAC data
lz <- seq(from = 0,to=ceactop,length.out = 1000)
ceacd <- list()
for(cn in unique(T1$country)){
    tmp <- T1[country==cn,.(Q=dDALY,P=Dcost)]
    ceacd[[cn]] <- data.table(
        country=cn,x=lz,
        y=make.ceac(tmp,lz))
}
ceacd <- rbindlist(ceacd)

## make CEAC plot
CEAC <- make.ceac.plot(ceacd,xpad=50) +
  xlab('Cost-effectiveness threshold (USD per DALY averted)')
if(!shell) CEAC

fn1 <- gh('graphs/CEAC') + SAT + '.png'
## fn2 <- gh('graphs/CEAC') + SAT + '.pdf'
ggsave(CEAC,file=fn1,w=7,h=7); ## ggsave(CEAC,file=fn2,w=7,h=7)

# labelled CEAC
CEACL <- ceacd %>%
  mutate(label_hjust = rep(c(0.85, 0.85, 0.25, 0.4, 0.95,0.25, 0.25, 0.25, 0.85), each=nrow(ceacd)/length(unique(country)))) %>%
  ggplot() +
  geomtextpath::geom_textline(aes(x,y,
                                  col=country,
                                  label=country,
                                  hjust = label_hjust),
                              size=4, vjust = 0.3, show.legend = FALSE, text_only = TRUE, text_smoothing = 40) + 
  theme_classic() +
  guides(col = guide_legend(order = 1, nrow = 2, title = 'Country')) +
  theme(legend.position = 'top',legend.title = element_blank())+
  # theme(legend.text=element_text(size=rel(0.8)))+
  ggpubr::grids()+
  scale_x_continuous(label=comma,
                     breaks = seq(0,ceactop,250)) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73",
                               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                               "darkorchid1"))+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD per DALY averted)')

CEACL
ggsave(CEACL,file=gh('graphs/CEACL') + SA + '.png',w=7,h=7)

## output where things X 50%
tmp <- ceacd[y>0.5 & abs(y-0.5)<5e-2]
tmp[,err:=abs(y-0.5)]
tmp[,ermin:=min(err),by=country]
tmp <- tmp[err==ermin]
tmp <- tmp[,.(country,ceac50=round(x))]
tmp

fn1 <- gh('outdata/CEAC50') + SA + '.csv'
fwrite(tmp,file=fn1)

## ICERs by country
ice <- T1[,
  .(cost.soc=mean(cost.soc), #costs
    cost.soc.lo=lof(cost.soc),
    cost.soc.hi=hif(cost.soc),
    cost.int=mean(cost.int),
    cost.int.lo=lof(cost.int),
    cost.int.hi=hif(cost.int),
    Dcost=mean(Dcost),
    Dcost.lo=lof(Dcost),
    Dcost.hi=hif(Dcost),
    tx=mean(tx),
    tx.lo=lof(tx),
    tx.hi=hif(tx),
    ## deaths
    deaths.soc=mean(deaths.soc),
    deaths.soc.lo=lof(deaths.soc),
    deaths.soc.hi=hif(deaths.soc),
    deaths.int=mean(deaths.int),
    deaths.int.lo=lof(deaths.int),
    deaths.int.hi=hif(deaths.int),
    Ddeaths=mean(deaths.int-deaths.soc),
    Ddeaths.lo=lof(deaths.int-deaths.soc),
    Ddeaths.hi=hif(deaths.int-deaths.soc),
    LS=mean(LS),LS.lo=lof(LS),LS.hi=hif(LS),
    ## dalys
    dDALY0=mean(dDALY0),
    dDALY=mean(dDALY),
    dDALY.nohiv=mean(dDALY.nohiv),
    dDALY.hi=hif(dDALY),
    dDALY.lo=lof(dDALY),
    dDALY0.hi=hif(dDALY0),
    dDALY0.lo=lof(dDALY0),
    ## ICER
    ICER=mean(Dcost)/mean(dDALY)
    ),
  by=country]

## multiply most by 100
vec <- names(ice)
vec <- vec[!vec %in% c('country','ICER')] #all but country and ICER
ice[,(vec):=lapply(.SD, function(x) 100*x), .SDcols = vec]

icer <- ice[,.(country=country,
               treated.soc=100,         #NOTE normalized to 100
               cost.soc = bracket(cost.soc,cost.soc.lo,cost.soc.hi),
               treated.int=bracket(tx,tx.lo,tx.hi),
               cost.int = bracket(cost.int,cost.int.lo,cost.int.hi),
               treated.dif=bracket(tx-1e2,tx.lo-1e2,tx.hi-1e2), #NOTE also 100
               cost.dif=bracket(Dcost,Dcost.lo,Dcost.hi),
               deaths.dif=bracket(-LS,-LS.hi,-LS.lo),
               LY0.dif=bracket(dDALY0,dDALY0.lo,dDALY0.hi),
               LY.dif=bracket(dDALY,dDALY.lo,dDALY.hi),
               ICER = format(round(ICER),big.mark = ',')
               )]
icer

fn1 <- gh('outdata/ICERatt') + SA + '.csv'
fwrite(icer,file=fn1)

## --table 2 format
Table2ATT <- icer[,.(country,treated.soc,cost.soc,treated.int,cost.int,
                     diff.PT=0,diff.incTB=0,
                     diff.ATT=treated.dif,
                     diff.deaths=deaths.dif,
                     diff.LYS=LY0.dif,diff.dLYS=LY.dif,
                     diff.cost=cost.dif,ICER)]

fn1 <- gh('outdata/Table2ATT') + SAT + '.Rdata'
save(Table2ATT,file=fn1)

T2$country <- gsub("ote","ôte",T2$country)

## --- ICER tables  by age (as above)
## ICERs by country & age
iceage <- T2[,
             .(cost.soc=mean(cost.soc), #costs
             cost.soc.lo=lof(cost.soc),
             cost.soc.hi=hif(cost.soc),
             cost.int=mean(cost.int),
             cost.int.lo=lof(cost.int),
             cost.int.hi=hif(cost.int),
             Dcost=mean(Dcost),
             Dcost.lo=lof(Dcost),
             Dcost.hi=hif(Dcost),
             tx=mean(tx),
             tx.lo=lof(tx),
             tx.hi=hif(tx),
             ## deaths
             deaths.soc=mean(deaths.soc),
             deaths.soc.lo=lof(deaths.soc),
             deaths.soc.hi=hif(deaths.soc),
             deaths.int=mean(deaths.int),
             deaths.int.lo=lof(deaths.int),
             deaths.int.hi=hif(deaths.int),
             Ddeaths=mean(deaths.int-deaths.soc),
             Ddeaths.lo=lof(deaths.int-deaths.soc),
             Ddeaths.hi=hif(deaths.int-deaths.soc),
             LS=mean(LS),LS.lo=lof(LS),LS.hi=hif(LS),
             ## dalys
             dDALY0=mean(dDALY0),
             dDALY=mean(dDALY),
             dDALY.nohiv=mean(dDALY.nohiv),
             dDALY.hi=hif(dDALY),
             dDALY.lo=lof(dDALY),
             dDALY0.hi=hif(dDALY0),
             dDALY0.lo=lof(dDALY0),
             ## ICER
             ICER=mean(Dcost)/mean(dDALY)
             ),
          by=.(country,age)]

## multiply most by 100
iceage[,(vec):=lapply(.SD, function(x) 100*x), .SDcols = vec]

icers <- iceage[,.(country=country,age,
                   treated.soc=100,         #NOTE normalized to 100
                   cost.soc = bracket(cost.soc,cost.soc.lo,cost.soc.hi),
                   treated.int=bracket(tx,tx.lo,tx.hi),
                   cost.int = bracket(cost.int,cost.int.lo,cost.int.hi),
                   treated.dif=bracket(tx-1e2,tx.lo-1e2,tx.hi-1e2), #NOTE also 100
                   cost.dif=bracket(Dcost,Dcost.lo,Dcost.hi),
                   deaths.dif=bracket(-LS,-LS.hi,-LS.lo),
                   LY0.dif=bracket(dDALY0,dDALY0.lo,dDALY0.hi),
                   LY.dif=bracket(dDALY,dDALY.lo,dDALY.hi),
                   ICER = format(round(ICER),big.mark = ',')
               )]
icers #inspect

icers[age=='0-4',.(country,treated.int,ICER)] #inspect

fn1 <- gh('outdata/ICERSatt') + SA + '.csv'
fwrite(icers,file=fn1)

## ---- drivers by variable

## frac, RR, ratios, LE
## X <- S[age=='0-4',.(iso3,frac)]         #frac
## X <- merge(X,CK,by='iso3')
X <- ASM[qty=='tx' & period=='Baseline',
         .(country,frac=`0-4`/(`0-4`+`5-14`))]
tmp <- T[,.(RR=mean(RR)),by=.(country,age)]
tmp <- dcast(tmp,country~age)
names(tmp)[2:3] <- paste0('RR',names(tmp)[2:3])
X <- merge(X,tmp,by='country')          #RR
tmp <- unique(LYK[,.(iso3,age,LYS)])
tmp <- dcast(tmp,iso3~age,value.var = 'LYS')
names(tmp)[2:3] <- paste0('LYS',names(tmp)[2:3])
tmp <- merge(tmp,CK,by='iso3')
X <- merge(X,tmp,by='country')          #LYS
tmp <- dcast(KA[,.(iso3,metric,ratio,frac)],
             iso3~metric,value.var = c('ratio','frac'))
X <- merge(X,tmp,by='iso3')          #cascade data

## make outcome and covariate data
Y <- ice[,.(country,ICER,Dcost,dDALY)]
XM <- melt(X,id=c('iso3','country'))
XYc <- merge(XM,Y,by='country',all.x=TRUE) #merge together
XYc$country <- gsub("ote","ôte",XYc$country)

## make and save plots

## cost
GP <- ggplot(XYc,aes(value,Dcost,col=iso3)) +
  geom_point() +
  facet_wrap(~variable,scales = 'free')
if(!shell) GP

ggsave(GP,file=gh('graphs/drivers_att_cost2.png'),h=10,w=10)
## ggsave(GP,file=gh('graphs/drivers_att_cost2.pdf'),h=10,w=10)

## DALYs
GP <- ggplot(XYc,aes(value,dDALY,col=iso3)) +
  geom_point() +
  facet_wrap(~variable,scales = 'free')
if(!shell) GP

ggsave(GP,file=gh('graphs/drivers_att_DALY.png'),h=10,w=10)
## ggsave(GP,file=gh('graphs/drivers_att_DALY.pdf'),h=10,w=10)


## make CEACY (0-4 years) data
lz <- seq(from = 0,to=ceactop,length.out = 1000)
ceacdy <- list()
for(cn in unique(T2$country)){
  tmp <- T2[age=='0-4',]
  tmp <- tmp[country==cn,.(Q=dDALY,P=Dcost)]
  ceacdy[[cn]] <- data.table(
    country=cn,x=lz,
    y=make.ceac(tmp,lz))
}
ceacdy <- rbindlist(ceacdy)

## make CEACO plot
CEACY <- make.ceac.plot(ceacdy,xpad=50)  +
  xlab('Cost-effectiveness threshold (USD per DALY averted)') +
  theme(legend.position="top", legend.title = element_blank())
if(!shell) CEACY

fn1 <- glue(here('graphs/CEACY')) + SAT + '.png'
## fn2 <- glue(here('graphs/CEACY')) + SAT + '.pdf'
ggsave(CEACY,file=fn1,w=7,h=7); ## ggsave(CEAC,file=fn2,w=7,h=7)

# labelled CEAC
CEACYL <- ceacdy %>%
  mutate(label_hjust = rep(c(0.85, 0.85, 0.25, 0.4, 0.6,0.4, 0.25, 0.4, 0.85), each=nrow(ceacd)/length(unique(country)))) %>%
  ggplot() +
  geomtextpath::geom_textline(aes(x,y,
                                  col=country,
                                  label=country,
                                  hjust = label_hjust),
                              size=4, vjust = 0.3, show.legend = FALSE, text_only = TRUE, text_smoothing = 40) + 
  theme_classic() +
  guides(col = guide_legend(order = 1, nrow = 2, title = 'Country')) +
  theme(legend.position = 'top',legend.title = element_blank())+
  # theme(legend.text=element_text(size=rel(0.8)))+
  ggpubr::grids()+
  scale_x_continuous(label=comma,
                     breaks = seq(0,ceactop,250)) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73",
                               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                               "darkorchid1"))+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD per DALY averted)')

CEACYL
ggsave(CEACYL,file=gh('graphs/CEACYL') + SA + '.png',w=7,h=7)

## make CEACO (5-14 years) data
lz <- seq(from = 0,to=ceactop,length.out = 1000)
ceacdo <- list()
for(cn in unique(T2$country)){
  tmp <- T2[age=='5-14',]
  tmp <- tmp[country==cn,.(Q=dDALY,P=Dcost)]
  ceacdo[[cn]] <- data.table(
    country=cn,x=lz,
    y=make.ceac(tmp,lz))
}
ceacdo <- rbindlist(ceacdo)

## make CEACO plot
CEACO <- make.ceac.plot(ceacdo,xpad=50)  +
  xlab('Cost-effectiveness threshold (USD per DALY averted)') +
  theme(legend.position="top", legend.title = element_blank())
if(!shell) CEACO

fn1 <- glue(here('graphs/CEACO')) + SAT + '.png'
## fn2 <- glue(here('graphs/CEACO')) + SAT + '.pdf'
ggsave(CEACO,file=fn1,w=7,h=7); ## ggsave(CEACO,file=fn2,w=7,h=7)

# labelled CEAC
CEACOL <- ceacdo %>%
  mutate(label_hjust = rep(c(0.85, 0.98, 0.25, 0.4, 0.95,0.25, 0.25, 0.25, 0.95), each=nrow(ceacd)/length(unique(country)))) %>%
  ggplot() +
  geomtextpath::geom_textline(aes(x,y,
                                  col=country,
                                  label=country,
                                  hjust = label_hjust),
                              size=4, vjust = 0.3, show.legend = FALSE, text_only = TRUE, text_smoothing = 40) + 
  theme_classic() +
  guides(col = guide_legend(order = 1, nrow = 2, title = 'Country')) +
  theme(legend.position = 'top',legend.title = element_blank())+
  # theme(legend.text=element_text(size=rel(0.8)))+
  ggpubr::grids()+
  scale_x_continuous(label=comma,
                     breaks = seq(0,ceactop,250)) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73",
                               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                               "darkorchid1"))+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD per DALY averted)')

CEACOL
ggsave(CEACOL,file=gh('graphs/CEACOL') + SA + '.png',w=7,h=7)

## part2
## ================= PT component ============================
## needs edat, PSA, INT
PT <- edat[quant=='px']

## format relevant PSA data by age
## NOTE check whether HHhivprev could be removed
PSAage <- PSA[,.(id,prog04,prog514,HHhivprev04,HHhivprev514,
                 LTBI04,LTBI514)]
PSAage <- melt(PSAage,id='id')
PSAage[,age:='5-14']
PSAage[grepl('04',variable),age:='0-4']
PSAage[,variable:=gsub('04','',variable)]
PSAage[,variable:=gsub('514','',variable)]
PSAage <- dcast(PSAage,id+age~variable)

## relevant PSA data with no age dependence
PSApt <- PSA[,.(id,hivpi,artp,iptRRtstpos,iptRRhivpos)]
PSApt <- merge(PSAage,PSApt,by='id')

## merge parameters onto effect data
PT <- merge(PT,PSApt,by=c('id','age')) #progression/PT etc
PT <- merge(PT,CFRdatam,by=c('id','age')) #CFRs

## merge in LEs
PT <- merge(PT,LYK,by=c('country','age'),all.x=TRUE)

## make PSA for country CDRs
CDRs <- CDR[qty=='cdr']
CDRs <- CDRs[rep(1:nrow(CDRs),each=max(PSA$id))]
CDRs[,id:=rep(1:max(PSA$id),nrow(CDRs)/max(PSA$id))]
CDRs[,sz:=value*(1-value)/cdr.v-1]
CDRs[,c('a','b'):=.(sz*value,sz*(1-value))]
CDRs[,cdr0:=rbeta(nrow(CDRs),a,b)]
CDRs[,cdr:=cdr0] #basecase
if(SA=='cdr'){   #sensitivity analysis
    CDRs[,cdr:=runif(nrow(CDRs))]
    CDRs[,cdr:=cdr0*(1-cdr) + cdr] #interpolate between cdr0 & 1
}
CDRs[,summary(cdr)]

## merge in CDRs
PT <- merge(PT,CDRs[,.(iso3,age,id,cdr)],by=c('iso3','age','id'))


## age & HIV-route splits for PT
## NOTE uncertainty probably not necessary due to large numbers
## age splits
tmp <- INT[metric %in% c("PThhcu5pc","PTHIVentryu5pc")]
tmp <- melt(tmp,id='metric')
tmp <- dcast(tmp[,.(metric,country=variable,value)],
             country~metric,value='pc')
tmp <- tmp[country %in% CK$country]
hag <- merge(tmp,PTFH[,.(country,ptinhiv)],by='country') #both splits
hag[,heu5:=ptinhiv*PTHIVentryu5pc];hag[,heo5:=ptinhiv*(1-PTHIVentryu5pc)]
hag[,hcu5:=(1-ptinhiv)*PThhcu5pc];hag[,hco5:=(1-ptinhiv)*(1-PThhcu5pc)]
hag[,heu5+heo5+hcu5+hco5]               #check
## combined age/entry-point splits in right format
hag <- melt(hag[,.(country,heu5,heo5,hcu5,hco5)],id='country')
hag[,age:='0-4']; hag[grepl('o5',variable),age:='5-14'] #age var
hag[,route:='hhc']; hag[grepl('he',variable),route:='hentry']
hag <- dcast(hag[,.(country,value,route,age)],
             country+age~route,value.var='value')

## outputting cascade data
hago <- dcast(hag,country~age,value.var = c('hentry','hhc'))
hago <- merge(hago,PTC,by='country')
vec <- names(hago)[2:5]
hago[,(vec):=lapply(.SD, function(x) 100*x), .SDcols = vec]
vec <- names(hago)[2:6]
hago[,(vec):=lapply(.SD, function(x) round(x,2)), .SDcols = vec]

fwrite(hago,file=gh('outdata/PTC.csv')) #save as output too

## reformatting into parts of Table 1?


## merge into PT data
PT <- merge(PT,hag,by=c('country','age'),all.x = TRUE)
## NOTE this is normalized over route x age - age-stratified results will need renorm'n:
PT[,ptentry:=hentry+hhc]  #denominator for age-stratified calx

## ARI assumption for risks in HIV entry-point cohort
PT <- merge(PT,PSA[,.(id,ari)],by='id',all.x=TRUE)

## === outcomes
## NOTE now normalizing by age, for consistency with ATT

## NOTE include split to hep = HIV-clinic entry-point
## --- no PT outcomes on average
## cases - stratified by HIV for CFR
PT[,casesnoPT.hiv.hep:=
  (hentry/ptentry) * (               #HIV-entry
    ari*prog*                        #prog - bl
    hivpi*artp                       #HIV modification
  )
]
PT[,casesnoPT.hiv:=(hhc/ptentry)*(   #HH route
  (LTBI*prog)*                       #progn - bl
  (HHhivprev*hivpi*artp)             #HIV modification
) + casesnoPT.hiv.hep
]
PT[,casesnoPT.nohiv:=(hhc/ptentry)*(   #HH route
  (LTBI*prog)*                         #progn - bl
  (1-HHhivprev)                        #HIV modification
)
]
PT[,casesnoPT:=casesnoPT.hiv + casesnoPT.nohiv]

## treatment
PT[,attnoPT:=(casesnoPT.hiv+casesnoPT.nohiv) * cdr]
PT[,attnoPT.hep:=casesnoPT.hiv.hep * cdr]

## deaths
PT[,deathsnoPT:=(
  casesnoPT.hiv*(cdr*ontxHA+(1-cdr)*notxHA)+
  casesnoPT.nohiv*(cdr*ontx+(1-cdr)*notx)
)]
PT[,deathsnoPT.hep:=casesnoPT.hiv.hep*(cdr*ontxHA+(1-cdr)*notxHA)]


## ---  PT outcomes on average
PT <- merge(PT,ptsuccess,by='country')

## cases

## possibly account for changes in incomplete treatment
PT[,c('fs','fi'):=1.0] #proportion completing
if(SA=='txd'){ #sensitivity analysis around completion
    PT[,c('fs','fi'):=.(BL,INT)]
}

## NOTE these are separated SOC/INT to allow SA around completion improvement
## SOC
PT[,casesPT.hiv.soc.hep:=
  (hentry/ptentry) * (               #HIV-entry
    ari*prog*                        #prog - bl
    hivpi*artp                       #HIV modification
  )*(iptRRhivpos*fs+1-fs)                   #PT
]
PT[,casesPT.hiv.soc:=(hhc/ptentry)*(   #HH route
    (LTBI*prog)*         #bl progn 
    (HHhivprev*hivpi*artp)* #HIV/PT modification
    (iptRRtstpos*fs + 1-fs)  #PT modification
) + casesPT.hiv.soc.hep      #HEP
]
PT[,casesPT.nohiv.soc:=(hhc/ptentry)*(   #HH route
    (LTBI*prog)*                         #progn - bl
    (1-HHhivprev)*(iptRRtstpos*fs+1-fs)           #PT
)
]
PT[,casesPT.soc:=casesPT.hiv.soc + casesPT.nohiv.soc]

## treatment
PT[,attPT.soc:=(casesPT.hiv.soc+casesPT.nohiv.soc) * cdr]
PT[,attPT.soc.hep:=casesPT.hiv.soc.hep * cdr]

## deaths
PT[,deathsPT.soc:=(
  casesPT.hiv.soc*(cdr*ontxHA+(1-cdr)*notxHA)+
  casesPT.nohiv.soc*(cdr*ontx+(1-cdr)*notx)
)]
PT[,deathsPT.soc.hep:=casesPT.hiv.soc.hep*(cdr*ontxHA+(1-cdr)*notxHA)]

## INT
PT[,casesPT.hiv.int.hep:=
  (hentry/ptentry) * (               #HIV-entry
    ari*prog*                        #prog - bl
    hivpi*artp                       #HIV modification
  )*(iptRRhivpos*fi+1-fi)                        #PT
]
PT[,casesPT.hiv.int:=(hhc/ptentry)*(   #HH route
    (LTBI*prog)*         #bl progn 
    (HHhivprev*hivpi*artp)* ## HIV modification
    (iptRRtstpos*fi+1-fi)  #PT modification
) +
  casesPT.hiv.int.hep
]
PT[,casesPT.nohiv.int:=(hhc/ptentry)*(   #HH route
    (LTBI*prog)*                         #progn - bl
    (1-HHhivprev)*(iptRRtstpos*fi+1-fi)           #PT
)
]
PT[,casesPT.int:=casesPT.hiv.int + casesPT.nohiv.int]

## treatment
PT[,attPT.int:=(casesPT.hiv.int+casesPT.nohiv.int) * cdr]
PT[,attPT.int.hep:=casesPT.hiv.int.hep * cdr]

## deaths
PT[,deathsPT.int:=(
    casesPT.hiv.int*(cdr*ontxHA+(1-cdr)*notxHA)+
    casesPT.nohiv.int*(cdr*ontx+(1-cdr)*notx)
)]
PT[,deathsPT.int.hep:=casesPT.hiv.int.hep*(cdr*ontxHA+(1-cdr)*notxHA)]


##  --- costs
## merge in traced HHs per PT initiation
PT <- merge(PT,PTC,by='country',all.x=TRUE) #cascade
## calculate costs NOTE normalized for each row, not over ages

## turn cost data into PSA
nr <- nrow(CD)
CDlong <- CD[rep(1:nr,each=max(PT$id))] #extending to no replicates
CDlong[,id:=rep(1:max(PT$id),nr)]
CDlong[is.na(CDlong)] <- 0 #setting NA to 0
CDlong[,socu:=rgamma(nrow(CDlong),shape=(uc.soc/(uc.soc.sd+1e-6))^2,
                  scale=uc.soc.sd^2/(uc.soc+1e-6))]
CDlong[,intu:=socu + uc.int] #NOTE intu aren't incremental now
CDlong <- CDlong[,.(iso3,Activity,id,socu,intu)] #restrict
CDlong <- dcast(CDlong,iso3+id~Activity,value.var=c('socu','intu'))#shape
names(CDlong)

## merges cost data into activity data
PT <- merge(PT,CDlong,by=c('iso3','id'),all.x=TRUE)    #costs

## the split between facility-based vs community-based HH screening
## NOTE is this the same under soc/int? presume mainly int
SBEP[,.(iso3,CBhhcm,FBhhcm)] # raw numbers
CvF <- SBEP[rep(1:nrow(SBEP),each=max(PT$id)),.(iso3,CBhhcm,FBhhcm)]
CvF[,id:=rep(1:max(PT$id),nrow(SBEP))]
CvF[,propFB:=rbeta(nrow(CvF),shape1=FBhhcm,shape2=CBhhcm)] #use numbers in beta dist
PT <- merge(PT,CvF,by=c('iso3','id'))                      #merge in

## traceperhhcpt - is households screened per PT init
## (check same as hhci)
## NOTE unit costs are per HH
PT[,costPT.soc:= (
      (hhc/ptentry) * (
        (1-propFB)*traceperhhcpt*`socu_Community hhci`+   #comm CT
        (propFB)*traceperhhcpt*`socu_Facility hhci`       #facility CT
        ) +
      (hentry/ptentry) * `socu_Screening in HIV clinic`+  #HEP
      `socu_TPT treatment`                                #TPT
 )]

PT[,costPT.int:=
      (hhc/ptentry) * (
        (1-propFB)*traceperhhcpt*`intu_Community hhci`+   #comm CT
        (propFB)*traceperhhcpt*`intu_Facility hhci`       #facility CT
        ) +
      (hentry/ptentry) * `intu_Screening in HIV clinic`+  #HEP
      `intu_TPT treatment`                                #TPT
   ]


## --- ACF here
## should be getting ~6 per 100 index cases across ages
## NOTE be careful to get right numbers across ages:
## for a u5 PT -> traceperhhcpt -> u5 & o5 screens etc
## for a o5 PT -> traceperhhcpt -> u5 & o5 screens etc
## => NOTE need to sum & attach to both u5 & o5 before aggregating
## complication that age in data structure is about PT, need to add
## across ages for coprev stuff & attach to each row/PT-age
## needs care to sum before applying RR - effect pertains to PT-age
## not contact-age
## ptentry gives the proportion of PT for each age group
PT[,.(sum(ptentry),sum(hhc)),by=.(id,country)]

## total index cases, summing over ages:
## PT[,totindexscreen := sum(ptentry*traceperhhcpt),by=.(id,country)]
## include here fraction of PT that is HH route
PT[,totindexscreen := traceperhhcpt*(hhc/ptentry)] #think this correct - weighted @end

## merge in HHCM cascade in right form (this splits out by age again)
## this is where numbers screened come in (aggregate across countries)
hhcascade <- dcast(HHCM[age %in% c('0-4','5-14')],age~activity)
PT <- merge(PT,hhcascade,by='age',all.x=TRUE)


## --- outcomes on coprevalence and deaths

## coprevalent detectable children
PT[,coprev:=totindexscreen*Diagnosed] #numbers coprevalent
PT[,coprev:=sum(coprev),by=.(id,country)] #sum over ages - per PT
PT[,mean(coprev)]                         #NOTE includes non-hhc PT

## deaths RR included later
PT[,deathsPrev.soc:= totindexscreen*Diagnosed*(
    HHhivprev*(cdr*ontxHA+(1-cdr)*notxHA)+
    (1-HHhivprev)*(cdr*ontx+(1-cdr)*notx)
)] #using relevant background CDR
PT[,deathsPrev.int:= totindexscreen*Diagnosed*(
    HHhivprev*(1.0*ontxHA+(1-1.0)*notxHA)+
    (1-HHhivprev)*(1.0*ontx+(1-1.0)*notx)
)] #using 100% CDR
PT[,c('deathsPrev.soc','deathsPrev.int'):=.(sum(deathsPrev.soc),
                                            sum(deathsPrev.int)),
   by=.(id,country)] #sum over ages - per PT


## check names
grep('socu',names(PT),value=TRUE)

## cost (not including thct)
names(PT)


## --- costs for ACF component
## NOTE
PT[,costACF.soc.hh:= ## coprev HHCM'd
      totindexscreen*(Presumptive*`socu_Presumptive TB evaluation` +
                      Diagnosed*`socu_TB treatment`)]
PT[,costACF.soc.nhh:=
      ## coprev not HHCM'd - not including PHC screening cost
      (totindexscreen*cdr*(Presumptive*`socu_Presumptive TB evaluation` +
                           Diagnosed*`socu_TB treatment`))]
PT[,costACF.int.hh:= ## coprev HHCM'd
      totindexscreen*(Presumptive*`intu_Presumptive TB evaluation`+
                      Diagnosed*`socu_TB treatment`)]
PT[,c('costACF.soc.hh','costACF.soc.nhh','costACF.int'):=
        .(sum(costACF.soc.hh),
          sum(costACF.soc.nhh),
          sum(costACF.int.hh)),
   by=.(id,country)] #sum over ages - per PT


## keeping count of ATTs
PT[,ATTACF.hh:=totindexscreen*(Diagnosed)]
PT[,ATTACF.nhh:=totindexscreen*(Diagnosed)*cdr]
PT[,c('ATTACF.hh','ATTACF.nhh'):=.(sum(ATTACF.hh),
                                   sum(ATTACF.nhh)),
   by=.(id,country)] #sum over ages - per PT

## --- NOTE now applying RRs to aggregated totals to include
## both contact-age groups
PT[,costACF.soc:= 1*costACF.soc.hh + (RR-1)*costACF.soc.nhh]
PT[,costACF.int:= RR*costACF.int.hh + 0*0]
PT[,ATTACF.soc:= 1*ATTACF.hh + (RR-1)*ATTACF.nhh]
PT[,ATTACF.int:= RR*ATTACF.hh + 0*ATTACF.nhh]
## deaths - how much of which condition (no age sum as already done)
PT[,deathsACF.soc:=1*deathsPrev.int + (RR-1)*deathsPrev.soc]
PT[,deathsACF.int:=RR*deathsPrev.int + 0*deathsPrev.soc]
PT[,LSACF:=-(deathsACF.int-deathsACF.soc)]
PT[,LYDACF.int:=0.25*0.331*RR*ATTACF.hh]         # life years lived with disease in intervention
PT[,LYDACF.soc:=1*ATTACF.hh*0.25*0.331 + (RR-1)*ATTACF.nhh*0.5*0.331]         # life years lived with disease in soc
PT[,dYLDACF:=LYDACF.soc-LYDACF.int]

## add in DALYs etc
PT[,c('dDALYacf','dDALY0acf'):=.(LYS*LSACF,LYS0*LSACF)]

if(SA=='dalys') {
  PT[,c('dDALYacf','dDALY0acf'):=.(dDALYacf,dDALY0acf)]
} else 
{
  PT[,c('dDALYacf','dDALY0acf'):=.(dDALYacf+dYLDACF,dDALY0acf+dYLDACF)]
}

## NOTE assumption of same age split under SOC - tot pop RR
## outcomes:
## cost
PT[,cost.soc:=(1*costPT.soc+(RR-1)*0) + #PT cost, includes hhc or tbe
      (1*attPT.soc+(RR-1)*attnoPT*`socu_TB treatment`)] #ATT cost
PT[,cost.int:=(RR*costPT.int+0*0) + #PT cost, includes hhc or tbe
      (RR*attPT.int+0*attnoPT*`intu_TB treatment`)] #ATT cost
## cases
PT[,cases.soc:=1*casesPT.soc + (RR-1)*casesnoPT]
PT[,cases.int:=RR*casesPT.int + 0*casesnoPT]


## ATT
PT[,ATT.soc:=cases.soc*cdr]
PT[,ATT.int:=cases.int*cdr]

## life years lived with disease
PT[,YLD.soc:=cases.soc*0.25*0.331*cdr  + cases.soc*(1-cdr)*0.5*0.331]
PT[,YLD.int:=cases.int*0.25*0.331*cdr  + cases.int*(1-cdr)*0.5*0.331]
PT[,dYLD:=YLD.soc-YLD.int]

## deaths
PT[,deaths.soc:=1*deathsPT.soc + (RR-1)*deathsnoPT]
PT[,deaths.int:=RR*deathsPT.int + 0*deathsnoPT]
PT[,LS:=-(deaths.int-deaths.soc)]
## add in DALYs etc
PT[,c('dDALY','dDALY0'):=.(LYS*LS,LYS0*LS)]

## 
if(SA=='dalys') {
  PT[,c('dDALY','dDALY0'):=.(dDALY,dDALY0)]
} else{

  PT[,c('dDALY','dDALY0'):=.(dDALY+dYLD,dDALY0+dYLD)]
}

## looking at HEP separately
PT[,cost.soc.hep:=(`socu_Screening in HIV clinic`+  #HEP
                   `socu_TPT treatment`)+
      (1*attPT.soc.hep+(RR-1)*attnoPT.hep*`socu_TB treatment`)
   ]
PT[,cost.int.hep:= RR*(
  `intu_Screening in HIV clinic`+  #HEP
  `intu_TPT treatment`                                #TPT
) + #PT cost, includes hhc or tbe
  (RR*attPT.int.hep+0*attnoPT.hep*`intu_TB treatment`) #ATT cost
]
PT[,deaths.soc.hep:=1*deathsPT.soc.hep + (RR-1)*deathsnoPT.hep]
PT[,deaths.int.hep:=RR*deathsPT.int.hep + 0*deathsnoPT.hep]
PT[,LS.hep:=-(deaths.int.hep-deaths.soc.hep)]
## add in DALYs etc
PT[,c('dDALY.hep','dDALY0.hep'):=.(LYS*LS.hep,LYS0*LS.hep)]


## ACF cascade
## denominators haven't been summed over PT-age, hhc correct split
## coprev weighted by hh/ptentry - needs ptentry to weight PT-age
tmp <- PT[,.(cppt=(sum(coprev*ptentry)/sum(hhc)), #per HH PT initiation
             ## per index traced
             cpit=(sum(coprev*ptentry)/sum(hhc*traceperhhcpt)),
             ## per presumed
             cppe=sum(coprev*ptentry)/sum(hhc*traceperhhcpt*Presumptive),
             ## per screened
             cpcs=sum(coprev*ptentry)/sum(hhc*traceperhhcpt*Screened)
      ),
      by=.(id,country)]
clz <- grep('cp',names(tmp),value=TRUE)
acfcascade <- tmp[,lapply(.SD,function(x)1e2*mean(x)),
                  by=country,.SDcols=clz]

fn <- gh('outdata/ACFcascade') + SA + '.' + ACF + '.csv'
fwrite(acfcascade,file=fn)

## adding in ACF if included
if(ACF>0){
    PT[,cost.soc:=cost.soc + costACF.soc]
    PT[,cost.int:=cost.int + costACF.int]
    ## cases NOTE left same = incidence
    PT[,cases.soc:=cases.soc]
    PT[,cases.int:=cases.int]
    ## ATT
    PT[,ATT.soc:=ATT.soc + ATTACF.soc]
    PT[,ATT.int:=ATT.int + ATTACF.int]
    ## deaths
    PT[,deaths.soc:=deaths.soc + deathsACF.soc]
    PT[,deaths.int:=deaths.int + deathsACF.int]
    PT[,LS:=LS + LSACF]
    ## add in DALYs etc
    PT[,c('dDALY','dDALY0'):=.(dDALY+dDALYacf,dDALY0+dDALY0acf)]
}

## OK to multiply by ptentry - this is split across PT-ages

## aggregate results
PT1 <- PT[,.(cost.soc=sum(cost.soc*ptentry),
             cost.int=sum(cost.int*ptentry),
             cases.soc=sum(cases.soc*ptentry),
             cases.int=sum(cases.int*ptentry),
             ATT.soc=sum(ATT.soc*ptentry),
             ATT.int=sum(ATT.int*ptentry),
             numPT.int=sum(RR*ptentry),
             Dcost=sum(cost.int*ptentry)-sum(cost.soc*ptentry),
             LS=sum(LS*ptentry),
             deaths.soc=sum(deaths.soc*ptentry),
             deaths.int=sum(deaths.int*ptentry),
             dDALY0=sum(dDALY0*ptentry),
             dDALY=sum(dDALY*ptentry)),
          by=.(country,id)]
PT1 <- merge(PT1,CK,by='country') #country iso3 merged on

## same with age - basically just consistent renaming
PT2 <- PT[,.(cost.soc=(cost.soc),
             cost.int=(cost.int),
             cases.soc=(cases.soc),
             cases.int=(cases.int),
             ATT.soc=(ATT.soc),
             ATT.int=(ATT.int),
             numPT.int=(RR),
             Dcost=(cost.int)-(cost.soc),
             LS=(LS),
             deaths.soc=(deaths.soc),
             deaths.int=(deaths.int),
             dDALY0=(dDALY0),
             dDALY=(dDALY)),
          by=.(country,age,id)]
PT2 <- merge(PT2,CK,by='country') #country iso3 merged on

## country name fix
PT1$country <- gsub("ote","ôte",PT1$country)
PT2$country <- gsub("ote","ôte",PT2$country)

## --- CEA and CEAC plots ---

## CEA plot
GP <- ggplot(PT1,aes(dDALY,Dcost)) +
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    geom_point(alpha=0.1,shape=1) +
    geom_abline(data=CETM[threshold %in% c("0.5x GDP","1x GDP")],
                aes(intercept=0,slope=value,col=threshold))+
    facet_wrap(~country) +
    scale_y_continuous(label=comma) +
    xlab('Incremental discounted life-years saved')+
    ylab('Incremental cost (USD)')+
    theme(legend.position = "top" )
if(!shell) GP

## save out
fn1 <- gh('graphs/CEallPT') + SA + '.' + ACF + '.png'
## fn2 <- gh('graphs/CEallPT') + SA + '.' + ACF + '.pdf'
ggsave(GP,file=fn1,w=10,h=10); ## ggsave(GP,file=fn2,w=10,h=10)


## make CEAC data
lz <- seq(from = 0,to=ceactop,length.out = 1000)
pceacd <- list()
for(cn in unique(PT1$country)){
    tmp <- PT1[country==cn,.(Q=dDALY,P=Dcost)]
    pceacd[[cn]] <- data.table(
        country=cn,x=lz,
        y=make.ceac(tmp,lz))
}
pceacd <- rbindlist(pceacd)

## make CEAC plot
PCEAC <- make.ceac.plot(pceacd,xpad=50) +
  xlab('Cost-effectiveness threshold (USD per DALY averted)')
if(!shell) PCEAC

## save out
fn1 <- gh('graphs/CEACpt') + SA + '.' + ACF + '.png'
## fn2 <- gh('graphs/CEACpt') + SA + '.' + ACF + '.pdf'
ggsave(PCEAC,file=fn1,w=10,h=10); ## ggsave(PCEAC,file=fn2,w=10,h=10)

# labelled CEAC
PCEACL <- pceacd %>%
  mutate(label_hjust = rep(c(0.4, 0.15, 0.3, 0.3, 0.2,0.4, 0.3, 0.3, 0.4), each=nrow(ceacd)/length(unique(country)))) %>%
  ggplot() +
  geomtextpath::geom_textline(aes(x,y,
                                  col=country,
                                  label=country,
                                  hjust = label_hjust),
                              size=4, vjust = 0.3, show.legend = FALSE, text_only = TRUE, text_smoothing = 40) + 
  theme_classic() +
  guides(col = guide_legend(order = 1, nrow = 2, title = 'Country')) +
  theme(legend.position = 'top',legend.title = element_blank())+
  # theme(legend.text=element_text(size=rel(0.8)))+
  ggpubr::grids()+
  scale_x_continuous(label=comma,
                     breaks = seq(0,ceactop,250)) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73",
                               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                               "darkorchid1"))+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD per DALY averted)')

PCEACL
ggsave(PCEACL,file=gh('graphs/CEACptL') + SA + '.png',w=7,h=7)

## make PCEACY (0-4 years) data
lz <- seq(from = 0,to=ceactop,length.out = 1000)
pceacdy <- list()
for(cn in unique(PT2$country)){
  tmp <- PT2[age=='0-4',]
  tmp <- tmp[country==cn,.(Q=dDALY,P=Dcost)]
  pceacdy[[cn]] <- data.table(
    country=cn,x=lz,
    y=make.ceac(tmp,lz))
}
pceacdy <- rbindlist(pceacdy)

## make PCEACY plot
PCEACY <- make.ceac.plot(pceacdy,xpad=50)  +
  xlab('Cost-effectiveness threshold (USD per DALY averted)') +
  guides(col = guide_legend(title = 'Country'))
if(!shell) PCEACY

## save out
fn1 <- glue(here('graphs/CEACptY')) + SA + '.' + ACF + '.png'
## fn2 <- glue(here('graphs/CEACpt')) + SA + '.' + ACF + '.pdf'
ggsave(PCEACY,file=fn1,w=10,h=10); ## ggsave(PCEACY,file=fn2,w=10,h=10)

# labelled CEAC
PCEACYL <- pceacdy %>%
  mutate(label_hjust = rep(c(0.35, 0.2, 0.3, 0.3, 0.2,0.2, 0.3, 0.3, 0.25), each=nrow(ceacd)/length(unique(country)))) %>%
  ggplot() +
  geomtextpath::geom_textline(aes(x,y,
                                  col=country,
                                  label=country,
                                  hjust = label_hjust),
                              size=4, vjust = 0.3, show.legend = FALSE, text_only = TRUE, text_smoothing = 40) + 
  theme_classic() +
  guides(col = guide_legend(order = 1, nrow = 2, title = 'Country')) +
  theme(legend.position = 'top',legend.title = element_blank())+
  # theme(legend.text=element_text(size=rel(0.8)))+
  ggpubr::grids()+
  scale_x_continuous(label=comma,
                     breaks = seq(0,ceactop,250)) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73",
                               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                               "darkorchid1"))+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD per DALY averted)')

PCEACYL
ggsave(PCEACYL,file=gh('graphs/CEACYptL') + SA + '.png',w=7,h=7)

## make PCEACY (5-14 years) data
lz <- seq(from = 0,to=ceactop,length.out = 1000)
pceacdo <- list()
for(cn in unique(PT2$country)){
  tmp <- PT2[age=='5-14',]
  tmp <- tmp[country==cn,.(Q=dDALY,P=Dcost)]
  pceacdo[[cn]] <- data.table(
    country=cn,x=lz,
    y=make.ceac(tmp,lz))
}
pceacdo <- rbindlist(pceacdo)

## make PCEACO plot
PCEACO <- make.ceac.plot(pceacdo,xpad=50)  +
  xlab('Cost-effectiveness threshold (USD per DALY averted)') +
  guides(col = guide_legend(title = 'Country'))
if(!shell) PCEACO

## save out
fn1 <- glue(here('graphs/CEACptO')) + SA + '.' + ACF + '.png'
## fn2 <- glue(here('graphs/CEACpt')) + SA + '.' + ACF + '.pdf'
ggsave(PCEACO,file=fn1,w=10,h=10); ## ggsave(PCEACO,file=fn2,w=10,h=10)

# labelled CEAC
PCEACOL <- pceacdo %>%
  mutate(label_hjust = rep(c(0.95, 0.3, 0.95, 0.1, 0.2,0.2, 0.35, 0.3, 0.3), each=nrow(ceacd)/length(unique(country)))) %>%
  ggplot() +
  geomtextpath::geom_textline(aes(x,y,
                                  col=country,
                                  label=country,
                                  hjust = label_hjust),
                              size=4, vjust = 0.3, show.legend = FALSE, text_only = TRUE, text_smoothing = 40) + 
  theme_classic() +
  guides(col = guide_legend(order = 1, nrow = 2, title = 'Country')) +
  theme(legend.position = 'top',legend.title = element_blank())+
  # theme(legend.text=element_text(size=rel(0.8)))+
  ggpubr::grids()+
  scale_x_continuous(label=comma,
                     breaks = seq(0,ceactop,250)) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73",
                               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                               "darkorchid1"))+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD per DALY averted)')

PCEACOL
ggsave(PCEACOL,file=gh('graphs/CEACOptL') + SA + '.png',w=7,h=7)

## output where things X 50%
tmp <- pceacd[y>0.5 & abs(y-0.5)<5e-2] #NOTE may need tuning if sample changes
tmp[,err:=abs(y-0.5)]
tmp[,ermin:=min(err),by=country]
tmp <- tmp[err==ermin]
tmp <- tmp[,.(country,ceac50=round(x))]
tmp
fn <- gh('outdata/CEAC50pt') + SA + '.' + ACF + '.csv'
fwrite(tmp,file=fn)

## ICERs by country
pice <- PT1[is.finite(cost.soc),.(numPT.soc=1, #TODO remove condition
               cost.soc=mean(cost.soc),
               cost.soc.lo=lof(cost.soc),
               cost.soc.hi=hif(cost.soc),
               ## int
               numPT.int=mean(numPT.int),
               numPT.int.lo=lof(numPT.int),
               numPT.int.hi=hif(numPT.int),
               cost.int=mean(cost.int),
               cost.int.lo=lof(cost.int),
               cost.int.hi=hif(cost.int),
               ## dif
               DnumPT=mean(numPT.int-1), #numbers PT
               DnumPT.lo=lof(numPT.int-1),
               DnumPT.hi=hif(numPT.int-1),
               Dcases=mean(cases.int-cases.soc), #TB cases
               Dcases.lo=lof(cases.int-cases.soc), #TB cases
               Dcases.hi=hif(cases.int-cases.soc), #TB cases
               Datt=mean(ATT.int-ATT.soc),   #ATT
               Datt.lo=lof(ATT.int-ATT.soc),   #ATT
               Datt.hi=hif(ATT.int-ATT.soc),   #ATT
               Ddeaths=-mean(LS), #deaths
               Ddeaths.lo=-hif(LS), #deaths
               Ddeaths.hi=-lof(LS), #deaths
               Dcost=mean(Dcost),
               Dcost.lo=lof(Dcost),
               Dcost.hi=hif(Dcost),
               ## dalys
               dDALY0=mean(dDALY0),
               dDALY0.hi=hif(dDALY0),
               dDALY0.lo=lof(dDALY0),
               dDALY=mean(dDALY),
               dDALY.hi=hif(dDALY),
               dDALY.lo=lof(dDALY),
               ## ICER
               ICER=mean(Dcost)/mean(dDALY)
               ),
            by=country]

## NOTE before x100, output SOC unit cost for table 1
Table1PTcost <- PT1[is.finite(cost.soc),
                   .(cost.soc=mean(cost.soc),
                     cost.soc.sd=sd(cost.soc)),by=country]
Table1PTcost[,txt:=paste0(round(cost.soc)," (",round(cost.soc.sd),")")]
Table1PTcost <- transpose(Table1PTcost[,.(country,txt)],
                          make.names = TRUE)
Table1PTcost[,c('condition',
              'variable'):=.('PT',
                             'Cost per PT initiation, US$ (SD)')]

save(Table1PTcost,file=gh('data/Table1PTcost.Rdata'))


## multiply most by 100
vec <- names(pice)
vec <- vec[!vec %in% c('country','ICER')] #all but country and ICER
pice[,(vec):=lapply(.SD, function(x) 100*x), .SDcols = vec]

## table output
picer <- pice[is.finite(cost.soc),.(country=country,
                    numPT.soc=paste0(numPT.soc),
                    cost.soc=bracket(cost.soc,cost.soc.lo,cost.soc.hi),
                    numPT.int=bracket(numPT.int,numPT.int.lo,
                                      numPT.int.hi),
                    cost.int=bracket(cost.int,cost.int.lo,cost.int.hi),
                    DnumPT=bracket(DnumPT,DnumPT.lo,DnumPT.hi),
                    Dcases=bracket(Dcases,Dcases.lo,Dcases.hi),
                    Datt=bracket(Datt,Datt.lo,Datt.hi),
                    Ddeaths=bracket(Ddeaths,Ddeaths.lo,Ddeaths.hi),
                    LY0.dif=bracket(dDALY0,dDALY0.lo,dDALY0.hi),
                    LY.dif=bracket(dDALY,dDALY.lo,dDALY.hi),
                    Dcost=bracket(Dcost,Dcost.lo,Dcost.hi),
                    ICER = format(round(ICER),big.mark = ',')
               )]
picer

fn <- gh('outdata/ICERpt') + SA +'.' + ACF + '.csv'
fwrite(picer,file=fn)

## --- PT stuff for Table 2

Table2PT <- picer[,.(country,
                      treated.soc=numPT.soc,cost.soc,
                      treated.int=numPT.int,cost.int,
                      diff.PT=DnumPT,diff.incTB=Dcases,
                      diff.ATT=Datt,
                      diff.deaths=Ddeaths,
                      diff.LYS=LY0.dif,diff.dLYS=LY.dif,
                      diff.cost=Dcost,
                      ICER)]

fn <- gh('outdata/Table2PT') + SA +'.' + ACF + '.Rdata'
save(Table2PT,file=fn)


piceage <- PT2[is.finite(cost.soc),.(numPT.soc=1,
               cost.soc=mean(cost.soc),
               cost.soc.lo=lof(cost.soc),
               cost.soc.hi=hif(cost.soc),
               ## int
               numPT.int=mean(numPT.int),
               numPT.int.lo=lof(numPT.int),
               numPT.int.hi=hif(numPT.int),
               cost.int=mean(cost.int),
               cost.int.lo=lof(cost.int),
               cost.int.hi=hif(cost.int),
               ## dif
               DnumPT=mean(numPT.int-1), #numbers PT
               DnumPT.lo=lof(numPT.int-1),
               DnumPT.hi=hif(numPT.int-1),
               Dcases=mean(cases.int-cases.soc), #TB cases
               Dcases.lo=lof(cases.int-cases.soc), #TB cases
               Dcases.hi=hif(cases.int-cases.soc), #TB cases
               Datt=mean(ATT.int-ATT.soc),   #ATT
               Datt.lo=lof(ATT.int-ATT.soc),   #ATT
               Datt.hi=hif(ATT.int-ATT.soc),   #ATT
               Ddeaths=-mean(LS), #deaths
               Ddeaths.lo=-hif(LS), #deaths
               Ddeaths.hi=-lof(LS), #deaths
               Dcost=mean(Dcost),
               Dcost.lo=lof(Dcost),
               Dcost.hi=hif(Dcost),
               ## dalys
               dDALY0=mean(dDALY0),
               dDALY0.hi=hif(dDALY0),
               dDALY0.lo=lof(dDALY0),
               dDALY=mean(dDALY),
               dDALY.hi=hif(dDALY),
               dDALY.lo=lof(dDALY),
               ## ICER
               ICER=mean(Dcost)/mean(dDALY)
               ),
            by=.(country,age)]

## multiply most by 100
vec <- names(piceage)
vec <- vec[!vec %in% c('country','age','ICER')] #all but country and ICER
piceage[,(vec):=lapply(.SD, function(x) 100*x), .SDcols = vec]

## table output
picers <- piceage[,.(country=country,
                     age=age,
                     numPT.soc=paste0(numPT.soc),
                     cost.soc=bracket(cost.soc,cost.soc.lo,cost.soc.hi),
                     numPT.int=bracket(numPT.int,numPT.int.lo,
                                       numPT.int.hi),
                     cost.int=bracket(cost.int,cost.int.lo,cost.int.hi),
                     DnumPT=bracket(DnumPT,DnumPT.lo,DnumPT.hi),
                     Dcases=bracket(Dcases,Dcases.lo,Dcases.hi),
                     Datt=bracket(Datt,Datt.lo,Datt.hi),
                     Ddeaths=bracket(Ddeaths,Ddeaths.lo,Ddeaths.hi),
                     LY0.dif=bracket(dDALY0,dDALY0.lo,dDALY0.hi),
                     LY.dif=bracket(dDALY,dDALY.lo,dDALY.hi),
                     Dcost=bracket(Dcost,Dcost.lo,Dcost.hi),
                     ICER = format(round(ICER),big.mark = ',')
               )]
picers

picers[age=='0-4',.(country,numPT.int,ICER)]

fn <- gh('outdata/ICERagept') + SA + '.' + ACF + '.csv'
fwrite(picers,file=fn)

BC$Country <- gsub("ote","ôte",BC$Country)

## ================= BOTH components ============================
## want weighting during baseline
PTT <- merge(T1[,.(country,iso3,id,Dcost.att=Dcost,dDALY.att=dDALY)],
             PT1[,.(country,iso3,id,Dcost.pt=Dcost,dDALY.pt=dDALY)],
             by=c('country','iso3','id'))
PTT <- merge(PTT,BC[,.(country=Country,BR)],by='country') #weight
PTT[,Dcost:=Dcost.att*BR/(1+BR) + Dcost.pt*1/(1+BR)] #weighted costs
PTT[,dDALY:=dDALY.att*BR/(1+BR) + dDALY.pt*1/(1+BR)] #weighted DALYs

# age splits 
PTT2 <- merge(T2[,.(country,iso3,age,id,Dcost.att=Dcost,dDALY.att=dDALY)],
              PT2[,.(country,iso3,age,id,Dcost.pt=Dcost,dDALY.pt=dDALY)],
              by=c('country','iso3','age','id'))
PTT2 <- merge(PTT2,BC[,.(country=Country,BR)],by='country') #weight
PTT2[,Dcost:=Dcost.att*BR/(1+BR) + Dcost.pt*1/(1+BR)] #weighted costs
PTT2[,dDALY:=dDALY.att*BR/(1+BR) + dDALY.pt*1/(1+BR)] #weighted DALYs

## CEA plot
GP <- ggplot(PTT,aes(dDALY,Dcost)) +
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)+
    geom_point(alpha=0.1,shape=1) +
    geom_abline(data=CETM[threshold %in% c("0.5x GDP","1x GDP")],
                aes(intercept=0,slope=value,col=threshold))+
    facet_wrap(~country) +
    scale_y_continuous(label=comma) +
    xlab('Disability-adjusted life-years averted')+
    ylab('Incremental cost (USD)')+
    theme(legend.position = "top" )
if(!shell) GP

## save out
fn1 <- gh('graphs/CEallALL') + SA + '.' + ACF + '.png'
## fn2 <- gh('graphs/CEallALL') + SA + '.' + ACF + '.pdf'
ggsave(GP,file=fn1,w=10,h=10); ## ggsave(GP,file=fn2,w=10,h=10)

## make CEAC data
lz <- seq(from = 0,to=ceactop,length.out = 1000)
bceacd <- list()
for(cn in unique(PT1$country)){
    tmp <- PTT[country==cn,.(Q=dDALY,P=Dcost)]
    bceacd[[cn]] <- data.table(
        country=cn,x=lz,
        y=make.ceac(tmp,lz))
}
bceacd <- rbindlist(bceacd)

## make CEAC plot
BCEAC <- make.ceac.plot(bceacd,xpad=50) +
  xlab('Cost-effectiveness threshold (USD per DALY averted)')
if(!shell) BCEAC

## save out
fn1 <- gh('graphs/CEACall') + SA + '.' + ACF + '.png'
## fn2 <- gh('graphs/CEACall') + SA + '.' + ACF + '.pdf'
ggsave(BCEAC,file=fn1,w=10,h=10); ## ggsave(PCEAC,file=fn2,w=10,h=10)

# labelled BCEAC
BCEACL <- bceacd %>%
  mutate(label_hjust = rep(c(0.85, 0.4, 0.25, 0.4, 0.95,0.25, 0.25, 0.15, 0.85), each=nrow(ceacd)/length(unique(country)))) %>%
  ggplot() +
  geomtextpath::geom_textline(aes(x,y,
                                  col=country,
                                  label=country,
                                  hjust = label_hjust),
                              size=4, vjust = 0.3, show.legend = FALSE, text_only = TRUE, text_smoothing = 40) + 
  theme_classic() +
  guides(col = guide_legend(order = 1, nrow = 2, title = 'Country')) +
  theme(legend.position = 'top',legend.title = element_blank())+
  # theme(legend.text=element_text(size=rel(0.8)))+
  ggpubr::grids()+
  scale_x_continuous(label=comma,
                     breaks = seq(0,ceactop,250)) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73",
                               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                               "darkorchid1"))+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD per DALY averted)')

BCEACL
ggsave(BCEACL,file=gh('graphs/CEACallL') + SA + '.png',w=7,h=7)

## make BCEACY (0-4 years) data
lz <- seq(from = 0,to=ceactop,length.out = 1000)
bceacdy <- list()
for(cn in unique(PTT2$country)){
  tmp <- PTT2[age=='0-4',]
  tmp <- tmp[country==cn,.(Q=dDALY,P=Dcost)]
  bceacdy[[cn]] <- data.table(
    country=cn,x=lz,
    y=make.ceac(tmp,lz))
}
bceacdy <- rbindlist(bceacdy)

## make CEAC plot
BCEACY <- make.ceac.plot(bceacdy,xpad=50)  +
  xlab('Cost-effectiveness threshold (USD per DALY averted)') +
  guides(col = guide_legend(title = 'Country'))
if(!shell) BCEACY

## save out
fn1 <- glue(here('graphs/CEACallY')) + SA + '.' + ACF + '.png'
## fn2 <- glue(here('graphs/CEACallY')) + SA + '.' + ACF + '.pdf'
ggsave(BCEACY,file=fn1,w=10,h=10); ## ggsave(BCEACY,file=fn2,w=10,h=10)

# labelled BCEAC
BCEACYL <- bceacdy %>%
  mutate(label_hjust = rep(c(0.95, 0.7, 0.2, 0.3, 0.9,0.4, 0.25, 0.1, 0.85), each=nrow(ceacd)/length(unique(country)))) %>%
  ggplot() +
  geomtextpath::geom_textline(aes(x,y,
                                  col=country,
                                  label=country,
                                  hjust = label_hjust),
                              size=4, vjust = 0.3, show.legend = FALSE, text_only = TRUE, text_smoothing = 40) + 
  theme_classic() +
  guides(col = guide_legend(order = 1, nrow = 2, title = 'Country')) +
  theme(legend.position = 'top',legend.title = element_blank())+
  # theme(legend.text=element_text(size=rel(0.8)))+
  ggpubr::grids()+
  scale_x_continuous(label=comma,
                     breaks = seq(0,ceactop,250)) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73",
                               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                               "darkorchid1"))+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD per DALY averted)')

BCEACYL
ggsave(BCEACYL,file=gh('graphs/CEACYallL') + SA + '.png',w=7,h=7)

## make BCEACO (5-14 years) data
lz <- seq(from = 0,to=ceactop,length.out = 1000)
bceacdo <- list()
for(cn in unique(PTT2$country)){
  tmp <- PTT2[age=='5-14',]
  tmp <- tmp[country==cn,.(Q=dDALY,P=Dcost)]
  bceacdo[[cn]] <- data.table(
    country=cn,x=lz,
    y=make.ceac(tmp,lz))
}
bceacdo <- rbindlist(bceacdo)

## make CEAC plot
BCEACO <- make.ceac.plot(bceacdo,xpad=50)  +
  xlab('Cost-effectiveness threshold (USD per DALY averted)') +
  guides(col = guide_legend(title = 'Country'))
if(!shell) BCEACO

## save out
fn1 <- glue(here('graphs/CEACallO')) + SA + '.' + ACF + '.png'
## fn2 <- glue(here('graphs/CEACallO')) + SA + '.' + ACF + '.pdf'
ggsave(BCEACO,file=fn1,w=10,h=10); ## ggsave(BCEACO,file=fn2,w=10,h=10)

# labelled BCEACO
BCEACOL <- bceacdo %>%
  mutate(label_hjust = rep(c(0.8, 0.3, 0.25, 0.5, 0.98,0.45, 0.25, 0.1, 0.85), each=nrow(ceacd)/length(unique(country)))) %>%
  ggplot() +
  geomtextpath::geom_textline(aes(x,y,
                                  col=country,
                                  label=country,
                                  hjust = label_hjust),
                              size=4, vjust = 0.3, show.legend = FALSE, text_only = TRUE, text_smoothing = 40) + 
  theme_classic() +
  guides(col = guide_legend(order = 1, nrow = 2, title = 'Country')) +
  theme(legend.position = 'top',legend.title = element_blank())+
  # theme(legend.text=element_text(size=rel(0.8)))+
  ggpubr::grids()+
  scale_x_continuous(label=comma,
                     breaks = seq(0,ceactop,250)) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73",
                               "#F0E442", "#0072B2","#D55E00", "#CC79A7",
                               "darkorchid1"))+
  ylab('Probability cost-effective')+
  xlab('Cost-effectiveness threshold (USD per DALY averted)')

BCEACOL
ggsave(BCEACOL,file=gh('graphs/CEACOallL') + SA + '.png',w=7,h=7)
##  combined ICERS
iceb <- PTT[,.(ICER=mean(Dcost)/mean(dDALY)),
          by=.(country,iso3)]

icebrr <- iceb[,.(country,iso3,
                  ICER = format(round(ICER),big.mark = ',') )]
icebrr

fn <- gh('outdata/ICERall') + SA + '.' + ACF + '.csv'
fwrite(icebrr,file=fn)

## output where things X 50%
tmp <- bceacd[y>0.5 & abs(y-0.5)<5e-2] #NOTE may need tuning if sample changes
tmp[,err:=abs(y-0.5)]
tmp[,ermin:=min(err),by=country]
tmp <- tmp[err==ermin]
tmp <- tmp[,.(country,ceac50=round(x))]
tmp
fn <- gh('outdata/CEAC50all') + SA + '.' + ACF + '.csv'
fwrite(tmp,file=fn)
PTT

## --- make a larger set of outputs from both components, and weight appropriately
A <- T1[,.(country,iso3,id,
           cost.soc.ATT=cost.soc,
           cost.int.ATT=cost.int,
           tx=1,tx.int=tx,
           Dcost.ATT=Dcost,
           Ddeaths.ATT=-LS,
           deaths.soc.ATT=deaths.soc,
           deaths.int.ATT=deaths.int,
           dDALY0.ATT=dDALY0,
           dDALY.ATT=dDALY)]
B <- PT1[,.(country,iso3,id,
            cases.soc,cases.int,
            ATT.soc,ATT.int,
            numPT.soc=1,numPT.int,
            cost.soc.TPT=cost.soc,
            cost.int.TPT=cost.int,
            Ddeaths.TPT=-LS,
            deaths.soc.TPT=deaths.soc,
            deaths.int.TPT=deaths.int,
            dDALY0.TPT=dDALY0,
            Dcost.TPT=Dcost,
            dDALY.TPT=dDALY)]

## fix country name
BC$Country <- gsub("ote","ôte",BC$Country)

## merge in
PTA <- merge(A,B,by=c('country','iso3','id'))
PTA <- merge(PTA,BC[,.(country=Country,BR)],by='country') #weight

## multiply by relevant factors
vec.att <- setdiff(names(A),c('country','iso3','id'))
vec.tpt <- setdiff(names(B),c('country','iso3','id'))
PTA[,(vec.att):=lapply(.SD, function(x) 100*x* BR/(1+BR)), .SDcols = vec.att]
PTA[,(vec.tpt):=lapply(.SD, function(x) 100*x* 1 /(1+BR)), .SDcols = vec.tpt]

## -- combined variables for table 2 3rd part
## SOC
PTA[,r.start.soc:=tx+numPT.soc] #Started treatment
PTA[,r.cost.soc:=cost.soc.ATT+cost.soc.TPT] #Cost
## INT
PTA[,r.start.int:=tx.int+numPT.int] #Started treatment
PTA[,r.cost.int:=cost.int.ATT+cost.int.TPT] #Cost
## differences
PTA[,r.started.TPT:=numPT.int-numPT.soc]
PTA[,r.incTB:=cases.int-cases.soc]
PTA[,r.started.ATT:=(tx.int+ATT.int) - (tx+ATT.soc)]
PTA[,r.TBdeaths:=Ddeaths.ATT+Ddeaths.TPT]
PTA[,r.LYS:=dDALY0.ATT+dDALY0.TPT]
PTA[,r.dLYS:=dDALY.ATT+dDALY.TPT]
PTA[,r.Dcost:=Dcost.ATT+Dcost.TPT]


## ----- table 2 output for combined intervention
bice <- PTA[,.(
  r.start.soc=mean(r.start.soc),r.start.soc.lo=lof(r.start.soc),r.start.soc.hi=hif(r.start.soc),
  r.cost.soc=mean(r.cost.soc),r.cost.soc.lo=lof(r.cost.soc),r.cost.soc.hi=hif(r.cost.soc),
  r.start.int=mean(r.start.int),r.start.int.lo=lof(r.start.int),r.start.int.hi=hif(r.start.int),
  r.cost.int=mean(r.cost.int),r.cost.int.lo=lof(r.cost.int),r.cost.int.hi=hif(r.cost.int),
  r.started.TPT=mean(r.started.TPT),r.started.TPT.lo=lof(r.started.TPT),r.started.TPT.hi=hif(r.started.TPT),
  r.incTB=mean(r.incTB),r.incTB.lo=lof(r.incTB),r.incTB.hi=hif(r.incTB),
  r.started.ATT=mean(r.started.ATT),r.started.ATT.lo=lof(r.started.ATT),r.started.ATT.hi=hif(r.started.ATT),
  r.TBdeaths=mean(r.TBdeaths),r.TBdeaths.lo=lof(r.TBdeaths),r.TBdeaths.hi=hif(r.TBdeaths),
  r.LYS=mean(r.LYS),r.LYS.lo=lof(r.LYS),r.LYS.hi=hif(r.LYS),
  r.dLYS=mean(r.dLYS),r.dLYS.lo=lof(r.dLYS),r.dLYS.hi=hif(r.dLYS),
  r.Dcost=mean(r.Dcost),r.Dcost.lo=lof(r.Dcost),r.Dcost.hi=hif(r.Dcost),
  ## ICER
  ICER=mean(r.Dcost)/mean(r.dLYS)
),
by=country]

## table output
bicer <- bice[,.(country=country,
                 r.start.soc=bracket(r.start.soc,r.start.soc.lo,r.start.soc.hi),
                 r.cost.soc=bracket(r.cost.soc,r.cost.soc.lo,r.cost.soc.hi),
                 r.start.int=bracket(r.start.int,r.start.int.lo,r.start.int.hi),
                 r.cost.int=bracket(r.cost.int,r.cost.int.lo,r.cost.int.hi),
                 r.started.TPT=bracket(r.started.TPT,r.started.TPT.lo,r.started.TPT.hi),
                 r.incTB=bracket(r.incTB,r.incTB.lo,r.incTB.hi),
                 r.started.ATT=bracket(r.started.ATT,r.started.ATT.lo,r.started.ATT.hi),
                 r.TBdeaths=bracket(r.TBdeaths,r.TBdeaths.lo,r.TBdeaths.hi),
                 r.LYS=bracket(r.LYS,r.LYS.lo,r.LYS.hi),
                 r.dLYS=bracket(r.dLYS,r.dLYS.lo,r.dLYS.hi),
                 r.Dcost=bracket(r.Dcost,r.Dcost.lo,r.Dcost.hi),
                 ICER = format(round(ICER),big.mark = ',')
               )]
bicer

bicer[,r.start.soc:='100'] #cosmetic/clarity correction

Table2both <- bicer

fn <- gh('outdata/Table2both') + SA +'.' + ACF + '.Rdata'
save(Table2both,file=fn)

## ---- combined CEAC plot ---
## CEAC #from ceacd
## PCEAC #from pceacd
## BCEAC #from bceacd
## TODO country names not iso3

ttls <- c('Intensified case-finding intervention',
          'Household contact management intervention',
          'Combined CaP TB intervention package')
CAP <- list(CEAC,PCEAC,BCEAC)
for(i in 1:3) CAP[[i]] <- CAP[[i]] + ggtitle(ttls[i])

CAPP <- ggarrange(plotlist = CAP,
                  labels=paste0(letters[1:3],")"),
                  common.legend = TRUE,legend="bottom",ncol=1)

fn1 <- glue(here('graphs/allCEACs')) + SA + '.' + ACF + '.png'
## fn2 <- glue(here('graphs/allCEACs')) + SA + '.' + ACF + '.pdf'
ggsave(CAPP,file=fn1,w=8,h=15); #ggsave(CAPP,file=fn2,w=8,h=15);

# labelled CEACs
CAPL <- list(CEACL,PCEACL,BCEACL)
for(i in 1:3) CAPL[[i]] <- CAPL[[i]] + ggtitle(ttls[i])

CAPPL <- ggarrange(plotlist = CAPL,
                  labels=paste0(letters[1:3],")"),
                  common.legend = TRUE,legend="bottom",ncol=1)

fn1 <- glue(here('graphs/allCEACLs')) + SA + '.' + ACF + '.png'
fn2 <- glue(here('graphs/allCEACLs')) + SA + '.' + ACF + '.pdf'
ggsave(CAPPL,file=fn1,w=8,h=15); 
ggsave(CAPPL,file=fn2,w=8,h=15);

# 0-4 years
CAP <- list(CEACYL,PCEACYL,BCEACYL)
for(i in 1:3) CAP[[i]] <- CAP[[i]] + ggtitle(ttls[i])

CAPPY <- ggarrange(plotlist = CAP,
                   labels=paste0(letters[1:3],")"),
                   common.legend = TRUE,legend="bottom",ncol=1)

fn1 <- glue(here('graphs/allCEACLsY')) + SA + '.' + ACF + '.png'
fn2 <- glue(here('graphs/allCEACLsY')) + SA + '.' + ACF + '.pdf'
ggsave(CAPPY,file=fn1,w=8,h=15); ggsave(CAPPY,file=fn2,w=8,h=15);

# 5-14 years
CAP <- list(CEACOL,PCEACOL,BCEACOL)
for(i in 1:3) CAP[[i]] <- CAP[[i]] + ggtitle(ttls[i])

CAPPO <- ggarrange(plotlist = CAP,
                   labels=paste0(letters[1:3],")"),
                   common.legend = TRUE,legend="bottom",ncol=1)

fn1 <- glue(here('graphs/allCEACLsO')) + SA + '.' + ACF + '.png'
fn2 <- glue(here('graphs/allCEACLsO')) + SA + '.' + ACF + '.pdf'
ggsave(CAPPO,file=fn1,w=8,h=15); ggsave(CAPPO,file=fn2,w=8,h=15);

# Age groups 
ttls <- c('0-4 years',
          '5-14 years',
          '0-14 years')
CAP <- list(PCEACYL,PCEACOL,PCEACL)
for(i in 1:3) CAP[[i]] <- CAP[[i]] + ggtitle(ttls[i])

CAPPA <- ggarrange(plotlist = CAP,
                   labels=paste0(letters[1:3],")"),
                   common.legend = TRUE,legend="bottom",ncol=1)

fn1 <- glue(here('graphs/allCEACsA')) + SA + '.' + ACF + '.png'
## fn2 <- glue(here('graphs/allCEACsO')) + SA + '.' + ACF + '.pdf'
ggsave(CAPPA,file=fn1,w=8,h=15); #ggsave(CAPPO,file=fn2,w=8,h=15);
## =============================================

## NOTE
## SHARED tippi folder:
## https://drive.google.com/drive/folders/18i-KlJtvHVUjcFlHSFPsU6NIKbH2wVgq
