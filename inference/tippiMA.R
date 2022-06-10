args <- commandArgs()
qty <- args[3]
page <- as.integer(args[4]) #age
print(qty)
print(page)

## ## ## stop()
## qty <- 'px'
## page <- 2 #NOTE change this to change age group


## frequentist MA of results
library(here)
library(rstanarm)
library(lme4)
library(metafor)

## ========== functions ===========
source(here('../dataprep/tippifunctions.R'))


## load data
if(qty=='tx'){ #select qty
    D <- fread(here('../dataprep/outdata/T.csv'))
} else{
    D <- fread(here('../dataprep/outdata/P.csv'))
}
D <- D[age==shhs[page,aged]] #select age
D[,Facility:=site]
setkey(D,Country,Facility)

## additional variables
D[,IRR:=Intervention.Rate/Baseline.Rate]
D[,facility:=as.integer(factor(Facility))]
D[,country:=as.integer(factor(Country))]
str(D)

## summary
SMY <- D[,.(IRR=(sum(Intervention.Num)/sum(Intervention.FT))/
              (sum(Baseline.Num)/sum(Baseline.FT))),
         by=Country]

DM <- melt(D[,.(country,site,
                Baseline.Num,Intervention.Num,
                Baseline.FT,Intervention.FT)],
           id=c('country','site'))
DM[,c('period','qty'):=tstrsplit(variable,split="\\.")]
DM <- dcast(DM,country+site+period ~ qty)
DM$Num <- as.integer(DM$Num)
DM$site <- as.integer(DM$site)
str(DM)


mm <- glmer(formula = Num ~ (period|country/site),
            offset=log(FT),
            data=DM,
            family = poisson(link = "log"))


print(mm)
SMY
exp(coef(mm)$country)

## freq smy
fsmy <- data.table(country=D[,unique(Country)],
                   RR.mid=exp(ranef(mm)$country)[,2])

save(fsmy,file=gh('outdata/fsmy_{qty}_{shhs[page,age]}.Rdata'))



smm <- stan_glmer(formula = Num ~ (period|country/site),
                  offset=log(FT),
                  data=DM,
                  chains=2,cores=1,
                  family = "poisson")

## MT <- as.matrix(smm,regex_pars = c('periodIntervention'))
## str(MT)
## colnames(MT)


nmz <- c(
  "b[periodIntervention country:1]",
  "b[periodIntervention country:2]",
  "b[periodIntervention country:3]",
  "b[periodIntervention country:4]",
  "b[periodIntervention country:5]",
  "b[periodIntervention country:6]",
  "b[periodIntervention country:7]",
  "b[periodIntervention country:8]",
  "b[periodIntervention country:9]"
)

MT <- as.matrix(smm,pars = nmz)
head(MT)
MT <- exp(MT)
head(MT)
save(MT,file=gh('outdata/MT_{qty}_{shhs[page,age]}.Rdata'))


## bayes smy
bsmy <- data.table(
  country=D[,unique(Country)],
  RR.lo=apply(MT,2,lof),
  RR.hi=apply(MT,2,hif),
  RR.mid=apply(MT,2,mean)
)

save(bsmy,file=gh('outdata/bsmy_{qty}_{shhs[page,age]}.Rdata'))

## RMA for each country separately
RMAR <- list()
for(cn in D[,unique(Country)]){
  print(cn)
  tmp <- D[Country==cn]
  modcn <- rma.glmm(measure = "IRR",
                    model='CM.EL',
                    data = tmp,
                    x2i = Baseline.Num, t2i = Baseline.FT,
                    x1i = Intervention.Num, t1i = Intervention.FT)
  bz <- c(modcn$b,modcn$ci.lb,modcn$ci.ub)
  bz <- exp(bz)
  RMAR[[cn]] <- data.table(country=cn,RR=bz[1],RR.lo=bz[2],RR.hi=bz[3],lRR.se=modcn$se)
}
RMAR <- rbindlist(RMAR)

save(RMAR,file=gh('outdata/RMAR_{qty}_{shhs[page,age]}.Rdata'))


## empirical data
## sites
siteeffects <- D[,.(site.effect = (Intervention.Num/Intervention.FT)/
                        (Baseline.Num/Baseline.FT),
                    int.number = Intervention.Num,
                    country=Country)]
## country
countryeffects <- D[,.(country.effect =
                          (sum(Intervention.Num)/sum(Intervention.FT))/
                           (sum(Baseline.Num)/sum(Baseline.FT)),
                       int.number=sum(Intervention.Num) ),
                    by=.(country=Country)]
## MC$country <- factor(MC$country,levels=rev(MC$country),ordered = TRUE)

## harmonization
siteeffects[,c('RR.mid','RR.lo','RR.hi'):=NA]
countryeffects[,c('RR.mid','RR.lo','RR.hi'):=NA]
fsmy[,c('RR.lo','RR.hi'):=NA]

## compare methods
ceall <- rbindlist(list(countryeffects[,.(country,RR.mid=country.effect,RR.lo,RR.hi,type='empirical')],
                        RMAR[,.(country,RR.mid=RR,RR.lo,RR.hi,type='countrywise RMA')],
                        fsmy[,.(country,RR.mid,RR.lo,RR.hi,type='frequentist MLM')],
                        bsmy[,.(country,RR.mid,RR.lo,RR.hi,type='Bayesian MLM')]))


ttl <- glue('{qty}: age {shhs[page,aged]}')
psn <- position_dodge(0.1)
GP <- ggplot(ceall,aes(country,y=RR.mid,ymin=RR.lo,ymax=RR.hi,col=type,shape=type))+
  geom_point(position=psn)+geom_errorbar(width=0,position=psn)+
  xlab('Country')+ylab('Rate ratio')+
  coord_flip() + theme_classic()+ggpubr::grids()+
  theme(legend.position = 'top')+
  ggtitle(ttl)+scale_y_sqrt(limits = c(0,min(100,max(ceall$RR.hi,na.rm = TRUE))))
## GP

ggsave(GP,file=gh('graphs/compare_{qty}_{shhs[page,age]}.pdf'),w=8,h=8)


## reorder
lvl <- unique(as.character(D$Country))
lvl <- rev(lvl)
siteeffects$country <- factor(siteeffects$country,levels=lvl)
countryeffects$country <- factor(countryeffects$country,levels=lvl)
fsmy$country <- factor(fsmy$country,levels=lvl)
bsmy$country <- factor(bsmy$country,levels=lvl)
clz <- RColorBrewer::brewer.pal(length(lvl),'Paired')


MAP <- 
  ggplot(data=bsmy,aes(x=country,
                       y=RR.mid,ymin=RR.lo,ymax=RR.hi,
                       col=country)) +
  geom_point(size=2) +
  geom_point(data=fsmy,size=2,col=2) +
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
  theme_classic() + ggpubr::grids()+ ggtitle(shhs[page,aged])+
  labs(size='Number (intervention)')+
  guides(colour="none") + 
  theme(legend.position = 'bottom')
MAP


fn <- glue(here('graphs/')) + qty + 'MA3' + shhs[page,age] + '.pdf'
fn2 <- glue(here('graphs/')) + qty + 'MA3' + shhs[page,age] + '.png'
ggsave(MAP,file=fn,w=7,h=7)
ggsave(MAP,file=fn2,w=7,h=7)
fn <- glue(here('data/')) + qty + 'MA3' + shhs[page,age] + '.Rdata'
save(MAP,file=fn)




## https://mc-stan.org/rstanarm/reference/as.matrix.stanreg.html

## library(lattice)
## xyplot(data=DM,Num/FT ~ period | country,
##        groups = site,
##        pch=19, lwd=2, type=c('p'))

## https://stats.stackexchange.com/questions/122009/extracting-slopes-for-cases-from-a-mixed-effects-model-lme4

## test <- model.matrix(object = Num ~ FT * (site + (period|country/site)),
##                      data=DM)

## test <- model.matrix(object = Num ~ FT * (site + period),
##                      data=DM)
## head(test)


## dm <- melt(D[,.(Country,Facility,
##                 Baseline.Num,Intervention.Num,
##                 Baseline.FT,Intervention.FT)],
##            id=c('Country','Facility'))
## dm[,c('period','qty'):=tstrsplit(variable,split="\\.")]
## dm <- dcast(dm,Country+Facility+period ~ qty)
## dm$Num <- as.integer(dm$Num)
## dm[,Country:=factor(Country)]
## dm[,Facility:=factor(Facility)]


## test <- model.matrix(object = Num ~ FT * (Facility + period),
##                      data=dm)
## test[1:6,1:10]






## str(dat.konstantopoulos2011)

## ## warm-up modelling everything for single country
## cn <- 'Uganda'
## tmp <- D[Country==cn]
## modcn <- rma.glmm(measure = "IRR",model='CM.EL',
##                    data = tmp,
##                    x2i = Baseline.Num, t2i = Baseline.FT,
##                    x1i = Intervention.Num, t1i = Intervention.FT)

## print(modcn)
## predict(modcn, transf=exp, digits=2)
## ## forest(modcn,transf = exp)
## ti <- glue(cn) + ', ' + shhs[page,aged]  + ', ' + qty
## forest(modcn,main=ti)#,transf = exp)
## tmp


## ## loop over all countries
## for(cn in D[,unique(Country)]){
##   print(cn)
##   tmp <- D[Country==cn]
##   print(tmp[,(sum(Intervention.Num)/sum(Intervention.FT))/
##              (sum(Baseline.Num)/sum(Baseline.FT))])
##   modcn <- rma.glmm(measure = "IRR",model='CM.EL',
##                      data = tmp,
##                      x2i = Baseline.Num, t2i = Baseline.FT,
##                      x1i = Intervention.Num, t1i = Intervention.FT)
##   ti <- glue(cn) + ', ' + shhs[page,aged]  + ', ' + qty
##   fn <- gsub("[ |,|-]","",ti)
##   fn <- gh('graphs/f_{fn}.pdf')
##   pdf(file=fn)
##   forest(modcn,main=ti,transf = exp)
##   dev.off()
## }

## ## 3 level model
## modall <- rma.glmm(measure = "IRR",model='CM.EL',
##                    mods = ~ 1 | country/facility,
##                    data = D,
##                    x2i = Baseline.Num, t2i = Baseline.FT,
##                    x1i = Intervention.Num, t1i = Intervention.FT)


## modall <- rma.mv(measure = "IRR",model='CM.EL',
##                    mods = ~ 1 | country/facility,
##                    data = D,
##                    x2i = Baseline.Num, t2i = Baseline.FT,
##                    x1i = Intervention.Num, t1i = Intervention.FT)


## ## multi-levle poisson-normal?
## fitstats.rma(modall)

## print(modall)
## summary(modall)
## predict(modall, transf=exp, digits=2)
## ## forest(modcn,transf = exp)
## ti <- glue(cn) + ', ' + shhs[page,aged]  + ', ' + qty
## forest(modall)#,transf = exp)



## str(modall)


## library(merTools)
## predictInterval(mm)

##

## library(parameters)
## model_parameters(mm)
## model_parameters(mm,effects = "all")
