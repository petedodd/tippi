args <- commandArgs()
qty <- args[3]
page <- as.integer(args[4]) #age
print(qty)
print(page)

## ## ## stop()
## qty <- 'tx'
## page <- 2 #NOTE change this to change age group


## frequentist MA of results
library(here)
library(rstanarm)

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
DM[,fcountry:=factor(country)]


## fit model
smm <- stan_glmer(formula = Num ~ period * fcountry + (1|site),
                  offset=log(FT),
                  data=DM,
                  chains=2,cores=1,iter=6e3,warmup=1e3,
                  family = "poisson")

## MT <- as.matrix(smm,regex_pars = c('periodIntervention'))
## str(MT)
## colnames(MT)
## prior_summary(smm)
smt <- summary(smm)
nnz <- row.names(smt)
smt <- as.data.table(smt)
smt[,variable:=nnz]
smt <- smt[!variable %in% c('mean_PPD',
                            'log-posterior')]
save(smt,file=gh('outdata/smt_{qty}_{shhs[page,age]}.Rdata'))

## convergence diagnostics
cdx <- c(smt[,median(Rhat)],smt[,range(Rhat)],
         smt[,median(n_eff)],smt[,range(n_eff)])
cat(cdx,file=gh('outdata/cdx_{qty}_{shhs[page,age]}.txt'))

## parameters to extract
nmz <- c(
  "periodIntervention",
  "periodIntervention:fcountry2",
  "periodIntervention:fcountry3",
  "periodIntervention:fcountry4",
  "periodIntervention:fcountry5",
  "periodIntervention:fcountry6",
  "periodIntervention:fcountry7",
  "periodIntervention:fcountry8",
  "periodIntervention:fcountry9"
)


MT <- as.matrix(smm,pars = nmz)
MT[,2:ncol(MT)] <- MT[,2:ncol(MT)] + MT[,1] #add on ict
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
## harmonization
siteeffects[,c('RR.mid','RR.lo','RR.hi'):=NA]
countryeffects[,c('RR.mid','RR.lo','RR.hi'):=NA]

## compare methods
ceall <- rbindlist(list(countryeffects[,.(country,RR.mid=country.effect,RR.lo,RR.hi,type='empirical')],
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
bsmy$country <- factor(bsmy$country,levels=lvl)
clz <- RColorBrewer::brewer.pal(length(lvl),'Paired')


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



