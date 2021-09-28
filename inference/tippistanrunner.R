## capture commandline args
args <- commandArgs()
qty <- args[3]
page <- as.integer(args[4]) #age
print(qty)
print(page)

## stop()
## qty <- 'tx'
## page <- 2 #NOTE change this to change age group

## NOTE be careful with site names confidentiality!
## library
library(here)
library(rstan)
library(bayesplot)

options(mc.cores = 2) #NOTE
rstan_options(auto_write = TRUE)

## ========== functions ===========
source(here('../dataprep/tippifunctions.R'))


## compile model
hsm <- stan_model(file=here('stan/tippihm.stan')) #hierarchical model

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

## ===============
## country key
CK <- D[,.(fst=min(index),lst=max(index)),by=.(country=Country)]
CK[,cno:=1:nrow(CK)]
fn <- glue(here('data/')) + qty + 'CK.' + shhs[page,age] + '.Rdata'
save(CK,file=fn)


## data for stan
D4S <- list(
    ## number of countries
    NC = length(unique(D[,Country])),
    ## total number of sites
    NST = length(unique(D[,Facility])),
    ## site time in intervention
    STi = D[,Intervention.FT],
    ## site time in baseline
    STb = D[,Baseline.FT],
    ## number intervention
    DXi = D[,Intervention.Num],
    ## number baseline
    DXb = D[,Baseline.Num],
    ## start index for country
    st = D[,.(fst=min(index)),by=Country][,fst],
    ## end index for country
    nd = D[,.(lst=max(index)),by=Country][,lst],
    lcmsig = 10, #effect mean prior variance
    lcssig = 10, #effect variance prior variance
    ratessig = 10, #rate prior variance
    sigsig = 20    #variance prior variance
)


fn <- glue(here('data/')) + qty + 'stanin.' + shhs[page,age] + '.Rdata'
save(D4S,file=fn)

## sample from model
niter <- 1e4; nchains <- 2 #TODO change
samp0 <- sampling(hsm,
                  data = D4S,
                  ## control = list(max_treedepth = 15),
                  iter = niter,
                  chains = nchains,
                  verbose = FALSE)

## ## TODO needs some tuning for page==2

## plotting
to_plot <- c('lcrr',  'lcm')
TP <- traceplot(samp0, pars = to_plot) + rotx
fn <- glue(here('graphs/')) + qty + 'TP' + shhs[page,age] + '.png'
ggsave(TP,file=fn,w=7,h=5)

## X <- as.array(samp0)

## data from MCMC
MC <- summary(samp0,to_plot)$summary[,c('2.5%','50%','97.5%')]
MC <- exp(MC)
MC <- as.data.table(MC)
MC <- cbind(country=c(CK$country,'SUMMARY'),MC)
MC$country <- factor(MC$country,levels=MC$country)

## write out
fn <- glue(here('outdata/')) + qty + 'MC' + shhs[page,age] + '.csv'
fwrite(MC,file=fn)

to_get <- c('lcrr','sig','lcm','lcs')
SS <- extract(samp0,to_get)
SS <- as.data.table(SS$lcrr)
names(SS) <- CK$country
fn <- glue(here('outdata/')) + qty + 'SS' + shhs[page,age] + '.Rdata'
save(SS,file=fn)

## save out some diagnostic statistics
dxstats <- summary(samp0, pars = to_get)$summary[,c('n_eff','Rhat')]
nn <- rownames(dxstats)
nn[1:D4S$NC] <- as.character(MC$country[1:D4S$NC])
dxstats <- as.data.table(dxstats)
dxstats[,variable:=nn]
setcolorder(dxstats,c('variable','n_eff','Rhat'))
dxstats$Rhat <- paste0(round(dxstats$Rhat,4))
dxstats$n_eff <- paste0(round(dxstats$n_eff,0))

## write out
fn <- glue(here('outdata/')) + qty + 'MCdx' + shhs[page,age] + '.csv'
fwrite(dxstats,file=fn)

## empirical data
## sites
siteeffects <- D[,.(site.effect = (Intervention.Num/Intervention.FT)/
                        (Baseline.Num/Baseline.FT),
                    int.number = Intervention.Num,
                    country=Country,
                    index)]
## country
countryeffects <- D[,.(country.effect =
                          (sum(Intervention.Num)/sum(Intervention.FT))/
                           (sum(Baseline.Num)/sum(Baseline.FT)),
                       int.number=sum(Intervention.Num) ),
                    by=.(country=Country)]


MC$country <- factor(MC$country,levels=rev(MC$country),ordered = TRUE)

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
    scale_y_sqrt(limits=c(0,150)) + #NOTE one DRC site omitted
    ylab('Rate ratio (square root scale)')+
    xlab('Country')+
    coord_flip() +
    theme_classic() + ggpubr::grids()+ ggtitle(shhs[page,aged])+
    labs(size='Number (intervention)')+
    guides(colour="none") + 
    theme(legend.position = 'bottom')


fn <- glue(here('graphs/')) + qty + 'MA3' + shhs[page,age] + '.pdf'
fn2 <- glue(here('graphs/')) + qty + 'MA3' + shhs[page,age] + '.png'
ggsave(MAP,file=fn,w=7,h=7)
ggsave(MAP,file=fn2,w=7,h=7)
fn <- glue(here('data/')) + qty + 'MA3' + shhs[page,age] + '.Rdata'
save(MAP,file=fn)
