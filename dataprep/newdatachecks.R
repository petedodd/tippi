library(here)
library(glue)

LP <- glue(scan(here('graphs/petelocalpath.txt'),what="string"))

## clinics to drop for reasons of non comparability
dtdrop <- c("Beatrice Road Infectious Disease Hospital (BRIDH)",
            "Dulibadzimu Clinic","Hopley Clinic","Kuwadzana Polyclinic",
            "Mabvuku Polyclinic","Mbare Polyclinic","Shashe Clinic")

ptdrop <- c("Ntungamo Ngoma HC III",
            "Shashe Clinic")

## ========== functions ===========
source(here('tippifunctions.R'))
## NOTE: raw data kept out of project folder so no site names make it in

## ================= dx init ==================

## CMR
fn <- LP+"cameroon/TB Diagnosed up to 31Dec2020_Cam.xlsx"
## dx014
D014cmr <- rexel(fn,sheet=1,skip = 1)
## dx04
D04cmr <- rexel(fn,sheet=2,skip = 1)
## dx514
D514cmr <- rexel(fn,sheet=3,skip = 1)

## KEN
fn <- LP+"kenya/TB Diagnosed up to 31Dec2020_Ken.xlsx"
## dx014
D014ken <- rexel(fn,sheet=1,skip = 1)
## dx04
D04ken <- rexel(fn,sheet=2,skip = 1)
## dx514
D514ken <- rexel(fn,sheet=3,skip = 1)

## LSO
fn <- LP+"lesotho/TB Diagnosed up to 31Dec2020_Les.xlsx"
## dx014
D014lso <- rexel(fn,sheet=1,skip = 1)
## dx04
D04lso <- rexel(fn,sheet=2,skip = 1)
## dx514
D514lso <- rexel(fn,sheet=3,skip = 1)

## DRC
fn <- LP+"malawianddrcdata/TB Diagnosed up to 31Dec2020_DRC.xlsx"
## dx014
D014drc <- rexel(fn,sheet=1,skip = 1)
## dx04
D04drc <- rexel(fn,sheet=2,skip = 1)
## dx514
D514drc <- rexel(fn,sheet=3,skip = 1)

## MWI
fn <- LP+"malawianddrcdata/TB Diagnosed up to 31Dec2020_Malawi.xlsx"
## dx014
D014mwi <- rexel(fn,sheet=1,skip = 1)
## dx04
D04mwi <- rexel(fn,sheet=2,skip = 1)
## dx514
D514mwi <- rexel(fn,sheet=3,skip = 1)

## CIV
fn <- LP+"malawianddrcdata/TB Diagnosed up to 31Dec2020_CDI.xlsx"
## dx014
D014civ <- rexel(fn,sheet=1,skip = 1)
## dx04
D04civ <- rexel(fn,sheet=2,skip = 1)
## dx514
D514civ <- rexel(fn,sheet=3,skip = 1)

## UGA
fn <- LP+"malawianddrcdata/TB Diagnosed up to 31Dec2020_Ug.xlsx"
## dx014
D014uga <- rexel(fn,sheet=1,skip = 1)
## dx04
D04uga <- rexel(fn,sheet=2,skip = 1)
## dx514
D514uga <- rexel(fn,sheet=3,skip = 1)

## ZWE
fn <- LP+"malawianddrcdata/TB Diagnosed up to 31Dec2020_Zim.xlsx"
## dx014
D014zwe <- rexel(fn,sheet=1,skip = 1)
## dx04
D04zwe <- rexel(fn,sheet=2,skip = 1)
## dx514
D514zwe <- rexel(fn,sheet=3,skip = 1)

## TZA NOTE this is interim data
fn <- LP+"../JeffDatacheck/interimdata/TB Diagnosed up to 30June2020.xlsx"
## dx014
D014tza <- rexel(fn,sheet=1,skip = 1)[`Facility Hierarchy - Country`=='Tanzania']
## dx04
D04tza <- rexel(fn,sheet=2,skip = 1)[`Facility Hierarchy - Country`=='Tanzania']
## dx514
D514tza <- rexel(fn,sheet=3,skip = 1)[`Facility Hierarchy - Country`=='Tanzania']

## TODO remove with final data
## renaming for consistency with final data
D04tza <- D04tza[,.(`Facility Hierarchy - Country`,
                    `Facility Hierarchy - Facility Name`,
                    `Num Diagnosed - Baseline 2`=`Num Diagnosed - Baseline`,
                    `Num Months - Baseline`=`Num Diagnosed - Baseline`/
                      `average monthly rate - Baseline`,
                    `Baseline TB DX Monthly Rate 2`=`average monthly rate - Baseline`,
                    `Num Diagnosed - Intervention 2`=`Num Diagnosed - Intervention`,
                    `Num Months - Intervention`=`Num Diagnosed - Intervention`/
                      `average monthly rate - Intervention`,
                    `Intervention TB Diagnosed Average Monthly Rate`=`average monthly rate - Intervention`,
                    `Fold Increase TB DX`,
                    `% improvement`=NA)]
D014tza <- D014tza[,.(`Facility Hierarchy - Country`,
                    `Facility Hierarchy - Facility Name`,
                    `Num Diagnosed - Baseline 2`=`Num Diagnosed - Baseline`,
                    `Num Months - Baseline`=`Num Diagnosed - Baseline`/
                      `average monthly rate - Baseline`,
                    `Baseline TB DX Monthly Rate 2`=`average monthly rate - Baseline`,
                    `Num Diagnosed - Intervention 2`=`Num Diagnosed - Intervention`,
                    `Num Months - Intervention`=`Num Diagnosed - Intervention`/
                      `average monthly rate - Intervention`,
                    `Intervention TB Diagnosed Average Monthly Rate`=`average monthly rate - Intervention`,
                    `Fold Increase TB DX`,
                    `% improvement`=NA)]
D514tza <- D514tza[,.(`Facility Hierarchy - Country`,
                    `Facility Hierarchy - Facility Name`,
                    `Num Diagnosed - Baseline 2`=`Num Diagnosed - Baseline`,
                    `Num Months - Baseline`=`Num Diagnosed - Baseline`/
                      `average monthly rate - Baseline`,
                    `Baseline TB DX Monthly Rate 2`=`average monthly rate - Baseline`,
                    `Num Diagnosed - Intervention 2`=`Num Diagnosed - Intervention`,
                    `Num Months - Intervention`=`Num Diagnosed - Intervention`/
                      `average monthly rate - Intervention`,
                    `Intervention TB Diagnosed Average Monthly Rate`=`average monthly rate - Intervention`,
                    `Fold Increase TB DX`,
                    `% improvement`=NA)]
D04tza[!is.finite(`Num Months - Baseline`),`Num Months - Baseline`:=12]
D014tza[!is.finite(`Num Months - Baseline`),`Num Months - Baseline`:=12]
D514tza[!is.finite(`Num Months - Baseline`),`Num Months - Baseline`:=12]

## renaming etc
D04 <- droprename2(rbindlist(list(D04civ,
                                  D04cmr,
                                  D04drc,
                                  D04ken,
                                  D04lso,
                                  D04mwi,
                                  D04uga,
                                  D04tza,
                                  D04zwe)))
D014 <- droprename2(rbindlist(list(D014civ,
                                  D014cmr,
                                  D014drc,
                                  D014ken,
                                  D014lso,
                                  D014mwi,
                                  D014uga,
                                  D014tza,
                                  D014zwe)))
D514 <- droprename2(rbindlist(list(D514civ,
                                   D514cmr,
                                   D514drc,
                                   D514ken,
                                   D514lso,
                                   D514mwi,
                                   D514uga,
                                   D514tza,
                                   D514zwe)))

## drop the prescribed missing data facilities
D04 <- D04[!Facility %in% dtdrop]
D514 <- D514[!Facility %in% dtdrop]
D014 <- D014[!Facility %in% dtdrop]

## 0-4
DM04 <- dataprep(D04)
DIP <- intplot(DM04)
ggsave(DIP,file=here('graphs/D04.pdf'),h=8.5,w=7.5)
## ggsave(DIP,file=here('graphs/D04.png'),h=8.5,w=7.5)

## 5-14
DM514 <- dataprep(D514)
DIP <- intplot(DM514)
ggsave(DIP,file=here('graphs/D514.pdf'),h=8.5,w=7.5)
## ggsave(DIP,file=here('graphs/D514.png'),h=8.5,w=7.5)

## 0-14
DM014 <- dataprep(D014)
DIP <- intplot(DM014)
ggsave(DIP,file=here('graphs/D014.pdf'),h=8.5,w=7.5)
## ggsave(DIP,file=here('graphs/D014.png'),h=8.5,w=7.5)


## ================= tx init ==================

## CMR
fn <- LP+"cameroon/TB Treatment up to 31Dec2020_Cam.xlsx"
## dx014
T014cmr <- rexel(fn,sheet=1,skip = 1)
## dx04
T04cmr <- rexel(fn,sheet=2,skip = 1)
## dx514
T514cmr <- rexel(fn,sheet=3,skip = 1)

## KEN
fn <- LP+"kenya/TB Treatment up to 31Dec2020_Ken.xlsx"
## dx014
T014ken <- rexel(fn,sheet=1,skip = 1)
## dx04
T04ken <- rexel(fn,sheet=2,skip = 1)
## dx514
T514ken <- rexel(fn,sheet=3,skip = 1)

## LSO
fn <- LP+"lesotho/TB Treatment up to 31Dec2020_Les.xlsx"
## dx014
T014lso <- rexel(fn,sheet=1,skip = 1)
## dx04
T04lso <- rexel(fn,sheet=2,skip = 1)
## dx514
T514lso <- rexel(fn,sheet=3,skip = 1)

## DRC
fn <- LP+"malawianddrcdata/TB Treatment up to 31Dec2020_DRC.xlsx"
## tx014
T014drc <- rexel(fn,sheet=1,skip = 1)
## tx04
T04drc <- rexel(fn,sheet=2,skip = 1)
## tx514
T514drc <- rexel(fn,sheet=3,skip = 1)

## MWI
fn <- LP+"malawianddrcdata/TB Treatment up to 31Dec2020_Malawi.xlsx"
## tx014
T014mwi <- rexel(fn,sheet=1,skip = 1)
## tx04
T04mwi <- rexel(fn,sheet=2,skip = 1)
## tx514
T514mwi <- rexel(fn,sheet=3,skip = 1)

## CIV
fn <- LP+"malawianddrcdata/TB Treatment up to 31Dec2020_CDI.xlsx"
## dx014
T014civ <- rexel(fn,sheet=1,skip = 1)
## dx04
T04civ <- rexel(fn,sheet=2,skip = 1)
## dx514
T514civ <- rexel(fn,sheet=3,skip = 1)

## UGA
fn <- LP+"malawianddrcdata/TB Treatment up to 31Dec2020_Ug.xlsx"
## dx014
T014uga <- rexel(fn,sheet=1,skip = 1)
## dx04
T04uga <- rexel(fn,sheet=2,skip = 1)
## dx514
T514uga <- rexel(fn,sheet=3,skip = 1)

## ZWE
fn <- LP+"malawianddrcdata/TB Treatment up to 31Dec2020_Zim.xlsx"
## dx014
T014zwe <- rexel(fn,sheet=1,skip = 1)
## dx04
T04zwe <- rexel(fn,sheet=2,skip = 1)
## dx514
T514zwe <- rexel(fn,sheet=3,skip = 1)

## TZA NOTE interim data
fn <- LP+"../JeffDatacheck/interimdata/TB TX Initiation up to 30June2020.xlsx"
## dx014
T014tza <- rexel(fn,sheet=1,skip = 1)[`Facility Hierarchy - Country`=='Tanzania']
## dx04
T04tza <- rexel(fn,sheet=2,skip = 1)[`Facility Hierarchy - Country`=='Tanzania']
## dx514
T514tza <- rexel(fn,sheet=3,skip = 1)[`Facility Hierarchy - Country`=='Tanzania']

## TODO remove with final data
## renaming for consistency with final data
T04tza <- T04tza[,.(`Facility Hierarchy - Country`,
                    `Facility Hierarchy - Facility Name`,
                    `Num TB TX Initiated - Baseline`,
                    `Num Months - Baseline`=`Num TB TX Initiated - Baseline`/
                      `Baseline TB TX Init Average Monthly Rate`,
                    `Baseline TB TX Init Average Monthly Rate`=`Baseline TB TX Init Average Monthly Rate`,
                    `Num TB TX Initiated - Intervention`,
                    `Num Months - Intervention`=`Num TB TX Initiated - Intervention`/
                      `Intervention TB Tx Init Average monthly rate`,
                    `Intervention TB TX Init Average Monthly Rate`=
                      `Intervention TB Tx Init Average monthly rate`,
                    `Fold Increase TB TX Init`)]
T014tza <- T014tza[,.(`Facility Hierarchy - Country`,
                      `Facility Hierarchy - Facility Name`,
                      `Num TB TX Initiated - Baseline`,
                      `Num Months - Baseline`=`Num TB TX Initiated - Baseline`/
                        `Baseline TB TX Init Average Monthly Rate`,
                      `Baseline TB TX Init Average Monthly Rate`=`Baseline TB TX Init Average Monthly Rate`,
                      `Num TB TX Initiated - Intervention`,
                      `Num Months - Intervention`=`Num TB TX Initiated - Intervention`/
                        `Intervention TB Tx Init Average monthly rate`,
                      `Intervention TB TX Init Average Monthly Rate`=
                        `Intervention TB Tx Init Average monthly rate`,
                      `Fold Increase TB TX Init`)]
T514tza <- T514tza[,.(`Facility Hierarchy - Country`,
                      `Facility Hierarchy - Facility Name`,
                      `Num TB TX Initiated - Baseline`,
                      `Num Months - Baseline`=`Num TB TX Initiated - Baseline`/
                        `Baseline TB TX Init Average Monthly Rate`,
                      `Baseline TB TX Init Average Monthly Rate`=`Baseline TB TX Init Average Monthly Rate`,
                      `Num TB TX Initiated - Intervention`,
                      `Num Months - Intervention`=`Num TB TX Initiated - Intervention`/
                        `Intervention TB Tx Init Average monthly rate`,
                      `Intervention TB TX Init Average Monthly Rate`=
                        `Intervention TB Tx Init Average monthly rate`,
                      `Fold Increase TB TX Init`)]
T04tza[!is.finite(`Num Months - Baseline`),`Num Months - Baseline`:=12]
T014tza[!is.finite(`Num Months - Baseline`),`Num Months - Baseline`:=12]
T514tza[!is.finite(`Num Months - Baseline`),`Num Months - Baseline`:=12]

names(T04drc); names(T04civ); names(T04mwi); names(T04uga)
names(T04cmr); names(T04lso); names(T04ken); names(T04tza);

## fix needed
T04drc[,`% improvement`:=NULL]
T04mwi[,`% improvement`:=NULL]
T04zwe[,`% improvement`:=NULL]
T04uga[,`% improvement`:=NULL]
T04cmr[,`% improvement`:=NULL]
T04lso[,`% improvement`:=NULL]
T04ken[,`% improvement`:=NULL]
T04civ[,`% improvement`:=NULL]
## NOTE may need for TZA when using final data
T514drc[,`% improvement`:=NULL]
T514mwi[,`% improvement`:=NULL]
T514zwe[,`% improvement`:=NULL]
T514uga[,`% improvement`:=NULL]
T514cmr[,`% improvement`:=NULL]
T514lso[,`% improvement`:=NULL]
T514ken[,`% improvement`:=NULL]
T514civ[,`% improvement`:=NULL]
## last age gp
T014drc[,`% improvement`:=NULL]
T014mwi[,`% improvement`:=NULL]
T014zwe[,`% improvement`:=NULL]
T014uga[,`% improvement`:=NULL]
T014cmr[,`% improvement`:=NULL]
T014lso[,`% improvement`:=NULL]
T014ken[,`% improvement`:=NULL]
T014civ[,`% improvement`:=NULL]

## renaming etc
T04 <- droprename2(rbindlist(list(T04civ,
                                  T04cmr,
                                  T04drc,
                                  T04ken,
                                  T04lso,
                                  T04mwi,
                                  T04uga,
                                  T04tza,
                                  T04zwe)))
T514 <- droprename2(rbindlist(list(T514civ,
                                   T514cmr,
                                   T514drc,
                                   T514ken,
                                   T514lso,
                                   T514mwi,
                                   T514uga,
                                   T514tza,
                                   T514zwe)))
T014 <- droprename2(rbindlist(list(T014civ,
                                   T014cmr,
                                   T014drc,
                                   T014ken,
                                   T014lso,
                                   T014mwi,
                                   T014uga,
                                   T014tza,
                                   T014zwe)))

## drop the prescribed missing data facilities
T04 <- T04[!Facility %in% dtdrop]
T514 <- T514[!Facility %in% dtdrop]
T014 <- T014[!Facility %in% dtdrop]

## 0-4
TM04 <- dataprep(T04)
TIP <- intplot(TM04)
ggsave(TIP,file=here('graphs/T04.pdf'),h=8.5,w=7.5)
## ggsave(TIP,file=here('graphs/T04.png'),h=8.5,w=7.5)

## 5-14
TM514 <- dataprep(T514)
TIP <- intplot(TM514)
ggsave(TIP,file=here('graphs/T514.pdf'),h=8.5,w=7.5)
## ggsave(TIP,file=here('graphs/T514.png'),h=8.5,w=7.5)

## 0-14
TM014 <- dataprep(T014)
TIP <- intplot(TM014)
ggsave(TIP,file=here('graphs/T014.pdf'),h=8.5,w=7.5)
## ggsave(TIP,file=here('graphs/T014.png'),h=8.5,w=7.5)

## ================= pt init ==================

## CMR
fn <- LP+"cameroon/TPT initations up to 31Dec2020_Cam.xlsx"
## dx014
P014cmr <- rexel(fn,sheet=1,skip = 1)
## dx04
P04cmr <- rexel(fn,sheet=2,skip = 1)
## dx514
P514cmr <- rexel(fn,sheet=3,skip = 1)

## KEN
fn <- LP+"kenya/TPT initations up to 31Dec2020_Ken.xlsx"
## dx014
P014ken <- rexel(fn,sheet=1,skip = 1)
## dx04
P04ken <- rexel(fn,sheet=2,skip = 1)
## dx514
P514ken <- rexel(fn,sheet=3,skip = 1)

## LSO
fn <- LP+"lesotho/TPT initations up to 31Dec2020_Les.xlsx"
## dx014
P014lso <- rexel(fn,sheet=1,skip = 1)
## dx04
P04lso <- rexel(fn,sheet=2,skip = 1)
## dx514
P514lso <- rexel(fn,sheet=3,skip = 1)

## DRC
fn <- LP+"malawianddrcdata/TPT initations up to 31Dec2020_DRC.xlsx"
## tx014
P014drc <- rexel(fn,sheet=1,skip = 1)
## tx04
P04drc <- rexel(fn,sheet=2,skip = 1)
## tx514
P514drc <- rexel(fn,sheet=3,skip = 1)

## MWI
fn <- LP+"malawianddrcdata/TPT initations up to 31Dec2020_Malawi.xlsx"
## tx014
P014mwi <- rexel(fn,sheet=1,skip = 1)
## tx04
P04mwi <- rexel(fn,sheet=2,skip = 1)
## tx514
P514mwi <- rexel(fn,sheet=3,skip = 1)

## CIV
fn <- LP+"malawianddrcdata/TPT initations up to 31Dec2020_CDI.xlsx"
## dx014
P014civ <- rexel(fn,sheet=1,skip = 1)
## dx04
P04civ <- rexel(fn,sheet=2,skip = 1)
## dx514
P514civ <- rexel(fn,sheet=3,skip = 1)

## UGA
fn <- LP+"malawianddrcdata/TPT initations up to 31Dec2020_Ug.xlsx"
## dx014
P014uga <- rexel(fn,sheet=1,skip = 1)
## dx04
P04uga <- rexel(fn,sheet=2,skip = 1)
## dx514
P514uga <- rexel(fn,sheet=3,skip = 1)

## ZWE
fn <- LP+"malawianddrcdata/TPT initations up to 31Dec2020_Zim.xlsx"
## dx014
P014zwe <- rexel(fn,sheet=1,skip = 1)
## dx04
P04zwe <- rexel(fn,sheet=2,skip = 1)
## dx514
P514zwe <- rexel(fn,sheet=3,skip = 1)

## TZA NOTE interim data
fn <- LP+"../JeffDatacheck/interimdata/TPT Initiation up to 30June2020.xlsx"
## dx014
P014tza <- rexel(fn,sheet=1,skip = 1)[`Facility Hierarchy - Country`=='Tanzania']
## dx04
P04tza <- rexel(fn,sheet=2,skip = 1)[`Facility Hierarchy - Country`=='Tanzania']
## dx514
P514tza <- rexel(fn,sheet=3,skip = 1)[`Facility Hierarchy - Country`=='Tanzania']

## TODO remove with final data
## renaming for consistency with final data
P04tza <- P04tza[,.(`Facility Hierarchy - Country`,
                    `Facility Hierarchy - Facility Name`,
                    `Num PT Initiated - Baseline`,
                    `Num Months - Baseline`=`Num PT Initiated - Baseline`/
                      `average monthly rate - Baseline`,
                    `Baseline PT Init Average Monthly Rate`=
                      `average monthly rate - Baseline`,
                    `Num PT Initiated - Intervention`,
                    `Num Months - Intervention`=`Num PT Initiated - Intervention`/
                      `average monthly rate - Intervention`,
                    `Intervention PT Init Average Monthly Rate`=
                      `average monthly rate - Intervention`,
                    `Fold Increase PT Init`,
                    `% improvement`=`Pct Improvement`)]
P014tza <- P014tza[,.(`Facility Hierarchy - Country`,
                      `Facility Hierarchy - Facility Name`,
                      `Num PT Initiated - Baseline`,
                      `Num Months - Baseline`=`Num PT Initiated - Baseline`/
                        `average monthly rate - Baseline`,
                      `Baseline PT Init Average Monthly Rate`=
                        `average monthly rate - Baseline`,
                      `Num PT Initiated - Intervention`,
                      `Num Months - Intervention`=`Num PT Initiated - Intervention`/
                        `average monthly rate - Intervention`,
                      `Intervention PT Init Average Monthly Rate`=
                        `average monthly rate - Intervention`,
                      `Fold Increase PT Init`,
                      `% improvement`=`Pct Improvement`)]
P514tza <- P514tza[,.(`Facility Hierarchy - Country`,
                      `Facility Hierarchy - Facility Name`,
                      `Num PT Initiated - Baseline`,
                      `Num Months - Baseline`=`Num PT Initiated - Baseline`/
                        `average monthly rate - Baseline`,
                      `Baseline PT Init Average Monthly Rate`=
                        `average monthly rate - Baseline`,
                      `Num PT Initiated - Intervention`,
                      `Num Months - Intervention`=`Num PT Initiated - Intervention`/
                        `average monthly rate - Intervention`,
                      `Intervention PT Init Average Monthly Rate`=
                        `average monthly rate - Intervention`,
                      `Fold Increase PT Init`,
                      `% improvement`=`Pct Improvement`)]
P04tza[!is.finite(`Num Months - Baseline`),`Num Months - Baseline`:=12]
P014tza[!is.finite(`Num Months - Baseline`),`Num Months - Baseline`:=12]
P514tza[!is.finite(`Num Months - Baseline`),`Num Months - Baseline`:=12]

## renaming etc
P04 <- droprename2(rbindlist(list(P04civ,
                                  P04cmr,
                                  P04drc,
                                  P04ken,
                                  P04lso,
                                  P04mwi,
                                  P04uga,
                                  P04tza,
                                  P04zwe)))
P514 <- droprename2(rbindlist(list(P514civ,
                                   P514cmr,
                                   P514drc,
                                   P514ken,
                                   P514lso,
                                   P514mwi,
                                   P514uga,
                                   P514tza,
                                   P514zwe)))
P014 <- droprename2(rbindlist(list(P014civ,
                                   P014cmr,
                                   P014drc,
                                   P014ken,
                                   P014lso,
                                   P014mwi,
                                   P014uga,
                                   P014tza,
                                   P014zwe)))

## drop the prescribed missing data facilities
P04 <- P04[!Facility %in% ptdrop]
P514 <- P514[!Facility %in% ptdrop]
P014 <- P014[!Facility %in% ptdrop]

## 0-4
PM04 <- dataprep(P04)
PIP <- intplot(PM04)
ggsave(PIP,file=here('graphs/P04.pdf'),h=8.5,w=7.5)
## ggsave(PIP,file=here('graphs/P04.png'),h=8.5,w=7.5)

## 5-14
PM514 <- dataprep(P514)
PIP <- intplot(PM514)
ggsave(PIP,file=here('graphs/P514.pdf'),h=8.5,w=7.5)
## ggsave(PIP,file=here('graphs/P514.png'),h=8.5,w=7.5)

## 0-14
PM014 <- dataprep(P014)
PIP <- intplot(PM014)
ggsave(PIP,file=here('graphs/P014.pdf'),h=8.5,w=7.5)
## ggsave(PIP,file=here('graphs/P014.png'),h=8.5,w=7.5)

## ==== join and write out ============
## --- dx
DL <- list(copy(D04),copy(D514),copy(D014))
for(i in 1:3) {
    DL[[i]][,Facility:=factor(Facility)]
    DL[[i]][,site:=as.integer(Facility)]
    DL[[i]][,Facility:=NULL]
}
DL[[1]][,age:='0-4']; DL[[2]][,age:='5-14']; DL[[3]][,age:='0-14']
D <- rbindlist(DL)

fwrite(D,file=here('outdata/D.csv'))

## --- tx
TL <- list(copy(T04),copy(T514),copy(T014))
for(i in 1:3) {
    TL[[i]][,Facility:=factor(Facility)]
    TL[[i]][,site:=as.integer(Facility)]
    TL[[i]][,Facility:=NULL]
}
TL[[1]][,age:='0-4']; TL[[2]][,age:='5-14']; TL[[3]][,age:='0-14']
T <- rbindlist(TL)

fwrite(T,file=here('outdata/T.csv'))

## --- pt
PL <- list(copy(P04),copy(P514),copy(P014))
for(i in 1:3) {
    PL[[i]][,Facility:=factor(Facility)]
    PL[[i]][,site:=as.integer(Facility)]
    PL[[i]][,Facility:=NULL]
}
PL[[1]][,age:='0-4']; PL[[2]][,age:='5-14']; PL[[3]][,age:='0-14']
P <- rbindlist(PL)

fwrite(P,file=here('outdata/P.csv'))


## ================== descriptive statistics ===========

miqr <- function(x) paste0(round(median(x),2),' (',
                           round(quantile(x,.25),2)," - ",
                           round(quantile(x,.75),2),")")

nmz <- c("Number of sites",
         "Site time",
         "Median sites per country",
         "Median IQR site time",
         "Site time for BL",
         "Site time for INT",
         "Rate BL",
         "Rate INT",
         "Median IQR Rate BL",
         "Median IQR Rate INT",
         "Range country rates BL",
         "Range country rates INT",
         "Number of txs 0-4",
         "Number of txs 5-14",
         "Number of txs 0-14",
         "Number of txs 0-4",
         "Number of txs 5-14",
         "Number of txs 0-14",
         "Percent sites improving 0-4",
         "Percent sites improving 5-14",
         "Percent sites improving 0-14")

## ---- Dx
Dvls <- c(length(unique(c(D04$Facility,D014$Facility,D514$Facility))),
         round(sum(D$Intervention.FT)+sum(D$Baseline.FT),2),
         NA,## "Median sites per country",
         miqr(D[,Intervention.FT+Baseline.FT]),## "Median IQR site time",
         round(sum(D$Baseline.FT),2),## "Site time for BL",
         round(sum(D$Intervention.FT),2),## "Site time for INT",
         round(sum(D$Baseline.Num)/sum(D$Baseline.FT),2),## "Rate BL",
         round(sum(D$Intervention.Num)/sum(D$Intervention.FT),2),## "Rate INT",
         miqr(D[,Baseline.Rate]),## "Median IQR Rate BL",
         miqr(D[,Baseline.Rate]),## "Median IQR Rate INT",
         paste(round(range(D[,.(rate=sum(Baseline.Num)/sum(Baseline.FT)),
                             by=Country]$rate),2),collapse = ' - '),## "Range country rates BL",
         paste(round(range(D[,.(rate=sum(Intervention.Num)/sum(Intervention.FT)),
                             by=Country]$rate),2),collapse = ' - '),## "Range country rates INT",
         NA,## "Number of txs 0-4",
         NA,## "Number of txs 5-14",
         NA,## "Number of txs 0-14",
         NA,## "Number of txs 0-4",
         NA,## "Number of txs 5-14",
         NA,## "Number of txs 0-14", TODO?
         D[age=='0-4',round(1e2*mean(Baseline.Rate<Intervention.Rate),2)],## "Percent sites improving 0-4",
         D[age=='5-14',round(1e2*mean(Baseline.Rate<Intervention.Rate),2)],## "Percent sites improving 5-14",
         D[age=='0-14',round(1e2*mean(Baseline.Rate<Intervention.Rate),2)]## "Percent sites improving 0-14"
         )
Dvls

Dstats <- data.table(name=nmz,value=Dvls)
fwrite(Dstats,file=here('outdata/Dstats.csv'))


## ---- Tx
Tvls <- c(length(unique(c(T04$Facility,T014$Facility,T514$Facility))),
         round(sum(T$Intervention.FT)+sum(T$Baseline.FT),2),
         NA,## "Median sites per country",
         miqr(T[,Intervention.FT+Baseline.FT]),## "Median IQR site time",
         round(sum(T$Baseline.FT),2),## "Site time for BL",
         round(sum(T$Intervention.FT),2),## "Site time for INT",
         round(sum(T$Baseline.Num)/sum(T$Baseline.FT),2),## "Rate BL",
         round(sum(T$Intervention.Num)/sum(T$Intervention.FT),2),## "Rate INT",
         miqr(T[,Baseline.Rate]),## "Median IQR Rate BL",
         miqr(T[,Baseline.Rate]),## "Median IQR Rate INT",
         paste(round(range(T[,.(rate=sum(Baseline.Num)/sum(Baseline.FT)),by=Country]$rate),2),collapse = ' - '),## "Range country rates BL",
         paste(round(range(T[,.(rate=sum(Intervention.Num)/sum(Intervention.FT)),by=Country]$rate),2),collapse = ' - '),## "Range country rates INT",
         NA,## "Number of txs 0-4",
         NA,## "Number of txs 5-14",
         NA,## "Number of txs 0-14",
         NA,## "Number of txs 0-4",
         NA,## "Number of txs 5-14",
         NA,## "Number of txs 0-14", TODO?
         T[age=='0-4',round(1e2*mean(Baseline.Rate<Intervention.Rate),2)],## "Percent sites improving 0-4",
         T[age=='5-14',round(1e2*mean(Baseline.Rate<Intervention.Rate),2)],## "Percent sites improving 5-14",
         T[age=='0-14',round(1e2*mean(Baseline.Rate<Intervention.Rate),2)]## "Percent sites improving 0-14"
         )
Tvls


Tstats <- data.table(name=nmz,value=Tvls)
fwrite(Tstats,file=here('outdata/Tstats.csv'))

## ---- PT
Pvls <- c(length(unique(c(P04$Facility,P014$Facility,P514$Facility))),
         round(sum(P$Intervention.FT)+sum(P$Baseline.FT),2),
         NA,## "Median sites per country",
         miqr(P[,Intervention.FT+Baseline.FT]),## "Median IQR site time",
         round(sum(P$Baseline.FT),2),## "Site time for BL",
         round(sum(P$Intervention.FT),2),## "Site time for INT",
         round(sum(P$Baseline.Num)/sum(P$Baseline.FT),2),## "Rate BL",
         round(sum(P$Intervention.Num)/sum(P$Intervention.FT),2),## "Rate INT",
         miqr(P[,Baseline.Rate]),## "Median IQR Rate BL",
         miqr(P[,Baseline.Rate]),## "Median IQR Rate INT",
         paste(round(range(P[,.(rate=sum(Baseline.Num)/sum(Baseline.FT)),by=Country]$rate),2),collapse = ' - '),## "Range country rates BL",
         paste(round(range(P[,.(rate=sum(Intervention.Num)/sum(Intervention.FT)),by=Country]$rate),2),collapse = ' - '),## "Range country rates INT",
         NA,## "Number of txs 0-4",
         NA,## "Number of txs 5-14",
         NA,## "Number of txs 0-14",
         NA,## "Number of txs 0-4",
         NA,## "Number of txs 5-14",
         NA,## "Number of txs 0-14", TODO?
         P[age=='0-4',round(1e2*mean(Baseline.Rate<Intervention.Rate),2)],## "Percent sites improving 0-4",
         P[age=='5-14',round(1e2*mean(Baseline.Rate<Intervention.Rate),2)],## "Percent sites improving 5-14",
         P[age=='0-14',round(1e2*mean(Baseline.Rate<Intervention.Rate),2)]## "Percent sites improving 0-14"
         )
Pvls


Pstats <- data.table(name=nmz,value=Pvls)
fwrite(Pstats,file=here('outdata/Pstats.csv'))

