library(here)
library(data.table)
library(glue)
library(googlesheets4)


## NOTE new sheet
## create an ID to access the googlesheets results sheet
yourl <- "https://docs.google.com/spreadsheets/d/1yRUIErcA9H_1rbgG46WoKARt9-8t94tnJRS3X91gzFo/edit#gid=0"
shid <- as.character(as_sheets_id(yourl))


## utility function
upload.to.sheets <- function(basename,filename,sheetid
                             ){
  filename <- gsub("\\.csv$","",filename) #safety in case csv included at and
  fn <- glue(basename) + filename + ".csv"
  tmp <- fread(file=fn)
  write_sheet(tmp,sheetid,sheet=filename)
}

## upload relevant table data
upload.to.sheets(here('model/outdata/'),'cascadetab',shid) #first will need to re-authenticate

## rest can be run as block
flz1 <- c(
 "ACFcascade.1.csv",   "ACFcascadecdr.1.csv","ACFcascadehi.1.csv",
 "ACFcascadelo.1.csv", "ACFcascadetxd.1.csv","cascadetab.csv",
 "CEAC50.csv",         "CEAC50cdr.csv",      "CEAC50hi.csv",
 "CEAC50lo.csv",       "CEAC50pt.1.csv",     "CEAC50ptcdr.1.csv",
 "CEAC50pthi.1.csv",   "CEAC50ptlo.1.csv",   "CEAC50pttxd.1.csv",
 "CEAC50txd.csv",      "ICERagept.1.csv",    "ICERageptcdr.1.csv",
 "ICERagepthi.1.csv",  "ICERageptlo.1.csv",  "ICERagepttxd.1.csv",
 "ICERall.1.csv",      "ICERallcdr.1.csv",   "ICERallhi.1.csv",
 "ICERalllo.1.csv",    "ICERalltxd.1.csv")
for( fn in flz1)
  upload.to.sheets(here('model/outdata/'),fn,shid)

Sys.sleep(120) #wait a bit so as not to annoy google

flz2 <- c("ICERatt.csv",
 "ICERattcdr.csv",     "ICERatthi.csv",      "ICERattlo.csv",
 "ICERatttxd.csv",     "ICERpt.1.csv",       "ICERptcdr.1.csv",
 "ICERpthi.1.csv",     "ICERptlo.1.csv",     "ICERpttxd.1.csv",
 "ICERSatt.csv",       "ICERSattcdr.csv",    "ICERSatthi.csv",
 "ICERSattlo.csv",     "ICERSatttxd.csv",    "PTC.csv",
 "ptsuccess.csv",      "txsuccess.csv"
)
for( fn in flz2)
  upload.to.sheets(here('model/outdata/'),fn,shid)


## need article tables
yurl <- "https://docs.google.com/spreadsheets/d/1p8ZT0BP-lABM0W0z6V6I_4cg51ndPYyyUKn11Cq9f5M/edit#gid=0"
shidneat <- as.character(as_sheets_id(yurl))

## ---- Table 1 -------
## build 1st
load(here('model/data/Table1ATT.Rdata'))
load(here('model/data/Table1PT.Rdata'))
load(here('model/data/Table1PTcost.Rdata'))

setcolorder(Table1PTcost,names(Table1PT))

Table1 <- rbindlist(list(Table1ATT,Table1PT,Table1PTcost))

write_sheet(Table1,shidneat,sheet="Tab1RAW")

## ---- Table 2 -------

## ATT part
load(here('model/outdata/Table2ATT.Rdata'))
write_sheet(Table2ATT,shidneat,sheet="Tab2ATT")

## PT part
load(here('model/outdata/Table2PT.1.Rdata'))
write_sheet(Table2PT,shidneat,sheet="Tab2PT")

## combined intervention
load(here('model/outdata/Table2both.1.Rdata'))
write_sheet(Table2both,shidneat,sheet="Tab2combined")


## backgrounds stats for results
upload.to.sheets(here('dataprep/outdata/'),"Pstats.csv",shidneat)
upload.to.sheets(here('dataprep/outdata/'),"Dstats.csv",shidneat)
upload.to.sheets(here('dataprep/outdata/'),"Tstats.csv",shidneat)

