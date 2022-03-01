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
upload.to.sheets(here('model/outdata/'),"ACFcascade.1.csv",shid)
upload.to.sheets(here('model/outdata/'),'ICERatt',shid)
upload.to.sheets(here('model/outdata/'),'ICERSatt',shid)
upload.to.sheets(here('model/outdata/'),"ICERagept.1.csv",shid)
upload.to.sheets(here('model/outdata/'),"ICERall.1.csv",shid)
upload.to.sheets(here('model/outdata/'),"ICERpt.1.csv",shid)
upload.to.sheets(here('model/outdata/'),'PTC',shid)
upload.to.sheets(here('model/outdata/'),"CEAC50.csv",shid)
upload.to.sheets(here('model/outdata/'),"CEAC50pt.1.csv",shid)
upload.to.sheets(here('model/outdata/'),"ptsuccess.csv",shid)
upload.to.sheets(here('model/outdata/'),"txsuccess.csv",shid)


## need article tables
yurl <- "https://docs.google.com/spreadsheets/d/1p8ZT0BP-lABM0W0z6V6I_4cg51ndPYyyUKn11Cq9f5M/edit#gid=0"
shidneat <- as.character(as_sheets_id(yurl))

## ---- Table 1 -------
## build 1st
load(here('model/data/Table1ATT.Rdata'))
load(here('model/data/Table1PT.Rdata'))
load(here('model/data/Table1PTcost.Rdata'))
## TODO load cost


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
## TODO
