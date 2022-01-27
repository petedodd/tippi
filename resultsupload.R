library(here)
library(data.table)
library(glue)
library(googlesheets4)


## NOTE new sheet
## create an ID to access the googlesheets results sheet
yourl <- "https://docs.google.com/spreadsheets/d/1yRUIErcA9H_1rbgG46WoKARt9-8t94tnJRS3X91gzFo/edit#gid=0"
shid <- as.character(as_sheets_id(yourl))


## utility function
upload.to.sheets <- function(basename,filename,sheetid,
                             sap #pulmonary sa?
                             ){
  fn <- glue(basename) + filename + ".csv"
  tmp <- fread(file=fn)
  write_sheet(tmp,sheetid,sheet=filename)
}

## upload relevant table data
upload.to.sheets(here('model/outdata/'),'cascadetab',shid)
upload.to.sheets(here('model/outdata/'),'ICERatt',shid)
upload.to.sheets(here('model/outdata/'),'ICERSatt',shid)
upload.to.sheets(here('model/outdata/'),'PTC',shid)




