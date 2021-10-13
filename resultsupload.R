library(here)
library(data.table)
library(googlesheets4)



## create an ID to access the googlesheets results sheet
"https://docs.google.com/spreadsheets/d/1yRUIErcA9H_1rbgG46WoKARt9-8t94tnJRS3X91gzFo/edit#gid=0" %>%
  as_sheets_id() %>%
  as.character() -> shid


## utility functions
see <- function(x,ns=3)formatC(signif(x,ns),big.mark = ",",format='fg') #for reading big numbers
r1 <- function(x) sprintf("%.1f", round(x,1) )


## reformatting etc needed
cascadetab <- fread(file=here('model/outdata/cascadetab.csv'))
write_sheet(cascadetab,shid,sheet='test')


