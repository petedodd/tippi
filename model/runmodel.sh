#!/bin/bash
# NOTE (after changing shell flag in modeloutcomes.R to TRUE)
# arg1: sensitivity analysis: none, base/lo/hi dscr, cdr (higher cdr for incidence), txd (completion of ATT/TPT included)
R --slave --vanilla --args <modeloutcomes.R hi & R --slave --vanilla --args <modeloutcomes.R lo & R --slave --vanilla --args <modeloutcomes.R cdr & R --slave --vanilla --args <modeloutcomes.R txd
R --slave --vanilla --args <modeloutcomes.R none



