#!/bin/bash
# (after changing shell flag in modeloutcomes.R to TRUE)
# arg1: sensitivity analysis: none, hhc (10USD soc HHC cost), cdr (higher cdr for incidence), txd (completion of ATT/TPT included)
# arg2: ACF: 1/0
R --slave --vanilla --args <modeloutcomes.R none 1



