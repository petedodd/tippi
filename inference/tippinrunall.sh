#!/bin/bash
# stan inference for each age, PT and ATT
# NOTE now parallelized - only use on decent machine in this form
R --vanilla < tippiMA.R px 1 & R --vanilla < tippiMA.R px 2 & R --vanilla < tippiMA.R px 3 & R --vanilla < tippiMA.R tx 1 & R --vanilla < tippiMA.R tx 2 & R --vanilla < tippiMA.R tx 3

# then some composite graphs/outputs
# R --vanilla < tippigather.R


# R --vanilla < tippistanrunner.R px 1
# R --vanilla < tippistanrunner.R px 2
# R --vanilla < tippistanrunner.R px 3
# R --vanilla < tippistanrunner.R tx 1
# R --vanilla < tippistanrunner.R tx 2
# R --vanilla < tippistanrunner.R tx 3

# # then some composite graphs/outputs
# R --vanilla < tippigather.R
