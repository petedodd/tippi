#!/bin/bash
# stan inference for each age, PT and ATT
# NOTE now parallelized - only use on decent machine in this form
# R --vanilla < tippiMA.R px 1 & R --vanilla < tippiMA.R px 2 & R --vanilla < tippiMA.R px 3 & R --vanilla < tippiMA.R tx 1 & R --vanilla < tippiMA.R tx 2 & R --vanilla < tippiMA.R tx 3

# less parallel for laptop
R --vanilla < tippiMA.R px 1 & R --vanilla < tippiMA.R px 2
R --vanilla < tippiMA.R px 3 & R --vanilla < tippiMA.R tx 1
R --vanilla < tippiMA.R tx 2 & R --vanilla < tippiMA.R tx 3


# then some composite graphs/outputs
R --vanilla < tippigather.R


echo "done!"
