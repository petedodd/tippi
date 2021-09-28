#!/bin/bash
# stan inference for each age, PT and ATT
R --vanilla < tippistanrunner.R px 1
R --vanilla < tippistanrunner.R px 2
R --vanilla < tippistanrunner.R px 3
R --vanilla < tippistanrunner.R tx 1
R --vanilla < tippistanrunner.R tx 2
R --vanilla < tippistanrunner.R tx 3

# then some composite graphs/outputs
R --vanilla < tippigather.R
