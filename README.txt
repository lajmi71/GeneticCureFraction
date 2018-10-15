

This is R code to perform an association test between a genetic regiron and a survival outcome with a cure fraction.

There are 4 files:
- Main.R : this is the main file. This file calls all other files.
- Sim.R : this file contains R code for the generation of genetic data and a survival outcome with a cure fraction.
- Score.R : this file contains R code to compute the test statistics and the p.values for the association tests developed in Lakhal-Chaieb et al. (2018).
- Functions.R : this file contains various R functions needed in Sim.R and in Score.R