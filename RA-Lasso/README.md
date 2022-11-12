INTRODUCTION
------------
 * This .zip file contains two parts: (1) Code for simulation studies; (2) Dada and code for real data example. 

SIMULATION STUDIES
------------------
 * Code for simulation studies are under the folder “simulation”.

 * “generate_data.R” gives the R code for generating the validation and test datasets (.csv files). 

 * “main_simulation.m” gives the Matlab code for running our method.

REAL DATA EXAMPLE
-----------------
 * Data and code for the real data example are under the folder “realdata”.

 * Gene “TLR8” is used as response and stored in “y_real.csv”. The other 464 genes are used as covariates and stored in “x_real.csv”. 

 * The 3 methods we compared are written in 3 functions: “L1PenL2.m” (Lasso), “L1PenL1Single.m” (R-Lasso) and “L1PenHuber.m” (RA-Lasso).
