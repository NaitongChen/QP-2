#####################################
# load libraries
#####################################
suppressWarnings(library(CVXR, warn.conflicts=FALSE))
library(glmnet)
library(doParallel)
library(foreach)
library(pense)

#####################################
# set seed and constants
#####################################
set.seed(2022)
n = 100;  # number of observations
p = 400;  # number of variables
CASE = 5;  # 5 scenarios for error distribution
LAMAX = 1.5;  # lambda_max
LAMAX_lasso = 10;
AMAX = 1.0;  # alpha_max
LENGTH_la = 15;  # number of lambdas used in validation
LENGTH_la_lasso = 15;
LENGTH_a = 15;  # number of alphas used in validation
N = 20;  # number of validation datasets
beta_0 = as.matrix(c(3*rep(1,20), rep(0,380)))  # true beta

# lamb = seq(0.0001,LAMAX,length.out = LENGTH_la);
# alp = seq(0.0001,AMAX,length.out = LENGTH_a);
# lamb_lasso = rev(seq(0.00001,LAMAX_lasso,length.out = LENGTH_la_lasso));
# lamb_pense = rev(seq(0.00001,LAMAX_lasso,length.out = LENGTH_la_lasso));

lamb = exp(seq(-10,2, length.out = LENGTH_la));
alp = exp(seq(-10,2, length.out = LENGTH_a));
lamb_lasso = rev(exp(seq(-15,2, length.out = LENGTH_la_lasso)));
lamb_pense = rev(exp(seq(-15,2, length.out = LENGTH_la_lasso)));

#####################################
# read outputs from simulation
#####################################
l2loss = read.csv("l2loss.csv")
l1loss = read.csv("l1loss.csv")
FP = read.csv("FP.csv")
FN = read.csv("FN.csv")
RA_1 = read.csv("RA_1.csv")
RA_2 = read.csv("RA_2.csv")

l2loss_pense = read.csv("l2loss_pense.csv")
l1loss_pense = read.csv("l1loss_pense.csv")
FP_pense = read.csv("FP_pense.csv")
FN_pense = read.csv("FN_pense.csv")
RA_1_pense = read.csv("RA_1_pense.csv")
RA_2_pense = read.csv("RA_2_pense.csv")

l2loss_lasso = read.csv("l2loss_lasso.csv")
l1loss_lasso = read.csv("l1loss_lasso.csv")
FP_lasso = read.csv("FP_lasso.csv")
FN_lasso = read.csv("FN_lasso.csv")

#####################################
# compute metrics
#####################################
for (i in 1:CASE) {
  round(mean(l1loss[,i]), 2)
  round(quantile(l1loss[,i], probs = c(0.05, 0.95)), 2)
  round(mean(l2loss[,i]), 2)
  round(quantile(l2loss[,i], probs = c(0.05, 0.95)), 2)
  round(mean(FP[,i]), 2)
  round(quantile(FP[,i], probs = c(0.05, 0.95)), 2)
  round(mean(FN[,i]), 2)
  round(quantile(FN[,i], probs = c(0.05, 0.95)), 2)
  round(mean(RA_1[,i]), 2)
  round(quantile(RA_1[,i], probs = c(0.05, 0.95)), 2)
  round(mean(RA_2[,i]), 2)
  round(quantile(RA_2[,i], probs = c(0.05, 0.95)), 2)
  
  round(mean(l1loss_lasso[,i]), 2)
  round(quantile(l1loss_lasso[,i], probs = c(0.05, 0.95)), 2)
  round(mean(l2loss_lasso[,i]), 2)
  round(quantile(l2loss_lasso[,i], probs = c(0.05, 0.95)), 2)
  round(mean(FP_lasso[,i]), 2)
  round(quantile(FP_lasso[,i], probs = c(0.05, 0.95)), 2)
  round(mean(FN_lasso[,i]), 2)
  round(quantile(FN_lasso[,i], probs = c(0.05, 0.95)), 2)
  
  round(mean(l1loss_pense[,i]), 2)
  round(quantile(l1loss_pense[,i], probs = c(0.05, 0.95)), 2)
  round(mean(l2loss_pense[,i]), 2)
  round(quantile(l2loss_pense[,i], probs = c(0.05, 0.95)), 2)
  round(mean(FP_pense[,i]), 2)
  round(quantile(FP_pense[,i], probs = c(0.05, 0.95)), 2)
  round(mean(FN_pense[,i]), 2)
  round(quantile(FN_pense[,i], probs = c(0.05, 0.95)), 2)
  round(mean(RA_1_pense[,i]), 2)
  round(quantile(RA_1_pense[,i], probs = c(0.05, 0.95)), 2)
  round(mean(RA_2_pense[,i]), 2)
  round(quantile(RA_2_pense[,i], probs = c(0.05, 0.95)), 2)
}