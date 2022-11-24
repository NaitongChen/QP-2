#####################################
# load libraries
#####################################
suppressWarnings(library(CVXR, warn.conflicts=FALSE))
library(glmnet)
library(doParallel)
library(foreach)

#####################################
# set seed and constants
#####################################
set.seed(2022)
n = 100;  # number of observations
p = 400;  # number of variables
CASE = 3;  # 3 scenarios for error distribution
LAMAX = 1.5;  # lambda_max
LAMAX_lasso = 10;
AMAX = 1.0;  # alpha_max
LENGTH_la = 30;  # number of lambdas used in validation
LENGTH_la_lasso = 200;
LENGTH_a = 15;  # number of alphas used in validation
N = 100;  # number of validation datasets
beta_0 = as.matrix(c(3*rep(1,20), rep(0,380)))  # true beta

lamb = seq(0.0001,LAMAX,length.out = LENGTH_la);
lamb_lasso = rev(seq(0.00001,LAMAX_lasso,length.out = LENGTH_la_lasso));
alp = seq(0.0001,AMAX,length.out = LENGTH_a);

#####################################
# read outputs from simulation
# and load actual data (not validation)
#####################################
hyperparams = read.csv('hyperparams.csv');
X = read.csv('x_ho.csv');
Y = read.csv('y_ho.csv');
Lambda.Alpha = hyperparams[,1:2]
Lambda = hyperparams[,3]

#####################################
# fit models using tuned hyper-param
#####################################
K=100;
CASE=3;
l2loss = matrix(0,K,CASE);
l1loss = matrix(0,K,CASE);
FP = matrix(0,K,CASE);
FN = matrix(0,K,CASE);
RA_1 = matrix(0,K,CASE);
RA_2 = matrix(0,K,CASE);

l2loss_lasso = matrix(0,K,CASE);
l1loss_lasso = matrix(0,K,CASE);
FP_lasso = matrix(0,K,CASE);
FN_lasso = matrix(0,K,CASE);


# loop through repetition
for (k in 1:K) {
  Betar = matrix(0,p,CASE);
  Beta_lasso = matrix(0,p,CASE);
  # loop through scenarios
  for (j in 1:CASE) {
    print(c(k,j))
    
    # RA lasso
    Betar[,j] = L1PenHuber(Y[(n*(k-1)+1):(n*k),j], X[(n*(k-1)+1):(n*k),], Lambda.Alpha[j,1], Lambda.Alpha[j,2]);
    Betar[,j] = Betar[,j] * (abs(Betar[,j]) > 1e-04);
    
    # lasso
    a = glmnet(x = X[(n*(k-1)+1):(n*k),], y = Y[(n*(k-1)+1):(n*k), j],
               family = "gaussian", alpha = 1, lambda = c(Lambda[j], 0), intercept = FALSE);
    Beta_lasso[,j] = a$beta[,1]
    Beta_lasso[,j] = Beta_lasso[,j] * (abs(Beta_lasso[,j]) > 1e-04);
    
    # compute metrics
    l2loss[k,j] = norm(Betar[,j] - beta_0, "2");
    l1loss[k,j] = norm(Betar[,j] - beta_0, "1");
    
    l2loss_lasso[k,j] = norm(Beta_lasso[,j] - beta_0, "2");
    l1loss_lasso[k,j] = norm(Beta_lasso[,j] - beta_0, "1");
    
    temp = Betar[,j];
    FP[k,j] = sum(temp[which(beta_0 == 0)] != 0);
    FN[k,j] = sum(temp[which(beta_0 != 0)] == 0);
    
    temp = Beta_lasso[,j];
    FP_lasso[k,j] = sum(temp[which(beta_0 == 0)] != 0);
    FN_lasso[k,j] = sum(temp[which(beta_0 != 0)] == 0);
    
    # oracle estimator for relative gain
    df = as.data.frame(cbind(Y[(n*(k-1)+1):(n*k),j], X[(n*(k-1)+1):(n*k),1:20]))
    colnames(df)[1] = "Y"
    fit = lm(formula = Y ~ . - 1, data = df)
    beta_oracle = fit$coefficients
    beta_oracle = c(beta_oracle, rep(0, 380))
    RA_1[k,j] = (norm(Beta_lasso[,j] - beta_0, "2") - norm(beta_oracle - beta_0, "2")) / (norm(Betar[,j] - beta_0, "2") - norm(beta_oracle - beta_0, "2"))
    RA_2[k,j] = (norm(Beta_lasso[,j] - beta_0, "1") - norm(beta_oracle - beta_0, "1")) / (norm(Betar[,j] - beta_0, "1") - norm(beta_oracle - beta_0, "1"))
  }
}

#####################################
# compute metric statistics
#####################################
mean(l1loss[,1])
quantile(l1loss[,1], probs = c(0.05, 0.95))
mean(l2loss[,1])
quantile(l2loss[,1], probs = c(0.05, 0.95))
mean(FP[,1]) # out of 400
quantile(FP[,1], probs = c(0.05, 0.95))
mean(FN[,1])
quantile(FN[,1], probs = c(0.05, 0.95))

mean(l1loss_lasso[,1])
quantile(l1loss_lasso[,1], probs = c(0.05, 0.95))
mean(l2loss_lasso[,1])
quantile(l2loss_lasso[,1], probs = c(0.05, 0.95))
mean(FP_lasso[,1])
quantile(FP_lasso[,1], probs = c(0.05, 0.95))
mean(FN_lasso[,1])
quantile(FN_lasso[,1], probs = c(0.05, 0.95))

mean(RA_1[,1])
quantile(RA_1[,1], probs = c(0.05, 0.95))
mean(RA_2[,1])
quantile(RA_2[,1], probs = c(0.05, 0.95))

mean(l1loss[,2])
quantile(l1loss[,2], probs = c(0.05, 0.95))
mean(l2loss[,2])
quantile(l2loss[,2], probs = c(0.05, 0.95))
mean(FP[,2])
quantile(FP[,2], probs = c(0.05, 0.95))
mean(FN[,2])
quantile(FN[,2], probs = c(0.05, 0.95))

mean(l1loss_lasso[,2])
quantile(l1loss_lasso[,2], probs = c(0.05, 0.95))
mean(l2loss_lasso[,2])
quantile(l2loss_lasso[,2], probs = c(0.05, 0.95))
mean(FP_lasso[,2])
quantile(FP_lasso[,2], probs = c(0.05, 0.95))
mean(FN_lasso[,2])
quantile(FN_lasso[,2], probs = c(0.05, 0.95))

mean(RA_1[,2])
quantile(RA_1[,2], probs = c(0.05, 0.95))
mean(RA_2[,2])
quantile(RA_2[,2], probs = c(0.05, 0.95))

mean(l1loss[,3])
quantile(l1loss[,3], probs = c(0.05, 0.95))
mean(l2loss[,3])
quantile(l2loss[,3], probs = c(0.05, 0.95))
mean(FP[,3])
quantile(FP[,3], probs = c(0.05, 0.95))
mean(FN[,3])
quantile(FN[,3], probs = c(0.05, 0.95))

mean(l1loss_lasso[,3])
quantile(l1loss_lasso[,3], probs = c(0.05, 0.95))
mean(l2loss_lasso[,3])
quantile(l2loss_lasso[,3], probs = c(0.05, 0.95))
mean(FP_lasso[,3])
quantile(FP_lasso[,3], probs = c(0.05, 0.95))
mean(FN_lasso[,3])
quantile(FN_lasso[,3], probs = c(0.05, 0.95))

mean(RA_1[,3])
quantile(RA_1[,3], probs = c(0.05, 0.95))
mean(RA_2[,3])
quantile(RA_2[,3], probs = c(0.05, 0.95))