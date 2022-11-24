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

Lambda.Alpha = matrix(0, CASE, 2);  # stores optimal pair of (lambda, alpha) for 4 scenarios
Lambda = matrix(0, CASE);

#####################################
# load data
#####################################
XV = read.csv('x_ho_v.csv');
YV = read.csv('y_ho_v.csv');

#####################################
# set up parallel computing
#####################################
n.cores <- parallel::detectCores() - 1;

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
);

doParallel::registerDoParallel(cl = my.cluster);

foreach::getDoParRegistered()
foreach::getDoParWorkers()

clusterExport(cl=my.cluster, list("n", "p", "CASE", "LAMAX", "LAMAX_lasso", "AMAX", 
                                  "LENGTH_la","LENGTH_la_lasso", "LENGTH_a", "N", "beta_0",
                                  "XV", "YV", "lamb", "lamb_lasso", "alp"),
              envir=environment())

#####################################
# helper: ra-lasso solver
#####################################
L1PenHuber <- function(Y, X, lambda, alpha) {
  np = dim(X);
  n = np[1]
  p = np[2]
  
  X = as.matrix(X)
  Y = as.matrix(Y)
  
  betaHat <- Variable(p,1)
  objective <- Minimize(sum(huber(Y - X %*% betaHat, 1/alpha)) + n*lambda*norm(betaHat,"1"))
  problem <- Problem(objective)
  result <- solve(problem, solver = "ECOS")
  
  return(result$getValue(betaHat))
}

#####################################
# helper: Validation for one scenario
#####################################
TuneHyperParam <- function(case) {
  
  print("RA lasso")
  loss1 = matrix(0, LENGTH_la, LENGTH_a);
  # for each of the validation set
  for (k in 1:N) {
    print(k)
    # compute validation error over 2d search grid
    loss1_ij = foreach (i = 1:LENGTH_la, .combine = 'cbind', .packages = "CVXR", .export = "L1PenHuber") %:%
      foreach (j = 1:LENGTH_a, .combine = 'c', .packages = "CVXR", .export ="L1PenHuber") %dopar% {
        betah = L1PenHuber(YV[(n*(k-1)+1):(n*k), case], XV[(n*(k-1)+1):(n*k),], lamb[i], alp[j]);
        betah = betah * (abs(betah) > 1e-04);
        betah = norm(betah-beta_0, "F");
      }
    loss1 = loss1 + t(loss1_ij)
  }
  
  print("RA lasso")
  loss1_lasso = matrix(0, LENGTH_la_lasso);
  # for each of the validation set
  for (k in 1:N) {
    print(k)
    a <- glmnet(x = XV[(n*(k-1)+1):(n*k),], y = YV[(n*(k-1)+1):(n*k), case],
                family = "gaussian", alpha = 1, lambda = lamb_lasso, intercept = FALSE)
    
    # compute validation error over search grid
    for (i in 1:LENGTH_la_lasso) {
      betah = a$beta[,i]
      betah = betah * (abs(betah) > 1e-04);
      loss1_lasso[i] = loss1_lasso[i] + norm(betah-beta_0, "F");
    }
  }
  
  # find index of optimal hyper parameter
  inds = which(loss1 == min(loss1), arr.ind = TRUE);
  q1 = inds[1]
  q2 = inds[2]
  
  ind = which(loss1_lasso == min(loss1_lasso))
  
  return(list(q1, q2, ind, loss1, loss1_lasso))
}

#####################################
# 1st scenario: N(0,4)
#####################################
inds = TuneHyperParam(1);
q1 = inds[[1]];
q2 = inds[[2]];
Lambda.Alpha[1,1] = lamb[q1];
Lambda.Alpha[1,2] = alp[q2];
Lambda[1] = lamb_lasso[inds[[3]]];

#####################################
# 2nd scenario: 2t_3
#####################################
inds = TuneHyperParam(2);
q1 = inds[[1]];
q2 = inds[[2]];
Lambda.Alpha[2,1] = lamb[q1];
Lambda.Alpha[2,2] = alp[q2];
Lambda[2] = lamb_lasso[inds[[3]]];

#####################################
# 3rd scenario: MixN
#####################################
inds = TuneHyperParam(3);
q1 = inds[[1]];
q2 = inds[[2]];
Lambda.Alpha[3,1] = lamb[q1];
Lambda.Alpha[3,2] = alp[q2];
Lambda[3] = lamb_lasso[inds[[3]]];

#####################################
# stop parallel cluster
#####################################
parallel::stopCluster(my.cluster)

#####################################
# write output
#####################################
output = as.data.frame(cbind(Lambda.Alpha, Lambda))
write.csv(output, file = "hyperparams.csv", row.names = FALSE)
