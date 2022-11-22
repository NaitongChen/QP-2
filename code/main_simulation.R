suppressWarnings(library(CVXR, warn.conflicts=FALSE))
library(glmnet)
library(doParallel)
library(foreach)

# set up parallel computing
n.cores <- parallel::detectCores() - 1;

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
);

doParallel::registerDoParallel(cl = my.cluster);

foreach::getDoParRegistered()
foreach::getDoParWorkers()


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

XV = read.csv('x_ho_v.csv');
YV = read.csv('y_ho_v.csv');
lamb = seq(0.0001,LAMAX,length.out = LENGTH_la);
lamb_lasso = seq(0.00001,LAMAX_lasso,length.out = LENGTH_la_lasso);
alp = seq(0.0001,AMAX,length.out = LENGTH_a);

Lambda.Alpha = matrix(0, CASE, 2);  # stores optimal pair of (lambda, alpha) for 4 scenarios
Lambda = matrix(0, CASE);

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

TuneHyperParam <- function(case) {
  loss1 = matrix(0, LENGTH_la, LENGTH_a);
  for (k in 1:N) {
    for (i in 1:LENGTH_la) {
      print(c(k,i))
      loss1_i = foreach (j = 1:LENGTH_a, .combine = 'c', .packages = "CVXR") %dopar% {
        betah = L1PenHuber(YV[(n*(k-1)+1):(n*k), case], XV[(n*(k-1)+1):(n*k),], lamb[i], alp[j]);
        betah = betah * (abs(betah) > 1e-04);
        betah = norm(betah-beta_0, "F");
      }
      loss1[i,] = loss1[i,] + loss1_i;
    }
  }
  
  loss1_lasso = matrix(0, LENGTH_la_lasso);
  for (k in 1:N) {
    a <- glmnet(x = XV[(n*(k-1)+1):(n*k),], y = YV[(n*(k-1)+1):(n*k), case],
                family = "gaussian", alpha = 1, lambda = lamb_lasso, intercept = FALSE)
    
    for (i in 1:LENGTH_la_lasso) {
      print(c(k,i))
      betah = a$beta[,i]
      betah = betah * (abs(betah) > 1e-04);
      loss1_lasso[i] = loss1_lasso[i] + norm(betah-beta_0, "F");
    }
  }
  
  inds = which(loss1 == min(loss1), arr.ind = TRUE);
  q1 = inds[1]
  q2 = inds[2]
  
  ind = which(loss1_lasso == min(loss1_lasso))
  
  return(c(q1, q2, ind))
}

# 1st scenario: N(0,4)
inds = TuneHyperParam(1);
q1 = inds[1];
q2 = inds[2];
Lambda.Alpha[1,1] = lamb[q1];
Lambda.Alpha[1,2] = alp[q2];
Lambda[1] = lamb_lasso[inds[3]]

# 2nd scenario: 2t_3
inds = TuneHyperParam(2);
q1 = inds[1];
q2 = inds[2];
Lambda.Alpha[2,1] = lamb[q1];
Lambda.Alpha[2,2] = alp[q2];
Lambda[2] = lamb_lasso[inds[3]]

# 3rd scenario: MixN
inds = TuneHyperParam(3);
q1 = inds[1];
q2 = inds[2];
Lambda.Alpha[3,1] = lamb[q1];
Lambda.Alpha[3,2] = alp[q2];
Lambda[3] = lamb_lasso[inds[3]]