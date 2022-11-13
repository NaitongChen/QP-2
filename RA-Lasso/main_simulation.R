suppressWarnings(library(CVXR, warn.conflicts=FALSE))

n = 100;  # number of observations
p = 400;  # number of variables
CASE = 5;  # 5 scenarios for error distribution
LAMAX = 1.5;  # lambda_max
AMAX = 1.0;  # alpha_max
LENGTH_la = 30;  # number of lambdas used in validation
LENGTH_a = 15;  # number of alphas used in validation
N = 100;  # number of validation datasets
beta_0 = as.matrix(c(3*rep(1,20), rep(0,380)))  # true beta

XV = read.csv('x_ho_v.csv');
YV = read.csv('y_ho_v.csv');
lamb = seq(0.0001,LAMAX,length.out = LENGTH_la);
alp = seq(0.0001,AMAX,length.out = LENGTH_a);

Lambda.Alpha = matrix(0, CASE, 2);  # stores optimal pair of (lambda, alpha) for 6 scenarios

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

# 1st scenario: N(0,4)
loss1 = matrix(0, LENGTH_la, LENGTH_a);
for (k in 1:N) {
  for (i in 1:LENGTH_la) {
    for (j in 1:LENGTH_a) {
      print(c(k,i,j))
      betah = L1PenHuber(YV[(n*(k-1)+1):(n*k),1], XV[(n*(k-1)+1):(n*k),], lamb[i], alp[j]);
      betah = betah * (abs(betah) > 1e-04);
      loss1[i,j] = loss1[i,j] + norm(betah-beta_0, "F");
    }
  }
}

inds = which(loss1 == min(loss1), arr.ind = TRUE);
q1 = inds[1]
q2 = inds[2]
minvalue = loss1[q1,q2]
Lambda.Alpha[1,1] = lamb[q1];
Lambda.Alpha[1,2] = alp[q2]