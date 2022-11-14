library(pense)

n = 100;  # number of observations
p = 400;  # number of variables
LAMAX = 1.5;  # lambda_max
LENGTH_la = 30;  # number of lambdas used in validation

XV = read.csv('x_ho_v.csv');
YV = read.csv('y_ho_v.csv');

lamb = seq(0.0001,LAMAX,length.out = LENGTH_la);

k=1
i=1
j=1

Y = YV[(n*(k-1)+1):(n*k),1]
X = as.matrix(XV[(n*(k-1)+1):(n*k),])
lambda = lamb[i]

result = pense(X, Y, alpha=1, lambda=lambda, intercept=FALSE)
betah = result$estimates[[1]]$beta
