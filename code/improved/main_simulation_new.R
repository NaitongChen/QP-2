suppressWarnings(library(CVXR, warn.conflicts=FALSE))
library(glmnet)
library(doParallel)
library(foreach)
library(pense)
set.seed(2022)

n = 100;  # number of observations
p = 400;  # number of variables
CASE = 5;  # 5 scenarios for error distribution
LAMAX = 1.5;  # lambda_max
LAMAX_lasso = 10;
AMAX = 1.0;  # alpha_max
LENGTH_la = 30;  # number of lambdas used in validation
LENGTH_la_lasso = 200;
LENGTH_a = 15;  # number of alphas used in validation
N = 20;  # number of validation datasets
beta_0 = as.matrix(c(3*rep(1,20), rep(0,380)))  # true beta

X_cv = read.csv('x_cv.csv');
X_lev1_cv = read.csv('x_lev1_cv.csv');
X_lev2_cv = read.csv('x_lev2_cv.csv');
Y_cv = read.csv('y_cv.csv');

lamb = seq(0.0001,LAMAX,length.out = LENGTH_la);
lamb_lasso = rev(seq(0.00001,LAMAX_lasso,length.out = LENGTH_la_lasso));
lamb_pense = rev(seq(0.00001,LAMAX_lasso,length.out = LENGTH_la_lasso));
alp = seq(0.0001,AMAX,length.out = LENGTH_a);

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

# set up parallel computing
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
                                  "X_cv", "X_lev1_cv", "X_lev2_cv", "Y_cv", "lamb", "lamb_lasso", "lamb_pense", "alp"),
              envir=environment())

TuneHyperParam <- function(case, nfold, rep) {
  if (case <= 3) {
    X_curr = X_cv[(n*(rep-1)+1):(n*rep),]
    Y_curr = Y_cv[(n*(rep-1)+1):(n*rep), case]
  } else if (case == 4) {
    X_curr = X_lev1_cv[(n*(rep-1)+1):(n*rep),]
    Y_curr = Y_cv[(n*(rep-1)+1):(n*rep), case]
  } else {
    X_curr = X_lev2_cv[(n*(rep-1)+1):(n*rep),]
    Y_curr = Y_cv[(n*(rep-1)+1):(n*rep), case]
  }
  
  # standardize data
  Y_curr = (Y_curr - mean(Y_curr)) / sd(Y_curr)
  for (j in 1:p) {
    X_curr[,j] = (X_curr[,j] - mean(X_curr[,j])) / sd(X_curr[,j])
  }
  
  # assign k folds
  ii <- (1:n) %% nfold + 1
  ii <- sample(ii)
  
  print("RA lasso")
  loss1 = matrix(0, LENGTH_la, LENGTH_a);
  for (k in 1:nfold) {
    print(k)
    loss1_ij = foreach (i = 1:LENGTH_la, .combine = 'cbind', .packages = c("CVXR", "Matrix"), .export = "L1PenHuber") %:%
      foreach (j = 1:LENGTH_a, .combine = 'c', .packages = c("CVXR", "Matrix"), .export ="L1PenHuber") %dopar% {
        betah = L1PenHuber(Y_curr[ii != k], X_curr[ii != k,], lamb[i], alp[j]);
        betah = norm(matrix(Y_curr[ii == k]) - as.matrix(X_curr[ii == k,]) %*% betah, "1");
      }
    loss1 = loss1 + t(loss1_ij)
  }
  
  print("lasso")
  loss1_lasso = matrix(0, LENGTH_la_lasso);
  for (k in 1:nfold) {
    print(k)
    a <- glmnet(x = as.matrix(X_curr[ii != k,]), y = matrix(Y_curr[ii != k]),
                family = "gaussian", alpha = 1, lambda = lamb_lasso, intercept = FALSE)
    
    for (i in 1:LENGTH_la_lasso) {
      betah = matrix(a$beta[,i])
      loss1_lasso[i] = loss1_lasso[i] + norm(matrix(Y_curr[ii == k]) - as.matrix(X_curr[ii == k,]) %*% betah, "1");
    }
  }
  
  print("pense")
  loss1_pense = matrix(0, LENGTH_la_lasso);
  for (k in 1:nfold) {
    print(k)
    Y_curr_k = matrix(Y_curr[ii == k])
    X_curr_k = as.matrix(X_curr[ii == k,])
    loss1_pense_k = foreach (i = 1:LENGTH_la_lasso, .combine = 'c', .packages = "pense") %dopar% {
      betah = pense(as.matrix(X_curr[ii != k,]), Y_curr[ii != k], alpha = 1, lambda = lamb_pense[i], intercept = FALSE)$estimates[[1]]$beta 
      betah = norm(Y_curr_k - X_curr_k %*% betah, "1");
    }
    
    loss1_pense = loss1_pense + loss1_pense_k
  }
  
  inds = which(loss1 == min(loss1), arr.ind = TRUE);
  q1 = inds[1]
  q2 = inds[2]
  
  ind_lasso = which(loss1_lasso == min(loss1_lasso))
  
  ind_pense = which(loss1_pense == min(loss1_pense))
  
  return(list(q1, q2, ind_lasso, ind_pense, loss1, loss1_lasso, loss1_pense))
}

# compute metrics
K=20;
CASE=5;

l2loss = matrix(0,K,CASE);
l1loss = matrix(0,K,CASE);
FP = matrix(0,K,CASE);
FN = matrix(0,K,CASE);
RA_1 = matrix(0,K,CASE);
RA_2 = matrix(0,K,CASE);

l2loss_pense = matrix(0,K,CASE);
l1loss_pense = matrix(0,K,CASE);
FP_pense = matrix(0,K,CASE);
FN_pense = matrix(0,K,CASE);
RA_1_pense = matrix(0,K,CASE);
RA_2_pense = matrix(0,K,CASE);

l2loss_lasso = matrix(0,K,CASE);
l1loss_lasso = matrix(0,K,CASE);
FP_lasso = matrix(0,K,CASE);
FN_lasso = matrix(0,K,CASE);

beta_0_copy = beta_0

for (k in 1:K) {
  for (j in 1:CASE) {
    print(c(k,j))
    if (j <= 3) {
      X_curr = X_cv[(n*(k-1)+1):(n*k),]
      Y_curr = Y_cv[(n*(k-1)+1):(n*k), j]
    } else if (j == 4) {
      X_curr = X_lev1_cv[(n*(k-1)+1):(n*k),]
      Y_curr = Y_cv[(n*(k-1)+1):(n*k), j]
    } else {
      X_curr = X_lev2_cv[(n*(k-1)+1):(n*k),]
      Y_curr = Y_cv[(n*(k-1)+1):(n*k), j]
    }
    
    df = as.data.frame(cbind(Y_curr, X_curr[,1:20]))
    colnames(df)[1] = "Y"
    fit = lm(formula = Y ~ . - 1, data = df)
    
    beta_oracle = rep(0, p+1)
    beta_oracle[2:21] = fit$coefficients
    beta_0 = c(0, beta_0_copy)
    
    # standardize data
    mean_sd = matrix(0, 2, p+1)
    
    mean_sd[1,1] = mean(Y_curr)
    mean_sd[2,1] = sd(Y_curr)
    Y_curr = (Y_curr - mean_sd[1,1]) / mean_sd[2,1]
    Y_curr = matrix(Y_curr)
    
    for (jj in 1:p) {
      mean_sd[1,jj+1] = mean(X_curr[,jj])
      mean_sd[2,jj+1] = sd(X_curr[,jj])
      X_curr[,jj] = (X_curr[,jj] - mean_sd[1,jj+1]) / mean_sd[2,jj+1]
    }
    X_curr = as.matrix(X_curr)
    
    inds = TuneHyperParam(j, 5, k)
    
    betah_ra = L1PenHuber(Y_curr, X_curr, lamb[inds[[1]]], alp[inds[[2]]]);
    betah_ra = betah_ra * (abs(betah_ra) > 1e-04);
    
    a = glmnet(x = X_curr, y = Y_curr,
               family = "gaussian", alpha = 1, lambda = c(lamb_lasso[inds[[3]]], 0), intercept = FALSE);
    betah_lasso = a$beta[,1]
    betah_lasso = betah_lasso * (abs(betah_lasso) > 1e-04);
    
    betah_pense = pense(X_curr, as.vector(Y_curr), alpha = 1, lambda = lamb_pense[inds[[4]]], intercept = FALSE)$estimates[[1]]$beta 
    betah_pense = betah_pense * (abs(betah_pense) > 1e-04);
    
    # de-standardize coefficients
    betah_ra = c(0, betah_ra)
    betah_lasso = c(0, betah_lasso)
    betah_pense = c(0, betah_pense)
    
    betah_ra[1] = mean_sd[1,1] - mean_sd[2,1] * sum(betah_ra[2:(p+1)] * (mean_sd[1,2:(p+1)] / mean_sd[2,2:(p+1)]))
    betah_lasso[1] = mean_sd[1,1] - mean_sd[2,1] * sum(betah_lasso[2:(p+1)] * (mean_sd[1,2:(p+1)] / mean_sd[2,2:(p+1)]))
    betah_pense[1] = mean_sd[1,1] - mean_sd[2,1] * sum(betah_pense[2:(p+1)] * (mean_sd[1,2:(p+1)] / mean_sd[2,2:(p+1)]))
    
    for (jj in 2:(p+1)) {
      betah_ra[jj] = mean_sd[2,1] / mean_sd[2,jj] * betah_ra[jj]
      betah_lasso[jj] = mean_sd[2,1] / mean_sd[2,jj] * betah_lasso[jj]
      betah_pense[jj] = mean_sd[2,1] / mean_sd[2,jj] * betah_pense[jj]
    }

    # compute metric
    l2loss[k,j] = norm(matrix(betah_ra - beta_0), "2");
    l1loss[k,j] = norm(matrix(betah_ra - beta_0), "1");
    temp = betah_ra;
    FP[k,j] = sum(temp[which(beta_0 == 0)] != 0);
    FN[k,j] = sum(temp[which(beta_0 != 0)] == 0);
    RA_1[k,j] = (norm(matrix(betah_lasso - beta_0), "2") - norm(matrix(beta_oracle) - beta_0, "2")) / (norm(matrix(betah_ra - beta_0), "2") - norm(matrix(beta_oracle - beta_0), "2"))
    RA_2[k,j] = (norm(matrix(betah_lasso - beta_0), "1") - norm(matrix(beta_oracle) - beta_0, "1")) / (norm(matrix(betah_ra - beta_0), "1") - norm(matrix(beta_oracle - beta_0), "1"))
    
    l2loss_pense[k,j] = norm(matrix(betah_pense - beta_0), "2");
    l1loss_pense[k,j] = norm(matrix(betah_pense - beta_0), "1");
    temp = betah_pense;
    FP_pense[k,j] = sum(temp[which(beta_0 == 0)] != 0);
    FN_pense[k,j] = sum(temp[which(beta_0 != 0)] == 0);
    RA_1_pense[k,j] = (norm(matrix(betah_lasso - beta_0), "2") - norm(matrix(beta_oracle - beta_0), "2")) / (norm(matrix(betah_pense - beta_0), "2") - norm(matrix(beta_oracle - beta_0), "2"))
    RA_2_pense[k,j] = (norm(matrix(betah_lasso - beta_0), "1") - norm(matrix(beta_oracle - beta_0), "1")) / (norm(matrix(betah_pense - beta_0), "1") - norm(matrix(beta_oracle - beta_0), "1"))
    
    l2loss_lasso[k,j] = norm(matrix(betah_lasso - beta_0), "2");
    l1loss_lasso[k,j] = norm(matrix(betah_lasso - beta_0), "1");
    temp = betah_lasso;
    FP[k,j] = sum(temp[which(beta_0 == 0)] != 0);
    FN[k,j] = sum(temp[which(beta_0 != 0)] == 0);
  }
}

parallel::stopCluster(my.cluster)

write.csv(l2loss, file = "l2loss.csv", row.names = FALSE)
write.csv(l1loss, file = "l1loss.csv", row.names = FALSE)
write.csv(FP, file = "FP.csv", row.names = FALSE)
write.csv(FN, file = "FN.csv", row.names = FALSE)
write.csv(RA_1, file = "RA_1.csv", row.names = FALSE)
write.csv(RA_2, file = "RA_2.csv", row.names = FALSE)

write.csv(l2loss_pense, file = "l2loss_pense.csv", row.names = FALSE)
write.csv(l1loss_pense, file = "l1loss_pense.csv", row.names = FALSE)
write.csv(FP_pense, file = "FP_pense.csv", row.names = FALSE)
write.csv(FN_pense, file = "FN_pense.csv", row.names = FALSE)
write.csv(RA_1_pense, file = "RA_1_pense.csv", row.names = FALSE)
write.csv(RA_2_pense, file = "RA_2_pense.csv", row.names = FALSE)

write.csv(l2loss_lasso, file = "l2loss_lasso.csv", row.names = FALSE)
write.csv(l1loss_lasso, file = "l1loss_lasso.csv", row.names = FALSE)
write.csv(FP_lasso, file = "FP_lasso.csv", row.names = FALSE)
write.csv(FN_lasso, file = "FN_lasso.csv", row.names = FALSE)
