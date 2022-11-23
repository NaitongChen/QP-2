hyperparams = read.csv('hyperparams.csv');
Lambda.Alpha = hyperparams[,1:2]
Lambda = hyperparams[,3]

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

X = read.csv('x_ho.csv');
Y = read.csv('y_ho.csv');

for (k in 1:K) {
  Betar = matrix(0,p,CASE);
  Beta_lasso = matrix(0,p,CASE);
  for (j in 1:CASE) {
    print(c(k,j))
    Betar[,j] = L1PenHuber(Y[(n*(k-1)+1):(n*k),j], X[(n*(k-1)+1):(n*k),], Lambda.Alpha[j,1], Lambda.Alpha[j,2]);
    Betar[,j] = Betar[,j] * (abs(Betar[,j]) > 1e-04);
    
    a = glmnet(x = X[(n*(k-1)+1):(n*k),], y = Y[(n*(k-1)+1):(n*k), j],
               family = "gaussian", alpha = 1, lambda = c(Lambda[j], 0), intercept = FALSE);
    Beta_lasso[,j] = a$beta[,1]
    Beta_lasso[,j] = Beta_lasso[,j] * (abs(Beta_lasso[,j]) > 1e-04);
    
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
    
    
    df = as.data.frame(cbind(Y[(n*(k-1)+1):(n*k),j], X[(n*(k-1)+1):(n*k),1:20]))
    colnames(df)[1] = "Y"
    fit = lm(formula = Y ~ . - 1, data = df)
    beta_oracle = fit$coefficients
    RA_1[k,j] = (norm(Beta_lasso[,j] - beta_0, "2") - norm(beta_oracle - beta_0, "2")) / (norm(Betar[,j] - beta_0, "2") - norm(beta_oracle - beta_0, "2"))
    RA_2[k,j] = (norm(Beta_lasso[,j] - beta_0, "1") - norm(beta_oracle - beta_0, "1")) / (norm(Betar[,j] - beta_0, "1") - norm(beta_oracle - beta_0, "1"))
  }
}

mean(l1loss)
quantile(l1loss, probs = c(0.05, 0.95))
mean(l2loss)
quantile(l2loss, probs = c(0.05, 0.95))
mean(FP)
quantile(FP, probs = c(0.05, 0.95))
mean(FN)
quantile(FN, probs = c(0.05, 0.95))

mean(l1loss_lasso)
quantile(l1loss_lasso, probs = c(0.05, 0.95))
mean(l2loss_lasso)
quantile(l2loss_lasso, probs = c(0.05, 0.95))
mean(FP_lasso)
quantile(FP_lasso, probs = c(0.05, 0.95))
mean(FN_lasso)
quantile(FN_lasso, probs = c(0.05, 0.95))

mean(RA_1)
quantile(RA_1, probs = c(0.05, 0.95))
mean(RA_2)
quantile(RA_2, probs = c(0.05, 0.95))