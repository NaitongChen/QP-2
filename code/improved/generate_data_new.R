# setwd("~/simulation")
N = 100
n = 100
contaminated = 10
p = 400
CASE = 5
beta_0 = c(rep(3,20), rep(0,380))

#########################################################
######            Homoscedastic Case               ######
#########################################################

generate_homo <- function(n,p,beta_0){  ## Generate data for the homoscedastic case
  X = matrix(rnorm(n*p), n, p)
  X_lev1 = X
  X_lev2 = X
  Y = error = matrix(0, n, CASE)
  #(1) N(0,4)
  error[,1] = rnorm(n, mean = 0, sd = 2)
  #(4) t_2
  error[,2] = rt(n, df = 2)
  #(3) MixN
  index = sample(c(1,2), size = n, replace = TRUE, prob = c(0.5, 0.5))
  error[which(index==1),3] = rnorm(length(which(index==1)), mean = -1, sd = 2)
  error[which(index==2),3] = rnorm(length(which(index==2)), mean = 8, sd = 1)
  error[,3] = (error[,3] - 3.5)
  #(5) High leverage 1
  for(j in 1:contaminated){
    a_t = sample(c(-1,1), p, replace = TRUE)
    a = a_t - (1/p) * sum(a_t)
    X_lev1[j,] = rnorm(p, sd = 0.1) + (1/norm(a, "2")) * a
  }
  error[,4] = rnorm(n, mean = 0, sd = 1)
  #(6) High leverage 2
  for(j in 1:contaminated){
    a_t = sample(c(-1,1), p, replace = TRUE)
    a = a_t - (1/p) * sum(a_t)
    X_lev1[j,] = rnorm(p, sd = 0.1) + (2/norm(a, "2")) * a
  }
  error[,5] = rnorm(n, mean = 0, sd = 1)
  for(j in 1:ncol(Y)){
    if (j != 4 & j != 5) {
      Y[,j] = X %*% beta_0 + error[,j] 
    } else if (j == 4) {
      Y[,j] = X_lev1 %*% beta_0 + error[,j] 
    } else {
      Y[,j] = X_lev2 %*% beta_0 + error[,j] 
    }
  }
  return(list(X = X, X_lev1 = X_lev1, X_lev2 = X_lev2, Y = Y))
}

totalX = matrix(0, n*N, p)
totalX_lev1 = matrix(0, n*N, p)
totalX_lev2 = matrix(0, n*N, p)
totalY = matrix(0, n*N, CASE)
for (k in 1:N){
  cat("k =", k, "\n")
  data = generate_homo(n,p,beta_0)
  totalX[(n*(k-1)+1):(n*k),] = data$X
  totalX_lev1[(n*(k-1)+1):(n*k),] = data$X_lev1
  totalX_lev2[(n*(k-1)+1):(n*k),] = data$X_lev2
  totalY[(n*(k-1)+1):(n*k),] = data$Y
} 

write.csv(totalX, file = "x_cv.csv", row.names = FALSE)
write.csv(totalX_lev1, file = "x_lev1_cv.csv", row.names = FALSE)
write.csv(totalX_lev2, file = "x_lev2_cv.csv", row.names = FALSE)
write.csv(totalY, file = "y_cv.csv", row.names = FALSE)