# setwd("~/simulation")
N = 100
n = 100
p = 400
CASE = 3
beta_0 = c(rep(3,20), rep(0,380))

#########################################################
######            Homoscedastic Case               ######
#########################################################

generate_homo <- function(n,p,beta_0){  ## Generate data for the homoscedastic case
  X = matrix(rnorm(n*p), n, p)
  Y = error = matrix(0, n, CASE)
  #(1) N(0,4)
  error[,1] = rnorm(n, mean = 0, sd = 2)
  error[,1] = (error[,1] - mean(error[,1]))
  #(2) 2t_3
  error[,2] = sqrt(2) * rt(n, df = 3)
  error[,2] = (error[,2] - mean(error[,2]))
  #(3) MixN
  index = sample(c(1,2), size = n, replace = TRUE, prob = c(0.5, 0.5))
  error[which(index==1),3] = -abs(rnorm(length(which(index==1)), mean = -1, sd = 2))
  error[which(index==2),3] = abs(rnorm(length(which(index==2)), mean = 8, sd = 1))
  error[,3] = (error[,3] - mean(error[,3]))
  for(j in 1:ncol(Y)){
    Y[,j] = X %*% beta_0 + error[,j]
  }
  return(list(X = X, Y = Y))
}

totalX_v = totalX = matrix(0, n*N, p)
totalY_v = totalY = matrix(0, n*N, CASE)
for (k in 1:N){
  cat("k =", k, "\n")
  data_v = generate_homo(n,p,beta_0)
  totalX_v[(n*(k-1)+1):(n*k),] = data_v$X
  totalY_v[(n*(k-1)+1):(n*k),] = data_v$Y
  data = generate_homo(n,p,beta_0)
  totalX[(n*(k-1)+1):(n*k),] = data$X
  totalY[(n*(k-1)+1):(n*k),] = data$Y
} 

write.csv(totalX_v, file = "x_ho_v.csv", row.names = FALSE)
write.csv(totalY_v, file = "y_ho_v.csv", row.names = FALSE)
write.csv(totalX, file = "x_ho.csv", row.names = FALSE)
write.csv(totalY, file = "y_ho.csv", row.names = FALSE)