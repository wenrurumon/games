
#####################################################
# Macro
#####################################################

#Robust Recovery

scale0 <- function(x){
  x <- pnorm(scale(x))
  x[is.na(x)] <- 0.5
  x
}

mc0 <- function(X,lambda=0.2,ifprint=FALSE){
  ### input the initial values
  m <- dim(X)[1]
  n <- dim(X)[2]                                  
  Z <- matrix(rep(0,n*n),nrow=n)
  J <- matrix(rep(0,n*n),nrow=n)
  E <- matrix(rep(0,m*n),nrow=m)
  Y1 <- matrix(rep(0,m*n),nrow=m)
  Y2 <- matrix(rep(0,n*n),nrow=n)
  mu <- 10^(-6) 
  mu_max <- 10^6
  rho <- 1.1
  epsilon <- 10^(-8)
  A <- X
  I <- diag(n)
  ### compute the while loop to run matrix completion
  while ((max(abs(X-A %*% Z -E)) > epsilon) | (max(abs(Z-J)) > epsilon)) {
    if(ifprint){print(paste((max(abs(Z-J))),(max(abs(X-A %*% Z -E)))))}
    ## update the J matrix
    U_J_svd  <-  svd(Z+Y2/mu)
    big_eps <- diag(U_J_svd$d - 1/mu)
    big_eps[big_eps<0] <- 0
    J <- U_J_svd$u %*% big_eps %*% t(U_J_svd$v)
    ## update the Z matrix
    I_A_svd <- svd(I + t(A) %*% A )
    Z <- I_A_svd$v %*% diag(1/(I_A_svd$d)) %*% t(I_A_svd$u) %*% (t(A)%*%(X-E)+J
                                                                 +(t(A) %*% Y1 - Y2)/mu)
    ## update the E matrix
    Q <- X - A %*% Z + Y1/mu
    Q_col_norm <- sqrt(colSums(Q^2))
    Q_mlambda_o_mu <- Q_col_norm - lambda/mu
    Q_mlambda_o_mu[Q_mlambda_o_mu < 0] <- 0 
    Q_c_coef <- matrix(rep((Q_mlambda_o_mu/Q_col_norm),m),nrow=m,byrow=TRUE)
    E <- Q_c_coef * Q
    ## update the multipliers
    Y1 <- Y1 + mu*(X-A%*%Z-E)
    Y2 <- Y2 + mu*(Z-J)
    ## update the parameter mu
    mu <- min(rho*mu,mu_max)
  }
  #rlt <- A
  rlt <- list(score=A,prop=cumsum(svd(X-E)$d/sum(svd(X-E)$d)))
  return(rlt)  
}

#####################################################
# Data Process
#####################################################

#Load data

setwd('C:\\Users\\zhu2\\Documents\\kaggle\\hprice')
raw <- read.csv('full.csv',row.names=1)
colnames(raw)
train <- raw[raw$SalePrice!='tbd',]
y <- train$SalePrice <- as.numeric(paste(train$SalePrice ))

#Preview Data

test <- function(i){
  summary(lm(y~train[,i]))$r.square
}
test <- sapply(2:ncol(train)-1,function(i){test(i)})
hist(test)

trans <- function(x){
  out <- outer(x,unique(x),'==')+0
  colnames(out) <- unique(x)
  out
}
x_sep <- lapply(2:ncol(raw)-1,function(i){
  if(is.numeric(raw[,i])){
    return((raw[,i,drop=F]))
  } else {
    return(trans(paste(raw[,i])))
  }
})

x_full <- scale0(do.call(cbind,x_sep))
x_mc <- mc(x_full,ifprint=TRUE)

