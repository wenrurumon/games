
#####################################################
# Macro
#####################################################

#Robust Recovery
scale0 <- function(x,to01=T,p=F){
  if(p){x <- pnorm(scale(x))}
  if(to01){
    x <- x-min(x,na.rm=TRUE)
    x <- x/max(x,na.rm=TRUE)
  }
  x[is.na(x)] <- 0.5
  as.numeric(x)
}
scale0_mat <- function(x,to01=T,p=F){
  x <- cbind(x)
  sapply(1:ncol(x),function(i){
    scale0(x[,i],to01,p)
  })
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
  rlt <- list(Z = A%*%Z, X=A, Y=Z ,prop=cumsum(svd(X-E)$d/sum(svd(X-E)$d)))
  return(rlt)  
}
#QPCA
qpca <- function(A,lambda=0,ifscale=TRUE){
  if(ifscale){
    A <- scale(as.matrix(A))
    A[is.na(A)] <- 0
  }else{
    A <- as.matrix(A)
  }
  A.svd <- svd(A)
  d <- A.svd$d-lambda*A.svd$d[1]
  d <- d[d > 1e-8]
  r <- length(d)
  prop <- d^2; info <- sum(prop)/sum(A.svd$d^2);prop <- cumsum(prop/sum(prop))
  d <- diag(d,length(d),length(d))
  u <- A.svd$u[,1:r,drop=F]
  v <- A.svd$v[,1:r,drop=F]
  x <- u%*%sqrt(d)
  y <- sqrt(d)%*%t(v)
  z <- x %*% y
  list(rank=r,X=x,Y=y,Z=x%*%y,prop=prop,info=info)
}
#PCA
pca <- function(X){
  X <- scale(as.matrix(X))
  m = nrow(X)
  n = ncol(X)
  X = scale(X)
  Xeigen <- svd(as.matrix(X))
  value <- (Xeigen$d)^2/m
  value <- cumsum(value/sum(value))
  score <- X %*% Xeigen$v
  mat <- Xeigen$v
  list(score=score,prop=value,mat=mat)
}
#Functional Variable Remodeling
library(fda)
library(flare)
functional_transform <- function(X){
  #calculate the distance for each variables
  X.dist <- dist(t(X))
  pos <- hclust(X.dist)$order
  X <- X[,pos]
  X.dist <- as.matrix(X.dist)[pos,pos]
  pos <- c(0,cumsum(diag(X.dist[-1,-length(pos)])))
  pos <- scale0(pos,T,T)
  #set fourier basis
  fbasis<-create.fourier.basis(c(0,1),nbasis=length(pos)*2-1)
  fphi <- eval.basis(pos,fbasis)
  fcoef <- ginv(t(fphi)%*%fphi)%*%t(fphi)%*%t(X)
  rlt <- t(fcoef-rowMeans(fcoef))/sqrt(nrow(X))
  return(rlt)
}

#####################################################
# Data Process
#####################################################

#Load data
setwd('C:\\Users\\zhu2\\Documents\\kaggle\\hprice')
raw <- read.csv('full.csv',row.names=1)
colnames(raw)
sel <- raw$SalePrice!='tbd'
train <- raw[sel,]
y <- train$SalePrice <- as.numeric(paste(train$SalePrice ))

#Preview Data
test <- function(i){
  summary(lm(y~train[,i]))$r.square
}
test <- sapply(2:ncol(train)-1,function(i){test(i)})
hist(test)

#Matrix Completion Process
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

#Generate Data
x_sep <- do.call(cbind,x_sep)
x_sep <- sapply(1:ncol(x_sep),function(i){scale0(as.numeric(paste(x_sep[,i])),T,F)})
x_mc <- mc0(x_sep,ifprint=TRUE)
# x_f <- functional_transform(x_mc$Z)
x_f <- x_mc$Z
x.qpca <- lapply((0:19)/20,function(l){qpca(x_f,l)})

#####################################################
# Modeling
#####################################################

#Train
models <- lapply(x.qpca,function(x){lm(y~x$X[sel,])})
(check_models <- sapply(models,function(x){summary(x)$r.square}))

x.all <- scale0_mat(x.qpca[[1]]$X)+0.0001
summary(lm(log(y)~x.all[sel,]))

#
x.temp <- x.all
for(i in 1:ncol(x.all)){
  #ori
    temp <- x.temp
    r.ori <- summary(lm(log(y)~x.temp[sel,]))$r.square
  #log
    temp <- x.temp
    temp[,i] <- log(temp[,i])
    r.log <- summary(lm(log(y)~temp[sel,]))$r.square
  #power
    temp <- x.temp
    temp[,i] <- temp[,i]^2
    r.2 <- summary(lm(log(y)~temp[sel,]))$r.square
  #sel
    r <- c(r.ori,r.log,r.2)
    r <- which(r==max(r))
    print(paste(r,i))
    if(r==2){
      x.temp[,i] <- log(x.temp[,i])
    } else if(r==3){
      x.temp[,i] <- x.temp[,i]^2
    }
}
summary(x.lm <- lm(log(y)~x.temp[sel,]))
plot.ts(y);lines(exp(predict(x.lm)),col=2)
qqplot(y,exp(predict(x.lm)))

#Predict
xtest <- cbind(1,x.temp[!sel,])
xcoef <- cbind(coef(x.lm)[1:ncol(xtest)])
ytest <- exp(xtest %*% xcoef)
rlt <- data.frame(Id=rownames(raw)[!sel],SalePrice=ytest)

#Output
write.csv(rlt,file='rlt20160915_1.csv',row.names=F)


