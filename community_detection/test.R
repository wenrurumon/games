
###############################################
# Test 1
###############################################

#Dummy data
x <- (make_full_graph(7)) + (make_full_graph(7) - path(1:5)) + (make_ring(15))
x <- as.matrix(get.adjacency(x)); plotnet(x)
set.seed(10);x[sample(which(x==0),nrow(x)/3)] <- runif(nrow(x)/3,0,0.5) ;diag(x) <- 0
dimnames(x) <- list(1:ncol(x),1:ncol(x))
x <- x + t(x); plotnet(x)

par(mfrow=c(2,4));igraph.options(vertex.size=10)
centers <- 3
plotclust(x,cutnet(x,0.51,T)$cluster,'Actual')
plotclust(x,fcclust(x)$cluster,'GREEDY')
plotclust(x,fcclust(x,TRUE)$cluster,'GREEDY_W')
plotclust(x,fcclusts(x,l=0,thres=0,w=TRUE)$cluster,'GREEDY_Ws')
plotclust(x,spclust(x,centers)$cluster,'SPCLUST')
plotclust(x,kmclust(x,centers)$cluster,'KMEANS')
plotclust(x,dkmclust(x,centers)$cluster,'DKMEANS')
plotclust(x,sparsehclust(x,centers)$cluster,'SparseHC')

par(mfrow=c(2,3))
plotclust(x,cutnet(x,0.51,T)$cluster,'Actual')
plotclust(x,fcclusts(x,l=1,thres=0,w=TRUE)$cluster,'GREEDY_Ws1')
plotclust(x,fcclusts(x,l=2,thres=0,w=TRUE)$cluster,'GREEDY_Ws2')
plotclust(x,fcclusts(x,l=3,thres=0,w=TRUE)$cluster,'GREEDY_Ws3')
plotclust(x,fcclusts(x,l=4,thres=0,w=TRUE)$cluster,'GREEDY_Ws4')
plotclust(x,fcclusts(x,l=0,thres=0,w=TRUE)$cluster,'GREEDY_Ws0')

par(mfrow=c(2,3))
plotclust(x,cutnet(x,0.51,T)$cluster,'Actual')
plotclust(x,fcclusts(x,l=1,thres=0,w=F)$cluster,'GREEDY_Ws1')
plotclust(x,fcclusts(x,l=2,thres=0,w=F)$cluster,'GREEDY_Ws2')
plotclust(x,fcclusts(x,l=3,thres=0,w=F)$cluster,'GREEDY_Ws3')
plotclust(x,fcclusts(x,l=4,thres=0,w=F)$cluster,'GREEDY_Ws4')
plotclust(x,fcclusts(x,l=0,thres=0,w=F)$cluster,'GREEDY_Ws0')

#Mixclust
par(mfrow=c(2,2))
plotclust(x,cutnet(x,0.51,T)$cluster,'Actual')
plotclust(x,mixclust(x,thres=0,w=T,b.rank=2.5,layer=Inf)$cluster,'mixclust 2.5')
plotclust(x,mixclust(x,thres=0,w=T,b.rank=2,layer=Inf)$cluster,'mixclust 2')
plotclust(x,mixclust(x,thres=0,w=T,b.rank=1.5,layer=Inf)$cluster,'mixclust 1.5')

###############################################
# Test 2
###############################################

#Complex Dummy Data
par(mfrow=c(2,4))
x <- erdos.renyi.game(100, 1/80) #Random Network
x <- as.matrix(get.adjacency(x));dimnames(x) <- list(1:100,1:100)
x[x>0] <- runif(sum(x>0),0.51,1) #Applying Weighting to the existed edges
x[sample(which(x==0),sum(x==0)/100)] <- runif(sum(x==0)/100,0,0.5); diag(x) <- 0 #Applying noice weights

centers <- 2
plotclust(x,cutnet(x,0.51,T)$cluster,'Actual')
plotclust(x,fcclust(x)$cluster,'GREEDY')
plotclust(x,fcclust(x,TRUE)$cluster,'GREEDY_W')
plotclust(x,spclust(x,centers)$cluster,'SPCLUST')
plotclust(x,kmclust(x,centers)$cluster,'KMEANS')
# plotclust(x,dkmclust(x,centers)$cluster,'DKMEANS')
plotclust(x,sparsehclust(x,centers)$cluster,'SparseHC')
plotclust(x,mixclust(x)$cluster,'Mixclust')
plotclust(x,fcclusts(x,l=0,thres=0,w=TRUE)$cluster,'GREEDY_Ws')

# par(mfrow=c(2,2))
# plotclust(x,cutnet(x,0.51,T)$cluster,'Actual')
# plotclust(x,dkmclust(x,centers,layers=1)$cluster,'DKMEANS1')
# plotclust(x,dkmclust(x,centers,layers=2)$cluster,'DKMEANS2')
# plotclust(x,dkmclust(x,centers,layers=3)$cluster,'DKMEANS3')

###############################################
# Test 3
###############################################

#Actual data
par(mfrow=c(2,4))
setwd('C:\\Users\\zhu2\\Documents\\getpathway');load('rlt_v4.rda')
groups=10
raw <- rlt[1:groups]
raw <- lapply(raw,function(x){x[[3]][[1]][[2]][,c(2,4,6),drop=F]})
forcheck <- lapply(raw,function(x){unique(c(x[,1],x[,2]))})
raw <- do.call(rbind,raw)
pathlist <- unique(as.vector(t(raw[,1:2])))
forcheck <- lapply(forcheck,function(x){match(x,pathlist)})
raw[,1] <- match(raw[,1],pathlist);raw[,2] <- match(raw[,2],pathlist)
raw <- cbind(response=as.numeric(raw[,1]),predict=as.numeric(raw[,2]),score=1-as.numeric(raw[,3]))
xbk <- x <- as.matrix(simple_triplet_matrix(c(raw[,1],max(raw[,2])),c(raw[,2],max(raw[,1])),c(raw[,3],0)))
set.seed(1000);x[sample(which(x==0),round(sum(x>0)/5))] <- runif(round(sum(x>0)/5),0,0.5);diag(x) <- 0
dimnames(x) <- list(1:ncol(x),1:nrow(x))

centers <- groups
plotclust(x,cutnet(x,0.51,T)$cluster,'Actual')
plotclust(x,fcclust(x)$cluster,'GREEDY')
plotclust(x,fcclust(x,TRUE)$cluster,'GREEDY_W')
plotclust(x,spclust(x,centers)$cluster,'SPCLUST')
plotclust(x,kmclust(x,centers)$cluster,'KMEANS')
plotclust(x,dkmclust(x,centers)$cluster,'DKMEANS')
plotclust(x,sparsehclust(x,centers)$cluster,'SparseHC')
plotclust(x,mixclust(x)$cluster,'Mix')
# plotclust(x,mixclust(x,b.rank=0)$cluster,'Mix,0')

par(mfrow=c(2,2))
plotclust(x,cutnet(x,0.51,T)$cluster,'Actual')
plotclust(x,mixclust(x,thres=0,w=T,b.rank=3,layer=Inf)$cluster,'mixclust 3')
plotclust(x,mixclust(x,thres=0,w=T,b.rank=2.5,layer=Inf)$cluster,'mixclust 2.5')
plotclust(x,mixclust(x,thres=0,w=T,b.rank=2.5,lambda=1,layer=Inf)$cluster,'mixclust 0')

##############################################
# Quadratic Transformation
##############################################

#Generate Data
par(mfrow=c(2,3))
setwd('C:\\Users\\zhu2\\Documents\\getpathway');load('rlt_v4.rda')
groups=8
raw <- rlt[1:groups]
raw <- lapply(raw,function(x){x[[3]][[1]][[2]][,c(2,4,6),drop=F]})
forcheck <- lapply(raw,function(x){unique(c(x[,1],x[,2]))})
raw <- do.call(rbind,raw)
pathlist <- unique(as.vector(t(raw[,1:2])))
forcheck <- lapply(forcheck,function(x){match(x,pathlist)})
raw[,1] <- match(raw[,1],pathlist);raw[,2] <- match(raw[,2],pathlist)
raw <- cbind(response=as.numeric(raw[,1]),predict=as.numeric(raw[,2]),score=1-as.numeric(raw[,3]))
xbk <- x <- as.matrix(simple_triplet_matrix(c(raw[,1],max(raw[,2])),c(raw[,2],max(raw[,1])),c(raw[,3],0)))
set.seed(10);x[sample(which(x==0),round(sum(x>0)/10))] <- runif(round(sum(x>0)/10),0,1);diag(x) <- 0;plotnet(x)
dimnames(xbk) <- dimnames(x) <- list(1:nrow(x),1:ncol(x))

#Test
z <- qpca(x,0.1)$Z
centers <- groups
plotclust(x,fcclust(xbk)$cluster,'answer')
plotclust(x,fcclust(x,w=FALSE)$cluster,'GREEDY')
plotclust(x,fcclust(x,w=TRUE)$cluster,'GREEDY_w')
plotclust(x,fcclust(z,w=TRUE)$cluster,'GREEDY_wz')
z[x==0] <- 0;plotclust(x,fcclust(z,w=TRUE)$cluster,'GREEDY_wz2')

###############################################
# Challenge
###############################################

setwd('C:\\Users\\zhu2\\Documents\\dreamer\\subchallenge1')
library(data.table)
# library(dplyr)
d <- 1
x <- fread(dir()[d])
v1 <- c(x$V1)+1
v2 <- c(x$V2)+1
v3 <- c(x$V3)
# for(i in 1:length(v1)){
#   temp <- sort(c(v1[i],v2[i]))
#   v1[i] <- temp[1]
#   v2[i] <- temp[2]
# }
v1 <- c(v1,max(v1,v2))
v2 <- c(v2,max(v1,v2))
v3 <- c(v3,0)
x <- slam::simple_triplet_matrix(v1,v2,v3)
x <- as.matrix(x)
x.raw <- x <- x+t(x)

system.time(rlt <- mixclust(x.raw,thres=0,w=T,b.rank=2.5,layer=Inf,lambda=0.5,maxrank=3))
sort(table(rlt$cluster))

#Review Result
x.cutnet <- rlt
x.cutnet[[1]] <- x.cutnet[[1]][order(-sapply(x.cutnet$subnets,ncol))]
x.cutnet[[2]] <- match(x.cutnet[[2]],as.numeric(names(table(x.cutnet[[2]])[order(-table(x.cutnet[[2]]))])))
x.cutnet[[3]] <- sapply(x.cutnet[[1]],mat.degree)

tempf2 <- function(x,membership=NULL,main='',cuts=0){
  G <- graph_from_adjacency_matrix(x>0)
  if(is.null(membership)){membership=rep(1,ncol(x))}
  plot(create.communities(G, membership), 
       as.undirected(G), 
       layout=layout.kamada.kawai(as.undirected(G)),
       main=main,
       edge.arrow.size=1,vertex.size=3,vertex.label.cex=.1,edge.width= .1
       )}
tempf <- function(sel){
  sel <- x.cutnet$cluster%in%sel
  xt <- x[sel,sel]
  ct <- x.cutnet[[2]][sel]
  tempf2(xt,ct)
}
head(-sort(-table(rlt$cluster)),10)
tempf(2:3)

i <- c(1,2)
par(mfrow=c(2,2))
for (l in 1:4){i[-1] <- i[-1]+1;tempf(i)}
par(mfrow=c(3,3))
for(l in 1:9){plotnet(x[x.cutnet[[2]]==l,x.cutnet[[2]]==l],.1,5,.1,.1,0)}
for(l in 1:9){plotnet(x[x.cutnet[[2]]==l,x.cutnet[[2]]==l],.1,5,.1,.1,1)}

out <- sapply(1:length(x.cutnet[[1]]),function(g){paste0(which(x.cutnet$cluster==g),collapse='\t')})
out <- paste(1:length(x.cutnet[[1]]),min(x.cutnet$score)/x.cutnet$score,out,sep='\t')
write(out,paste0('out_',dir()[d],'.txt'))
