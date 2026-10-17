.libPaths(c('/workspace/scratch/4668db2fbf57/validation/R-library', .libPaths()))
for (f in list.files('/workspace/scratch/4668db2fbf57/spca/R', full.names=TRUE)) source(f)
set.seed(802)
X <- matrix(rnorm(18*6),18,6)
X[,2:6] <- X[,2:6] + .4*X[,1:5]
B <- matrix(rnorm(12),6,2)/3
M <- crossprod(X,X%*%B)
cat('Cross-product singular values:',svd(M,nu=0,nv=0)$d,'\n')
value <- function(B) .5*sum(X^2)+.5*sum((X%*%B)^2)-sum(svd(crossprod(X,X%*%B),nu=0,nv=0)$d)
gradients <- lapply(c(1e-4,1e-5,1e-6),function(h){g<-B; for(j in seq_along(B)){plus<-minus<-B;plus[j]<-plus[j]+h;minus[j]<-minus[j]-h;g[j]<-(value(plus)-value(minus))/(2*h)};g})
cat('FD gradients maximum differences, h=1e-5 versus 1e-4/1e-6:',max(abs(gradients[[2]]-gradients[[1]])),max(abs(gradients[[2]]-gradients[[3]])),'\n')
L2 <- svd(X,nu=0,nv=0)$d[1]^2
step <- 1/(L2*1.2)
z <- B-step*(gradients[[2]]+.2*L2*B)
cat('Minimum soft-threshold distance:',min(abs(abs(z)-step*.035*L2)),'\n')
cat('Minimum hard-threshold cost distance:',min(abs(.5*z^2-step*.035*L2)),'\n')
args <- list(X=X,k=2,alpha=.035,beta=.2,center=FALSE,B_init=B,max_iter=1,verbose=FALSE)
for(penalty in c('elasticnet','l0')){args$penalty<-penalty;fit<-do.call(spca,args);expected<-if(penalty=='elasticnet')sign(z)*pmax(abs(z)-step*.035*L2,0)else z*(.5*z^2>step*.035*L2);cat(penalty,'max coefficient error',max(abs(fit$loadings-expected)),'\n')}
