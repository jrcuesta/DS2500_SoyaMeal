dup<-duplex(X=h.soja.snvdt1d,k=30,metric="mahal",pc=3)
par(mfrow=c(2,2),ps=14)
plot(dup$pc[,1],dup$pc[,2],xlab="PC1",ylab="PC2")
points(dup$pc[dup$model,1],dup$pc[dup$model,2],pch=19,col="red")
drawMahal(dup$pc[,1:2],center=apply(dup$pc[,1:2],2,mean),
          covariance=cov(dup$pc[,1:2]),quantile=0.975)
#Once remove the redundants, we can calculate the PCs
#But first we have to subset the 30 selected samples "hsoja.Xmod"
Xmod.prcomp<-prcomp(hsoja.Xmod)
plot(Xmod.prcomp$x[,1],Xmod.prcomp$x[,2],xlab="PC1",ylab="PC2")
library(chemometrics)
drawMahal(Xmod.prcomp$x[,1:2],center=apply(Xmod.prcomp$x[,1:2],2,mean),
          covariance=cov(Xmod.prcomp$x[,1:2]),quantile=0.975,col="red")
#To project the test set into this space, we have to calculate its 
#scores in this new PC space first.
res<-Moutlier(Xmod.prcomp$x[,1:2],quantile=0.975,plot=TRUE)
res<-Moutlier(Xmod.prcomp$x[,2:3],quantile=0.975,plot=TRUE)
