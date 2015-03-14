library(prospectr)
hsoja.X<-as.matrix(read.table("clipboard",header=FALSE))
hsoja.Y<-as.data.frame(read.table("clipboard",header=TRUE))
hsoja.Y<-edit(hsoja.Y)   #Poner NA a los valores cero

wavelength1<-seq(from=1100.0,to=2499.5, by=0.5)
X.ID<-seq(from=1,to=159, by=1)
colnames(hsoja.X)<-wavelength1
rownames(hsoja.X)<-X.ID
plot(as.numeric(colnames(hsoja.X)),hsoja.X[1,],type="l",xlab="wavelength",ylab="Absorbance")
matplot(wavelength1,t(hsoja.X),type="l")
hsoja.X.bin<-binning(hsoja.X,bin.size=4)
wavelength2<-seq(from=1100,to=2498, by=2)
colnames(hsoja.X.bin)<-wavelength2
matplot(wavelength2,t(hsoja.X.bin),type="l",ylab="",xlab="Wavelength")
points(as.numeric(colnames(hsoja.X.bin)),hsoja.X.bin[1,],pch=2)

######## Apply anti-scatter mathtreatments ###########################
matplot(wavelength1,t(hsoja.X),type="l",ylab="",xlab="Wavelength",col="black")
par(new=T)
#SNV
hsoja.Xsnv<-standardNormalVariate(X=hsoja.X)
matplot(wavelength1,t(hsoja.Xsnv),type="l",xaxt="n",yaxt="n",ylab="",xlab="",col="blue")
#SNV and DT
hsoja.Xsnvdt<-detrend(X=hsoja.X,wav=as.numeric(colnames(hsoja.X)))
par(new=T)
matplot(wavelength1,t(hsoja.Xsnvdt),type="l",xaxt="n",yaxt="n",xlab="",ylab="",col="red")

######## Apply derivatives ###########################
h.soja.snvdt1d<-t(diff(t(hsoja.Xsnvdt),differences=1,lag=16))
wavelength3<-seq(from=1108,to=2499.5, by=0.5)
par(new=T)
matplot(wavelength3,t(h.soja.snvdt1d),type="l",xaxt="n",yaxt="n",xlab="",ylab="",col="green")

#####################  All the sequence  ###########################3
matplot(wavelength1,t(hsoja.X),type="l",
        ylab="",xlab="Wavelength",col="black")
par(new=T)
matplot(wavelength1,t(hsoja.Xsnv),type="l",
        xaxt="n",yaxt="n",ylab="",xlab="",col=4)
par(new=T)
matplot(wavelength1,t(hsoja.Xsnvdt),type="l",
        xaxt="n",yaxt="n",xlab="",ylab="",col=2)
par(new=T)
matplot(wavelength3,t(h.soja.snvdt1d),type="l",xaxt="n",
        yaxt="n",xlab="",ylab="",col="green")
legend("topleft",legend=c("Raw","SNV","SNV+DT","1ª Deriv"),
       lty=c(1,1),col=c("black","blue","red","green"))

####################### Shenk  ##################################
library(prospectr)
shenk<-shenkWest(X=h.soja.snvdt1d,d.min=0.2,pc=6)
plot(shenk$pc[,1],shenk$pc[,2],xlab="PC1",ylab="PC2")
points(shenk$pc[shenk$model,1],shenk$pc[shenk$model,2],pch=19,col="red")
################   DUPLEX  #####################################
library(prospectr)
dup<-duplex(X=h.soja.snvdt1d,k=30,metric="mahal",pc=3)
par(mfrow=c(3,1),ps=14)
plot(dup$pc[,1],dup$pc[,2],xlab="PC1",ylab="PC2")
points(dup$pc[dup$model,1],dup$pc[dup$model,2],pch=19,col="red")
points(dup$pc[dup$test,1],dup$pc[dup$test,2],pch=19,col="blue")
plot(dup$pc[,1],dup$pc[,3],xlab="PC1",ylab="PC3")
points(dup$pc[dup$model,1],dup$pc[dup$model,3],pch=19,col="red")
points(dup$pc[dup$test,1],dup$pc[dup$test,3],pch=19,col="blue")
plot(dup$pc[,2],dup$pc[,3],xlab="PC2",ylab="PC3")
points(dup$pc[dup$model,2],dup$pc[dup$model,3],pch=19,col="red")
points(dup$pc[dup$test,2],dup$pc[dup$test,3],pch=19,col="blue")

class(dup$model)
hsoja.Ymod<-hsoja.Y[c(80,9,149,136,93,3,36,44,79,119,125,
                      112,100,146,2,30,87,42,106,147,
                      19,90,54,105,82,27,117,5,153,63),]
hsoja.Ytest<-hsoja.Y[c(101,47,141,145,83,14,16,150,134,72,148,
                       84,8,158,65,10,151,62,59,22,140,29, 
                       109,98,137,122,75,142,111,12),]

summary(hsoja.Y)
hist(hsoja.Y$Protein,col="green")
summary(hsoja.Ymod)
hist(hsoja.Ymod$Protein,col="blue")
summary(hsoja.Ytest)
hist(hsoja.Ytest$Protein,col="violet")
