#LAB values for Samples selected for Trainning Set:
hsoja.Ymod   #Lab data

#LAB values for Samples selected for Validation Set:
hsoja.Ytest
#Spectra of Samples selected for Trainning Set(dup$model)
dup$model
hsoja.Xmod<-h.soja.snvdt1d[c(80,9,149,136,93,3,36,44,79,119,125,
                      112,100,146,2,30,87,42,106,147,
                      19,90,54,105,82,27,117,5,153,63),]
matplot(wavelength3,t(hsoja.Xmod),type="l",xlab="Wavelength",ylab="Log 1/R",
                 col="blue")
class.mod<-rep("TRAIN",30)
#Spectra of Samples selected for Validation Set (dup$test)
dup$test
hsoja.Xtest<-h.soja.snvdt1d[c(101,47,141,145,83,14,16,150,134,72,148,
                       84,8,158,65,10,151,62,59,22,140,29, 
                       109,98,137,122,75,142,111,12),]
par(new=T)
matplot(wavelength3,t(hsoja.Xtest),type="l",xlab="Wavelength",ylab="Log 1/R",
        col="green")
class.test<-rep("TEST",30)
hsoja.mod<-data.frame(hsoja.Ymod,Class=I(class.mod),NIR=I(hsoja.Xmod))
hsoja.test<-data.frame(hsoja.Ytest,Class=I(class.test),NIR=I(hsoja.Xtest))
#Combining Training and Test Set into a Dataframe.
hsoja<-rbind(hsoja.mod,hsoja.test)

##################  Developing the model for Protein with the Training Set
library(pls)
mod1prot<-plsr(hsoja.mod$Protein~hsoja.mod$NIR,data=hsoja.mod,ncomp=5,validation="LOO")
summary(mod1prot)
#Checking performance with the Training Set (looking for outliers)
hsoja.modpred.prot<-as.numeric(predict(mod1prot,ncomp=5,newdata=hsoja.mod$NIR))
monitor10ftest(hsoja.Ymod$ID,hsoja.modpred.prot,hsoja.Ymod$Protein)
#Looking to the performance with the Test Set
hsoja.testpred.prot<-as.numeric(predict(mod1prot,ncomp=5,newdata=hsoja.test$NIR))
monitor10ftest(hsoja.Ytest$ID,hsoja.testpred.prot,hsoja.Ytest$Protein)
##################  Developing the model for Protein with the Training + Test Set
mod2prot<-plsr(hsoja$Protein~hsoja$NIR,data=hsoja,ncomp=8,validation="LOO")
summary(mod2prot)
#Validating with the rest of the samples
hsoja.rest<-h.soja.snvdt1d[-c(80,9,149,136,93,3,36,44,79,119,125,
                              112,100,146,2,30,87,42,106,147,
                              19,90,54,105,82,27,117,5,153,63,
                              101,47,141,145,83,14,16,150,134,
                              72,148,84,8,158,65,10,151,62,59,22,
                              140,29,109,98,137,122,75,142,111,12),]
par(new=T)
matplot(wavelength3,t(hsoja.rest),type="l",xlab="Wavelength",ylab="Log 1/R",
        col="orange")
hsoja.rest.pred<-as.numeric(predict(mod2prot,ncomp=8,newdata=hsoja.rest))
monitor<-monitor10ftest(hsoja.Yrest$ID,hsoja.rest.pred,
                        hsoja.Yrest$Protein)
#We use the output Selected for a new prediction
