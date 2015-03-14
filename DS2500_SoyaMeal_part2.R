#LAB values for Samples selected for Trainning Set:
hsoja.Ymod   #Lab data

#LAB values for Samples selected for Validation Set:
hsoja.Ytest
#Spectra of Samples selected for Trainning Set(dup$model)
dup$model
hsoja.Xmod<-h.soja.snvdt1d[c(80,9,149,136,93,3,36,44,79,119,125,
                      112,100,146,2,30,87,42,106,147,
                      19,90,54,105,82,27,117,5,153,63),]
class.mod<-rep("TRAIN",30)
#Spectra of Samples selected for Validation Set (dup$test)
dup$test
hsoja.Xtest<-h.soja.snvdt1d[c(101,47,141,145,83,14,16,150,134,72,148,
                       84,8,158,65,10,151,62,59,22,140,29, 
                       109,98,137,122,75,142,111,12),]
class.test<-rep("TEST",30)
hsoja.mod<-data.frame(hsoja.Ymod,Class=I(class.mod),NIR=I(hsoja.Xmod))
hsoja.test<-data.frame(hsoja.Ytest,Class=I(class.test),NIR=I(hsoja.Xtest))
#Combining Training and Test Set into a Dataframe.
hsoja<-rbind(hsoja.mod,hsoja.test)
##################  Developing the model for Protein
library(pls)
mod1prot<-plsr(hsoja.mod$Protein~hsoja.mod$NIR,data=hsoja.mod,ncomp=5,validation="LOO")
summary(mod1prot)
hsoja.testpred.prot<-as.numeric(predict(mod1prot,ncomp=5,newdata=hsoja.test$NIR))
monitor10ftest(hsoja.Ytest$ID,hsoja.testpred.prot,hsoja.Ytest$Protein)
##################  Developing the model for Moisture
mod1moi<-plsr(hsoja.mod$Moisture~hsoja.mod$NIR,data=hsoja.mod,ncomp=5,validation="LOO")
summary(mod1moi)
hsoja.testpred.moi<-as.numeric(predict(mod1moi,ncomp=5,newdata=hsoja.test$NIR))
monitor10ftest(hsoja.Ytest$ID,hsoja.testpred.moi,hsoja.Ytest$Moisture)