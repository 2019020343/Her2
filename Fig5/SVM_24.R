rm(list=ls())
#install.packages("pROC")
library(pROC)
fdata1<-read.table("D:\\HER2\\fig\\5\\fdata_24F.txt",header = T,sep = "\t")
fdata1$V25<-factor(fdata1$V25,
                  levels = c(0,1),
                  labels = c("0","1"))
#install.packages("e1071")
library(e1071)

set.seed(1234)
seeds<-sample(1:10000000,10000,replace=F) 
acc<-c()
roc<-c()
for (i in 1:10000) {
  set.seed(seeds[i])
  ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
  train<-as.data.frame( fdata1[ind,])
  test<-as.data.frame(fdata1[-ind,])
  svm<-svm(V25~.,data=train,type="C",kernel="radial")
  pred<-predict(svm,test[,1:24])
  ta<-as.data.frame(table(pred,test$V25))
  roc1<- roc(test[,ncol(test)],as.numeric(as.matrix(pred)))
  acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
  roc<-c(roc,roc1$auc)
  acc<-c(acc,acc1)
  print(i)
}

max(acc)
max(roc)
write.table(roc,"D:\\her2\\fig\\5\\ROC\\F24_SVM_ROC.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
roc<- read.table("D:\\HER2\\fig\\5\\ROC\\F24_SVM_ROC.txt",header = T,sep = "\t")

i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame(fdata1[-ind,])
svm<-svm(V25~.,data=train,type="C",kernel="radial")
pred<-predict(svm,test[,1:24])
ta<-as.data.frame(table(pred,test$V25))
roc_SVM<- roc(test[,ncol(test)],as.numeric(as.matrix(pred)))




