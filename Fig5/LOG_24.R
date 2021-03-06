rm(list=ls())
#install.packages("pROC")
library(pROC)
fdata1<-read.table("G:\\HER2\\fig\\5\\fdata_24F.txt",header = T,sep = "\t")
fdata1$V25<-factor(fdata1$V25,
                   levels = c(0,1),
                   labels = c("0","1"))
library(e1071)

set.seed(1234)
seeds<-sample(1:10000000,10000,replace=F) 
acc<-c()
roc<-c()
for (i in 1:10000) {
  set.seed(seeds[i])
  ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
  train<-as.data.frame( fdata1[ind,])
  test<-as.data.frame( fdata1[-ind,])
  model1<-glm(V25~.,data=train,family = binomial())
  prediction1 <- predict(model1, test[,1:24])
  prediction1<-ifelse(prediction1>=0,"1","0") 
  ta<-as.data.frame(table(prediction1,test$V25))
  acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
  roc1<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
  roc<-c(roc,roc1$auc)
  acc<-c(acc,acc1)
  print(i)
  
}

max(acc)
max(roc)

write.table(roc,"G:\\her2\\fig\\5\\ROC\\F8_LOG_ROC.txt",col.names = T, row.names = F,sep = "\t" ,append = FALSE, quote = F)
roc<- read.table("G:\\HER2\\fig\\5\\ROC\\F8_LOG_ROC.txt",header = T,sep = "\t")

i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
model1<-glm(V25~.,data=train,family = binomial())
prediction1 <- predict(model1, test[,1:24])
prediction1<-ifelse(prediction1>=0,"1","0") 
ta<-as.data.frame(table(prediction1,test$V25))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc_LOG<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))

