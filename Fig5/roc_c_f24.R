rm(list=ls())
#install.packages("pROC")
library(pROC)
fdata1<-read.table("D:\\HER2\\fig\\5\\fdata_24F.txt",header = T,sep = "\t")
fdata1$V25<-factor(fdata1$V25,
                   levels = c(0,1),
                   labels = c("0","1"))
#install.packages("e1071")
library(e1071)
library(kknn)
library(rpart)
library(nnet)
library(randomForest)
set.seed(1234)
seeds<-sample(1:10000000,10000,replace=F)  

roc<- read.table("D:\\HER2\\fig\\5\\ROC\\F24_SVM_ROC.txt",header = T,sep = "\t")

i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame(fdata1[-ind,])
svm<-svm(V25~.,data=train,type="C",kernel="radial")
pred<-predict(svm,test[,1:24])
ta<-as.data.frame(table(pred,test$V25))
roc(test[,25],as.numeric(as.matrix(pred)),plot=T,col="#5972AE",print.auc=T,print.auc.x=0.1,print.auc.y=0.01, main="SVM")


roc<- read.table("D:\\HER2\\fig\\5\\ROC\\F24_bayes_ROC.txt",header = T,sep = "\t")
i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
nb1 <- naiveBayes(V25 ~ ., data =train,laplace=0 )
prediction <- predict(nb1, test[,1:24])
ta<-as.data.frame(table(prediction,test$V25))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc_bayes<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction)))
roc(test$V25,as.numeric(as.matrix(prediction)),plot=T,col="#D97435",print.auc=T,print.auc.x=0.1,print.auc.y=0.01, main="bayes")



roc<- read.table("D:\\HER2\\fig\\5\\ROC\\F24_RF_ROC.txt",header = T,sep = "\t")
i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame(fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
nb1 <- randomForest(V25~.,data = train)
prediction1 <- predict(nb1, test[,1:24])
#prediction1 <- ifelse(prediction[,1]>0.5,0,1)
ta<-as.data.frame(table(prediction1,test$V25))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc1_RF<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
roc(test$V25,as.numeric(as.matrix(prediction1)),plot=T,col="#1C833D",print.auc=T,print.auc.x=0.1,print.auc.y=0.01, main="RandomForest")


roc<- read.table("D:\\HER2\\fig\\5\\ROC\\F24_knn_ROC.txt",header = T,sep = "\t")
i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame(fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
nb1 <- kknn(V25~.,train = train,test=test,distance = 1,kernel = "triangular")
prediction1<-fitted(nb1)
#prediction1 <- ifelse(prediction[,1]>0.5,0,1)
ta<-as.data.frame(table(prediction1,test$V25))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc_knn<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
  roc(test$V25,as.numeric(as.matrix(prediction1)),plot=T,col="#1476A6",print.auc=T,print.auc.x=0.1,print.auc.y=0.05, main="KNN")

############# c  ######################################################################
roc<- read.table("D:\\HER2\\fig\\5\\ROC\\F24_DTREE_ROC.txt",header = T,sep = "\t")

i=which(roc==max(roc))
set.seed(1234)
seeds<-sample(1:10000000,10000,replace=F)
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
nb1 <- rpart(V25~.,data = train, method = "class")
prediction <- predict(nb1, test[,1:24])
prediction1 <- ifelse(prediction[,1]>0.5,0,1)
ta<-as.data.frame(table(prediction1,test$V25))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc_DTREE<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
roc(test$V25,as.numeric(as.matrix(prediction[,1])),plot=T,col="#5FA09B",print.auc=T,print.auc.x=0.1,print.auc.y=0.05, main="Decision Tree")




roc<- read.table("D:\\HER2\\fig\\5\\ROC\\F24_LOG_ROC.txt",header = T,sep = "\t")

i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
model1<-glm(V25~.,data=train,family = binomial())
prediction <- predict(model1, test[,1:24])
prediction1<-ifelse(prediction>=0,"1","0") 
ta<-as.data.frame(table(prediction1,test$V25))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc_LOG<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
roc(test$V25,as.numeric(as.matrix(prediction)),plot=T,col="#C4532C",print.auc=T,print.auc.x=0.1,print.auc.y=0.05, main="Logistic")




roc<- read.table("D:\\HER2\\fig\\5\\ROC\\F24_ann_ROC.txt",header = T,sep = "\t")

i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame(fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
nb1 <-nnet(V25~.,data=train,  linout=F,size=5,decay=5e-4,maxit=100,trace=F,rang=0.7)
prediction1 <- predict(nb1, test[,1:24])

#prediction1 <- ifelse(prediction[,1]>1,1,0)
ta<-as.data.frame(table(prediction1,test$V25))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc1<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)),plot=T,col="#4C2B83",print.auc=T,print.auc.x=0.1,print.auc.y=0.05, main="ANN")



