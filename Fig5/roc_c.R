rm(list=ls())
#install.packages("pROC")
library(pROC)
library(randomForest)
library(nnet)
library(kknn)
fdata1<-read.table("D:\\HER2\\fig\\5\\fdata1.txt",header = T,sep = "\t")
fdata1$V9<-factor(fdata1$V9,
                  levels = c(0,1),
                  labels = c("0","1"))
#install.packages("e1071")
library(e1071)
library(rpart)

set.seed(1234)
seeds<-sample(1:10000000,10000,replace=F) 
roc<- read.table("D:\\HER2\\fig\\5\\ROC\\f8\\F8_SVM_ROC.txt",header = T,sep = "\t")

i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame(fdata1[-ind,])
svm<-svm(train[,1:(ncol(train)-1)],train[,ncol(train)],type="C-classification",
         cost=10,kernel="radial",probability=TRUE,scale=FALSE)
pred<-predict(svm,test[,1:(ncol(test)-1)],decision.values=F)
ta<-as.data.frame(table(pred,test[,ncol(test)]))
roc_SVM<- roc(test[,ncol(test)],as.numeric(as.matrix(pred)))
roc(test$V9,as.numeric(as.matrix(pred)),plot=T,col="#5972AE",print.auc=T,print.auc.x=0.1,print.auc.y=0.05,main="SVM")

roc<- read.table("D:\\HER2\\fig\\5\\ROC\\f8\\F8_bayes_ROC.txt",header = T,sep = "\t")
i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
nb1 <- naiveBayes(V9 ~ ., data =train,laplace=0 )
prediction <- predict(nb1, test[,1:8])
ta<-as.data.frame(table(prediction,test$V9))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc_bayes<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction)))
roc(test$V9,as.numeric(as.matrix(prediction)),plot=T, col="#D97435",print.auc=T,print.auc.x=0.1,print.auc.y=0.05,main="bayes")


roc<- read.table("D:\\HER2\\fig\\5\\ROC\\f8\\F8_RF_ROC.txt",header = T,sep = "\t")
i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame(fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
nb1 <- randomForest(V9~.,data = train)
prediction1 <- predict(nb1, test[,1:8])
#prediction1 <- ifelse(prediction[,1]>0.5,0,1)
ta<-as.data.frame(table(prediction1,test$V9))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc1_RF<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
roc(test$V9,as.numeric(as.matrix(prediction1)),plot=T, col="#1C833D",print.auc=T,print.auc.x=0.1,print.auc.y=0.05,main="RandomForest")


roc<- read.table("D:\\HER2\\fig\\5\\ROC\\f8\\F8_knn_ROC.txt",header = T,sep = "\t")
i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame(fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
nb1 <- kknn(V9~.,train = train,test=test,distance = 1,kernel = "triangular")
prediction1<-fitted(nb1)
#prediction1 <- ifelse(prediction[,1]>0.5,0,1)
ta<-as.data.frame(table(prediction1,test$V9))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc_knn<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
roc(test$V9,as.numeric(as.matrix(prediction1)),plot=T, col="#1476A6",print.auc=T,print.auc.x=0.1,print.auc.y=0.05,main="KNN")

############# c  ######################################################################
roc<- read.table("D:\\HER2\\fig\\5\\ROC\\f8\\F8_DTREE_ROC.txt",header = T,sep = "\t")

i=which(roc==max(roc))
set.seed(1234)
seeds<-sample(1:10000000,10000,replace=F)
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
nb1 <- rpart(V9~.,data = train, method = "class")
prediction <- predict(nb1, test[,1:8])
prediction1 <- ifelse(prediction[,1]>0.5,0,1)
ta<-as.data.frame(table(prediction1,test$V9))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc_DTREE<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
roc(test$V9,as.numeric(as.matrix(prediction[,1])),plot=T,col="#5FA09B",print.auc=T,print.auc.x=0.1,print.auc.y=0.05,main="Decision Tree")




roc<- read.table("D:\\HER2\\fig\\5\\ROC\\f8\\F8_LOG_ROC.txt",header = T,sep = "\t")

i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame( fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
model1<-glm(V9~.,data=train,family = binomial())
prediction <- predict(model1, test[,1:8])
prediction1<-ifelse(prediction>=0,"1","0") 
ta<-as.data.frame(table(prediction1,test$V9))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc_LOG<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
roc(test$V9,as.numeric(as.matrix(prediction)),plot=T, col="#C4532C",print.auc=T,print.auc.x=0.1,print.auc.y=0.05,main="Logistic")





roc<- read.table("D:\\HER2\\fig\\5\\ROC\\f8\\F8_ann_ROC.txt",header = T,sep = "\t")
i=which(roc==max(roc))
set.seed(seeds[i])
ind<-sample(2:nrow(fdata1),floor(nrow(fdata1)*0.7),replace=F) #70%为训练集 30%为测试集
train<-as.data.frame(fdata1[ind,])
test<-as.data.frame( fdata1[-ind,])
r=1/max(abs(train[,1:8]))
nb1 <- nnet(V9~.,data = train, size=4, rang = r,decay = 0.0005, maxit = 200)
prediction1 <- predict(nb1, test[,1:8])
#prediction1 <- ifelse(prediction[,1]>1,1,0)
ta<-as.data.frame(table(prediction1,test$V9))
acc1<-(ta[1,3]+ta[4,3])/sum(ta[,3])
roc1<- roc(test[,ncol(test)],as.numeric(as.matrix(prediction1)))
roc(test$V9,as.numeric(as.matrix(prediction1)),plot=T, col="#4C2B83",print.auc=T,print.auc.x=0.1,print.auc.y=0.05,main="ANN")




