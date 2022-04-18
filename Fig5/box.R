rm(list=ls())
setwd("D:\\her2\\fig\\5\\ROC")
a<-list.files()
c_name<-c("SVM","RF","knn","bayes","LOG","DTREE","ann")
library("RColorBrewer")
library("ggplot2")
library("ggpubr")
data_out<-c()
for (i in 1:4) {   ##lisanxing 
  #for (i in 5:7) {
  classifter_8<-cbind(read.table(paste("F8_",c_name[i],"_ROC.txt",sep = "") ,header = T,sep = "\t"),as.character(c_name[i]) ,"Features_8")
  classifter_24<-cbind(read.table(paste("F24_",c_name[i],"_ROC.txt",sep = "") ,header = T,sep = "\t"),as.character(c_name[i]) ,"Features_24")
 colnames(classifter_8)<-c("auc","classifier","Feature")
 colnames(classifter_24)<-c("auc","classifier","Feature")
 data<-rbind(classifter_8,classifter_24)
  data_out<-rbind(data_out,data)

}

ggplot(data = data_out, aes(x = classifier, y = auc, fill = Feature)) +
  geom_boxplot(alpha=1,) +
  #facet_wrap(~Feature)+
  scale_y_continuous(name = "auc")+
  scale_x_discrete(name = "classifiers") +
  ggtitle("Boxplot of classifiers") +
  theme_bw() +
 scale_fill_manual(values = c("#1E709D","#D07237"))+
  #scale_fill_manual(values = c("#5E9B97","#6B50A0"))+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) +
  stat_compare_means(aes(label = ..p.signif..),method = "t.test")
