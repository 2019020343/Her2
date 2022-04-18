rm(list=ls())
#install.packages("pROC")
library(pROC)
mydata<-read.table("D:\\HER2\\F-S_matrix.txt",header = T,sep = "\t")
neg_names<-read.table("D:\\HER2\\gene\\result\\pos\\GSE129559\\neg_names.txt",header = T,sep = "\t")
pos_names<-read.table("D:\\HER2\\gene\\result\\pos\\GSE129559\\pos_names.txt",header = T,sep = "\t")

#fname<-unique(unlist(lapply(pos_index1_2[,1], myfun1)))
fname<-c("original_glszm_SizeZoneNonUniformityNormalized","original_glszm_ZoneEntropy","original_glrlm_ShortRunHighGrayLevelEmphasis",
         "original_glszm_SizeZoneNonUniformity","original_glszm_SmallAreaEmphasis","original_glrlm_ShortRunEmphasis",
         "original_glrlm_RunLengthNonUniformityNormalized", "original_glszm_SmallAreaHighGrayLevelEmphasis")


f_index<-unlist(lapply(fname, function(x){which(mydata[,1]==x)}))
fdata<-as.data.frame(t(mydata[f_index,]))
colnames(fdata)<-fname
fdata<-fdata[-1,]

neg_index<-unlist(lapply(neg_names[,1], function(x){which(rownames(fdata)==x)}))
pos_index<-unlist(lapply(pos_names[,1], function(x){which(rownames(fdata)==x)}))
fdata[neg_index,ncol(fdata)+1]<-"0"
fdata[pos_index,ncol(fdata)]<-"1"

fdata1<-as.numeric(as.matrix(fdata))
fdata1<-as.data.frame(matrix(fdata1 ,ncol=9,nrow=320,byrow = F) ) 
fdata1$V9<-factor(fdata1$V9,
                  levels = c(0,1),
                  labels = c("0","1"))
write.table(fdata1,"D:\\HER2\\fig\\5\\fdata1.txt", sep = "\t", quote = F, row.names = F)
