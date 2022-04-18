rm(list=ls())
##############   GENE HEATMAP
library(gplots)
#library(ComplexHeatmap)
#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
setwd("D:\\her2\\fig")
FS_matrix<-read.table("F-S_matrix.txt",header=T,sep = "\t", quote = "")
#FS_matrix<-read.table("F-S_matrix_maxstandardization.txt",header=T,sep = "\t", quote = "")
posname<-read.table("pos_names.txt",header=T,sep = "\t", quote = "")
negname<-read.table("neg_names.txt",header=T,sep = "\t", quote = "")

posindex<-unlist(lapply(posname[,1],function(x){which(colnames(FS_matrix)==x)}))
negindex<-unlist(lapply(negname[,1],function(x){which(colnames(FS_matrix)==x)}))
DEGexp1<-FS_matrix[,c(posindex,negindex)]
rownames(DEGexp1)<-FS_matrix[,1]
for (j in 1:nrow(DEGexp1)) {
  U_test<-wilcox.test(as.numeric(DEGexp1[j,1:length(posindex)]),as.numeric(DEGexp1[j,(length(posindex)+1):ncol(DEGexp1)]))
  pvalue<-U_test[[3]]
  DEGexp1[j,321]<-pvalue
}
indexFO<-c()
indexglcm<-c()
indexglrlm<-c()
indexglszm<-c()
indexgldm<-c()
for(i in 1:nrow(DEGexp1)){
  type<-unlist(strsplit(rownames(DEGexp1)[i], "_"))[2]
  if(type=="firstorder"){
    indexFO<-c(indexFO,i)
  }else if(type=="glcm"){
    indexglcm<-c(indexglcm,i)
  }else if(type=="glrlm"){
    indexglrlm<-c(indexglrlm,i)
  }else if(type=="glszm"){
    indexglszm<-c(indexglszm,i)
  }else if(type=="gldm"){
    indexgldm<-c(indexgldm,i)
  }
}

DEGexp1[indexFO,322]<- "firstorder"
DEGexp1[indexglcm,322]<- "glcm"
DEGexp1[indexglrlm,322]<- "glrlm"
DEGexp1[indexglszm,322]<- "glszm"
DEGexp1[indexgldm,322]<- "gldm"
m=12
gene2<-DEGexp1[,c(order(DEGexp1[m,1:length(posindex)]) ,(order(DEGexp1[m,(length(posindex)+1):320]) +length(posindex)),321,322)]
DEGexp1<-gene2[order(gene2$V321),]


test<-as.matrix(DEGexp1[,1:320])
annotation_row = data.frame( GeneClass = factor(DEGexp1[,322]))
rownames(annotation_row) = rownames(DEGexp1)
colnames(annotation_row)=c("Class")
labels_row = ifelse(DEGexp1[,321]<0.05, "*", "")
Genus=rep(c("pos","neg"),c(length(posindex),length(negindex)))
group_genus=data.frame(Genus)
rownames(group_genus)=colnames(test)

colors=list(Genus=c(pos="#C96D6D", neg="#5FA09B"),
            Class = c(firstorder= "#DCD78D", 
                      glcm= "#CFDB79",
                      glrlm= "#9DC874", 
                      glszm= "#69A4A6", 
                      gldm= "#589DBF"))

bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01)) # ×öÈÈÍ¼:


p1 = pheatmap(test, 
              fontsize = 5,
              fontsize_row = 5, 
              cutree_rows = 2,  
              scale ="row",
              annotation_col=group_genus, 
              annotation_row=annotation_row, 
              annotation_colors=colors,
              show_rownames = T,
              show_colnames  = F,
              labels_row = labels_row, 
              cluster_rows = T,
              cluster_cols = F,
              color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                        colorRampPalette(colors = c("white","red"))(length(bk)/2)),
              legend_breaks=seq(-4,4,2), breaks=bk)
#color = colorRampPalette(c("blue", "white", "red"))(200))
test2<-as.matrix(DEGexp1[DEGexp1$V321<0.05,1:320])

p2 = pheatmap(test2, 
              fontsize = 5,
              fontsize_row = 5, 
              cutree_rows = 2,  
              scale ="row",
              annotation_col=group_genus, 
              annotation_row=annotation_row, 
              annotation_colors=colors,
              show_rownames = T,
              show_colnames  = F,
              labels_row = labels_row, 
              cluster_rows = T,
              cluster_cols = F,
              color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                        colorRampPalette(colors = c("white","red"))(length(bk)/2)),
              legend_breaks=seq(-4,4,2), breaks=bk)

###################   
#install.packages("ggplot2")
library("ggplot2")
DEfeature<-as.data.frame(t(DEGexp1[DEGexp1$V321<0.05,1:320]))

DEfeature1<-cbind(DEfeature, rep(c("pos","neg"),c(length(posindex),length(negindex))))
colnames(DEfeature1)<-c(colnames(DEfeature), "group")

for (m in 1:8) {
  pdf(paste(colnames(DEfeature1)[m],".pdf",sep = ""))
  data<-DEfeature1[,c(m,9)]
  colnames(data)<-c("V1","group")
  ggplot(data = data, aes(x = group, y =V1, fill = group)) +
    geom_violin()+
    ylab(colnames(DEfeature1)[m])+
    geom_boxplot(width=0.1,
                 position = position_dodge(0.9))+
    theme_classic()+
    scale_fill_manual(values = c(pos="#C96D6D",neg="#5FA09B"))
  dev.off()
}
