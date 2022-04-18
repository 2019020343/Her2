rm(list=ls())
###################   A" barplot  ###############
#install.packages("ggplot2")
library("ggplot2")
data1<-read.table("G:\\HER2\\gene\\gene-geo.txt",header=F,sep = "\t", quote = "\"'")
colnames(data1)<-c("type","value","name")
ggplot(data = data1,aes(x = type,y = value))+
  geom_bar(stat = "identity", aes(fill=name),
           position = 'dodge')+
  theme_classic()+
  scale_fill_manual(values = c("#5FA09B","#C96D6D"))
  
########### B heatmap: GSE45827 ############
library(gplots)
#library(ComplexHeatmap)
#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
posgene<-read.table("D:\\HER2\\gene\\result\\11\\posgene.txt",header=T,sep = "\t", quote = "")
GPL5<-read.table("G:\\HER2\\gene\\result\\pos\\GSE162228\\GPL570-55999.txt",header=T,sep = "\t", quote = "")
GPL5<-GPL5[,c(1,11)]
### GSE45827
geneexp<-read.table("G:\\HER2\\gene\\result\\pos\\GSE45827\\GSE45827_series_matrix.txt",header=T,sep = "\t", quote = "\"'")
GSE<-read.table("G:\\HER2\\gene\\result\\pos\\GSE45827\\GSE45827.txt",header=F,sep = "\t", quote = "\"'")
gsematch<-read.table("G:\\HER2\\gene\\result\\pos\\GSE45827\\gsematch.txt",header=F,sep = "\t", quote = "\"'")
### GSE162228
geneexp<-read.table("G:\\HER2\\gene\\result\\pos\\GSE162228\\GSE162228_series_matrix.txt",header=T,sep = "\t", quote = "\"'")
GSE<-read.table("G:\\HER2\\gene\\result\\pos\\GSE162228\\GSE162228.txt",header=F,sep = "\t", quote = "\"'")
gsematch<-read.table("G:\\HER2\\gene\\result\\pos\\GSE162228\\gsematch.txt",header=F,sep = "\t", quote = "\"'")
### GSE129559
GPL5<-read.table("G:\\HER2\\gene\\result\\pos\\GSE129559\\GPL96-57554.txt",header=T,sep = "\t", quote = "")
GPL5<-GPL5[,c(1,11)]
geneexp<-read.table("G:\\HER2\\gene\\result\\pos\\GSE129559\\GSE129559_series_matrix.txt",header=T,sep = "\t", quote = "\"'")
GSE<-read.table("G:\\HER2\\gene\\result\\pos\\GSE129559\\GSE129559.txt",header=F,sep = "\t", quote = "\"'")
gsematch<-read.table("G:\\HER2\\gene\\result\\pos\\GSE129559\\gsematch.txt",header=F,sep = "\t", quote = "\"'")
### GSE81538
GPL5<-read.csv("D:\\her2\\gene\\result\\pos\\val\\GSE81538\\GSE81538_map_transcriptID_geneSymbol.csv",header=T,sep = "\t", quote = "",fill = TRUE)
geneexp<-read.csv("D:\\her2\\gene\\result\\pos\\val\\GSE81538\\GSE81538_transcript_expression_405.csv",header=T,sep = ",", quote = "\"'")
GSE<-read.table("D:\\her2\\gene\\result\\pos\\val\\GSE81538\\GSE81538.txt",header=F,sep = "\t", quote = "\"'")
gsematch<-read.table("D:\\her2\\gene\\result\\pos\\val\\GSE81538\\gsematch.txt",header=F,sep = "\t", quote = "\"'")


#sample name1
posgse<-GSE[1,GSE[2,]=="pos"]
neggse<-GSE[1,GSE[2,]=="neg"]
#sample col_index
posgse_index<-unlist(lapply(posgse[1,], function(x){which(gsematch[1,]==as.character( x))})) 
neggse_index<-unlist(lapply(neggse[1,], function(x){which(gsematch[1,]==as.character( x))})) 
#sample gene expression :pos-neg
geneexp<-geneexp[,c(1,posgse_index,neggse_index)]
#GPL5 :ID  GENE_SYMBOL

#TOP1000 gene related with features 
myfun2<-function(x){unlist(strsplit(as.character(x),"-"))[2]}
gene<-as.data.frame(table(unlist(lapply(posgene[,1],myfun2))))
gene1<-gene[order(-gene$Freq)[1:100],]

#ID of TOP1000 gene 
ID_out<-unlist(lapply(gene1[,1], function(x){as.character(GPL5[grep(x,GPL5[,2]),1]) })) 
ID_out<-as.data.frame(unique(ID_out))
GPL5_gene1<-merge(ID_out,geneexp,by.x ="unique(ID_out)" ,by.y = "X")
#geneexp<-read.table("G:\\HER2\\gene\\result\\pos\\GSE45827\\GSE45827_mean.txt",header=T,sep = "\t", quote = "")
#geneexp<-read.table("G:\\HER2\\gene\\result\\pos\\GSE162228\\GSE45827_mean.txt",header=T,sep = "\t", quote = "")
#geneexp<-read.table("G:\\HER2\\gene\\result\\pos\\GSE129559\\GSE45827_mean.txt",header=T,sep = "\t", quote = "")
#geneexp_1000<-merge(geneexp,gene1,by.x ="V1" ,by.y = "Var1")
test<-as.matrix(GPL5_gene1[,2:ncol(GPL5_gene1)])
test<-test[unlist(apply(test, 1, function(x){length(which(x=="0"))!=401})),]
test<-scale(test)

Sample_type=rep(c("pos","neg"),c(ncol(posgse),ncol(neggse)))
Sample_type=data.frame(Sample_type)
rownames(Sample_type)=colnames(test)

colors=list(Sample_type=c(pos="#C96D6D", neg="#5FA09B"))
bk <- c(seq(-1,0,by=0.01),seq(0.1,24,by=0.01)) # ×öÈÈÍ¼:

pheatmap(test, 
              fontsize = 5,
              fontsize_row = 5, 
              cutree_rows = 2,  
              #scale ="row",
              annotation_col=Sample_type, 
              #annotation_row=annotation_row, 
              annotation_colors=colors,
              show_rownames = T,
              show_colnames  = F,
              #labels_row = labels_row, 
              cluster_rows = T,
              cluster_cols = F,
              color = colorRampPalette(c("blue", "white", "red"))(200))

         
         
              color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                                  colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                        legend_breaks=seq(-10,10,2), 
                        breaks=bk)
              
              
              
              

              #color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
              #          colorRampPalette(colors = c("white","red"))(length(bk)/2)),
              #legend_breaks=seq(-10,10,2), 
              #breaks=bk)
#color = colorRampPalette(c("blue", "white", "red"))(200),









