rm(list=ls())
library(gplots)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)
setwd("D:\\HER2\\gene\\result\\pos\\val")
posgene<-read.table("posgene.txt",header=T,sep = "\t", quote = "\"'")
neggene<-read.table("neggene.txt",header=T,sep = "\t", quote = "\"'")

myfun1<-function(x){unlist(strsplit(as.character(x),"-"))[1]}
myfun2<-function(x){unlist(strsplit(as.character(x),"-"))[2]}

feature<-unique(unlist(lapply(posgene[,1],myfun1)))

fname<-c("original_glszm_SizeZoneNonUniformityNormalized",
         "original_glszm_ZoneEntropy",
         "original_glrlm_ShortRunHighGrayLevelEmphasis",
         "original_glszm_SizeZoneNonUniformity",
         "original_glszm_SmallAreaEmphasis",
         "original_glrlm_ShortRunEmphasis",
         "original_glrlm_RunLengthNonUniformityNormalized",
         "original_glszm_SmallAreaHighGrayLevelEmphasis")
feature1<-unique(unlist(lapply(fname,function(x){which(feature==x)})))
ontology<-c("BP","MF","CC","kegg")
dataout_GR<-as.data.frame(matrix(0,ncol = 8,nrow = 4))
colnames(dataout_GR)<-fname
rownames(dataout_GR)<-ontology
dataout_p<-as.data.frame(matrix(0,ncol = 8,nrow = 4))
colnames(dataout_p)<-fname
rownames(dataout_p)<-ontology
ontology_i_out<-c()
for (i in 1:length(feature1)) {
  feature_i<-read.table(paste("D:\\HER2\\F", feature1[i], ".txt", sep = ""),header=T,sep = "\t", quote = "\"'")
  feature_bar<-feature_i[1:10,]
  feature_bar[,11]<--log10(feature_i$pvalue[1:10])
  ggplot(feature_bar, aes(V11,Description))+
    #geom_bar(stat = "identity",fill="#ABCECD")+ #green
    #geom_bar(stat = "identity",fill="#C09C60")+ #brown
    #geom_bar(stat = "identity",fill="#6398AD")+ #blue
    geom_bar(stat = "identity",fill="#AF4D47")+ #
    #coord_flip()+
    theme_classic()+
    scale_fill_manual()+
    labs(y=feature[feature1[i]])
  
  
  for (j in 1:length(ontology)) {
      ontology_i<- feature_i[feature_i$ONTOLOGY==ontology[j],] 
      if(nrow(ontology_i)!=0){
        dataout_p[j,i]<-(-log10(mean(ontology_i$pvalue)))
        dataout_GR[j,i]<-mean(ontology_i$Count)/as.numeric( strsplit(as.character( ontology_i$GeneRatio[1]),"/")[[1]][2])
        
        }
      
    }
  }

pheatmap(dataout_GR[1,], 
         color = colorRampPalette(c("#A9D3CF", "#5FA09B"))(200),#BP
         #color = colorRampPalette(c("#E8ACAC", "#C96D6D"))(200),#MF
         #color = colorRampPalette(c("#91D1E5", "#2A758C"))(200),#CC
         #color = colorRampPalette(c("#EFBFA5", "#D97435"))(200),#KEGG
         fontsize = 5,
         fontsize_row = 5, 
         show_rownames = T,
         show_colnames  = T,
         cluster_rows = F,
         cluster_cols = F)
library("reshape")
library("viridis")
library("ggridges")

dataout_GR<-cbind(rownames(dataout_GR),dataout_GR)
dataout_GR1<-melt(dataout_GR,measure.vars = fname, variable_name = "Features", value.name="value")
colnames(dataout_GR1)<-c("name","Features","value")

dataout_p<-cbind(rownames(dataout_p),dataout_p)
dataout_p1<-melt(dataout_p,measure.vars = fname, variable_name = "Features", value.name="value")
colnames(dataout_p1)<-c("name","Features","value")


ggplot(dataout_GR1, aes(x=value, y=name, fill=name))+
  geom_density_ridges(scale=4,alpha=0.8)+
  scale_fill_cyclical(values = c("#5FA09B", "#2A758C","#D97435","#C96D6D"))+
  theme_ridges(grid = FALSE)




