rm(list=ls())
library(gplots)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)
pos_index1_2 <- read.delim('D:\\her2\\gene\\result\\pos\\val\\pos_index1_2.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
myfun1<-function(x){unlist(strsplit(x,"-",fixed = T))[1]}
myfun2<-function(x){unlist(strsplit(x,"-",fixed = T))[2]}

#fname<-unique(unlist(lapply(pos_index1_2[,1], myfun1)))
fname<-c("original_glszm_SizeZoneNonUniformityNormalized",
         "original_glszm_ZoneEntropy",
         "original_glrlm_ShortRunHighGrayLevelEmphasis",
         "original_glszm_SizeZoneNonUniformity",
         "original_glszm_SmallAreaEmphasis",
         "original_glrlm_ShortRunEmphasis",
         "original_glrlm_RunLengthNonUniformityNormalized",
         "original_glszm_SmallAreaHighGrayLevelEmphasis")

fgname<-unique(pos_index1_2[unlist(lapply(fname, function(x){grep(x,pos_index1_2[,1])})),1])

#fgname<-unique(pos_index1_2[,1])
gname<-as.data.frame(table(unlist(lapply(fgname, myfun2))))
gname_100<-gname[order(-gname$Freq)[1:25],]

#re1 = melt(data = weight,id.vars=c("s"),variable.name="feature",value.name="y")
gname_number<-c()
FG_pr<-c()
for (i in 1:length(fname)) {
  feature_i<-pos_index1_2[unlist(lapply(pos_index1_2[,1], myfun1))==fname[i],]
  gname_i<-cbind(fname[i], length(unique(unlist(lapply(feature_i[,1], myfun2)))))
  gname_number<-c(gname_number,gname_i)
  for (j in 1:nrow(gname_100)) {
    gene_j<-feature_i[unlist(lapply(feature_i[,1], myfun2))==gname_100[j,1],]
    
    if(nrow(gene_j)!=0){
      p_index<-as.matrix(gene_j[,seq(3,51,2)]) 
      R_index<-as.matrix(gene_j[,seq(2,51,2)])
      P_mean<-mean(as.numeric(-log10(p_index[which(p_index!="NA")])))
      R_mean<-mean(as.numeric(abs(R_index[which(R_index!="NA")])))
    }else{
      P_mean<-0
      R_mean<-0
    }
    data<-cbind(fname[i], as.character(gname_100[j,1]), R_mean, P_mean)
    FG_pr<-rbind(FG_pr,data)
    print(paste(i,"-",j,sep = ""))
  }
}


FG_pr<-as.data.frame(FG_pr)
FG_pr[as.numeric(as.matrix( FG_pr$P_mean))<0.01,5]<-"most"
FG_pr[0.01<as.numeric(as.matrix( FG_pr$P_mean)) & as.numeric(as.matrix( FG_pr$P_mean))<0.03,5]<-"medium"
FG_pr[as.numeric(as.matrix( FG_pr$P_mean))>0.03,5]<-"least"

write.table(FG_pr,"G:\\HER2\\fig\\2\\FG_pr.txt", sep = "\t", quote = F, row.names = F)

ggplot(FG_pr, aes(V1,V2,col=as.numeric(as.matrix(R_mean)) ))+
  geom_tile(color="grey85", fill="white", size=1)+
  geom_point(aes(size = as.numeric(as.matrix(P_mean))))+
  scale_colour_gradient2(high = "#F27D58",mid = "blue",low = "blue")+coord_flip()+
  theme(axis.text.x=element_text(angle = 45))


gname_number<-matrix(gname_number,ncol = 2,nrow = 8,byrow = T)
gname_number<-as.data.frame(gname_number)



p2<-ggplot(gname_number, aes(V1,as.numeric( as.matrix(V2) ) ))+
  geom_bar(stat = "identity",fill="#ABCECD")+
  coord_flip()+
  theme_classic()+
  scale_fill_manual()

p2



p3<-ggplot(gname_100, aes(Var1,Freq))+
  geom_bar(stat = "identity",fill="#669EB5")+
  theme_classic()+
  scale_fill_manual()+
  theme(axis.text.x=element_text(angle = 45))

p3


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
  feature_i<-read.table(paste("D:\\her2\\gene\\result\\pos\\val\\feature_", feature1[i], "_pos.txt", sep = ""),header=T,sep = "\t", quote = "\"'")
  feature_bar<-feature_i
  feature_bar[,11]<--log10(feature_i$pvalue)
  #ggplot(feature_bar, aes(V11,Description))+
    #geom_bar(stat = "identity",fill="#ABCECD")+ #green
    #geom_bar(stat = "identity",fill="#C09C60")+ #brown
    #geom_bar(stat = "identity",fill="#6398AD")+ #blue
  #geom_bar(stat = "identity",fill="#AF4D47")+ #
    #coord_flip()+
  #theme_classic()+
  #scale_fill_manual()+
  #labs(y=feature[feature1[i]])
  
  
  for (j in 1:length(ontology)) {
    ontology_i<- feature_i[feature_i$ONTOLOGY==ontology[j],] 
    if(nrow(ontology_i)!=0){
      dataout_p[j,i]<-(-log10(mean(ontology_i$pvalue)))
      dataout_GR[j,i]<-mean(ontology_i$Count)/as.numeric( strsplit(as.character( ontology_i$GeneRatio[1]),"/")[[1]][2])
      
    }
    
  }
}

pheatmap(dataout_p[4,], 
         #color = colorRampPalette(c("#A9D3CF", "#5FA09B"))(200),#BP
         #color = colorRampPalette(c("#E8ACAC", "#C96D6D"))(200),#MF
         #color = colorRampPalette(c("#91D1E5", "#2A758C"))(200),#CC
         color = colorRampPalette(c("#EFBFA5", "#D97435"))(200),#KEGG
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

ggplot(dataout_p1, aes(x=value, y=name, fill=name))+
  geom_density_ridges(scale=4,alpha=0.8)+
  scale_fill_cyclical(values = c("#5FA09B", "#2A758C","#D97435","#C96D6D"))+
  theme_ridges(grid = FALSE)
