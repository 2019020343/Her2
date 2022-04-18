rm(list=ls())
library(gplots)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(ggplot2)
pos_index1_2 <- read.delim('D:\\HER2\\gene\\result\\11\\pos_index1_2.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
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
      p_index<-as.matrix(gene_j[,seq(3,151,2)]) 
      R_index<-as.matrix(gene_j[,seq(2,151,2)])
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

write.table(FG_pr,"G:\\HER2\\fig\\2\\FG_pr.txt", sep = "\t", quote = F, row.names = F)

p1<-ggplot(FG_pr, aes(V1,V2,col=as.numeric(as.matrix(R_mean)) ))+
  geom_tile(color="grey85", fill="white", size=1)+
  geom_point(aes(size = as.numeric(as.matrix(P_mean))))+
  scale_colour_gradient2(high = "#F27D58",mid = "blue",low = "blue")+coord_flip()+
  theme(axis.text.x=element_text(angle = 45))
V2  

gname_number<-matrix(gname_number,ncol = 2,nrow = 8,byrow = T)
gname_number<-as.data.frame(gname_number)



p2<-ggplot(gname_number, aes(V1,V2))+
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






