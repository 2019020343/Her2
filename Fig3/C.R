setwd("D:\\HER2\\gene\\result")
posgene<-read.table("D:\\HER2\\gene\\result\\pos\\val\\posgene.txt",header=T,sep = "\t", quote = "\"'")
myfun1<-function(x){unlist(strsplit(as.character(x),"-"))[1]}
feature<-unique(unlist(lapply(posgene[,1],myfun1)))
feature_CF_out<-c()
feature_CF_out_2<-as.data.frame(matrix( nrow = 86, ncol = 86))
colnames(feature_CF_out_2)<-feature
rownames(feature_CF_out_2)<-feature
for(i in 1:86){
  feature_i<-read.table(paste("D:\\HER2\\gene\\result\\feature_", i, "_pos.txt", sep = ""),header=T,sep = "\t", quote = "\"'")
  for (j in 1:86) {
    feature_j<-read.table(paste("D:\\HER2\\gene\\result\\feature_", j, "_pos.txt", sep = ""),header=T,sep = "\t", quote = "\"'")
    sim<-length(intersect(feature_i$ID, feature_j$ID) )/min(length(feature_i$ID), length(feature_j$ID))
    feature_CF_out1<-cbind(feature[i],feature[j],sim)
    feature_CF_out<-rbind(feature_CF_out,feature_CF_out1)
    feature_CF_out_2[i,j]<-sim
    }
}
write.table(feature_CF_out,"D:\\HER2\\gene\\result\\sim.txt", sep = "\t", quote = F, row.names = T)
write.table(feature_CF_out_2,"D:\\HER2\\gene\\result\\sim2.txt", sep = "\t", quote = F, row.names = T)

feature_CF_out<-read.table("D:\\HER2\\gene\\result\\sim.txt",header=T,sep = "\t", quote = "\"'")

library(ggplot2)
dftmp<-as.data.frame(feature_CF_out)
colnames(dftmp)<-c("x","y","value")
ggplot()+
  geom_tile(data=dftmp,aes(x,y),fill="white",color="grey")+
  geom_point(data=dftmp,aes(x,y,color=value),shape=15)+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  scale_x_discrete(position="top")+
  scale_y_discrete(labels =feature)+
  labs(x=NULL,y=NULL)+
  scale_colour_gradient(low = "white",high = "#C96D6D")

