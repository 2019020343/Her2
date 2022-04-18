rm(list = ls())
#install.packages("igraph")
library(igraph)
posgene<-read.table("D:\\HER2\\gene\\result\\pos\\val\\posgene.txt",header=T,sep = "\t", quote = "\"'")
GPL<-read.table("D:\\HER2\\gene\\result\\pos\\GSE45827\\GPL570-55999.txt",header=T,sep = "\t", quote = "\"'",fill = T)
GPL1<-GPL[,c(11,12)]
myfun1<-function(x){unlist(strsplit(as.character(x),"-"))[1]}
feature<-unique(unlist(lapply(posgene[,1],myfun1)))
myfun2<-function(x){unlist(strsplit(as.character(x),"/"))}
fname<-c("original_glszm_SizeZoneNonUniformityNormalized",
         "original_glszm_ZoneEntropy",
         "original_glrlm_ShortRunHighGrayLevelEmphasis",
         "original_glszm_SizeZoneNonUniformity",
         "original_glszm_SmallAreaEmphasis",
         "original_glrlm_ShortRunEmphasis",
         "original_glrlm_RunLengthNonUniformityNormalized",
         "original_glszm_SmallAreaHighGrayLevelEmphasis")
feature1<-unique(unlist(lapply(fname,function(x){which(feature==x)})))

dataout1<-c()
edges<-c()
nodes<-c()
nodes_gene_out<-c()
nodes_function_out<-c()
data1<-c()
data2<-c()
edges_gene_out<-c()
for(i in 1:8){
  feature_i<-read.table(paste("D:\\her2\\gene\\result\\pos\\val\\feature_", feature1[i], "_pos.txt", sep = ""),header=T,sep = "\t", quote = "\"'")

  if(nrow(feature_i)!=0){
    for (j in 1:nrow(feature_i)) {
      if(feature_i[j,1]!="kegg"){
        gene<-unique(unlist(lapply(feature_i$geneID[j],myfun2)))
        edges_gene_out<-unique(c(edges_gene_out,gene))
        G_G<-cbind(gene,as.character(feature_i$ID[j]), -log10(feature_i$p.adjust[j]))
        data1<-c(data1,gene)
        data2<-unique(rbind(data2,G_G)) 
      }else{
        gene<-unique(unlist(lapply(feature_i$geneID[j],myfun2)))
        index<-match(gene,GPL1[,2])
        gene<-c(as.character( GPL1[index[which(is.na(index)==0)],1]),gene[which(is.na(index)==1)] )
        edges_gene_out<-unique(c(edges_gene_out,gene))
        G_G<-cbind(gene,as.character(feature_i$ID[j]), -log10(feature_i$p.adjust[j]))
        data1<-c(data1,gene)
        data2<-unique(rbind(data2,G_G)) 
     }
    }
    
    dataout<-cbind(feature[feature1[i]], as.data.frame(table(data1)))
    colnames(dataout)<-c("from","to","fre")
    dataout1<-rbind(dataout1,dataout)
    colnames(dataout1)<-c("from","to","fre")
    colnames(data2)<-c("from","to","fre")
    edges<-rbind(dataout1,data2)
    
    nodes_gene_out<-unique(c(nodes_gene_out,edges_gene_out))
    nodes_function_out<-unique(c(nodes_function_out,data2[,2]))
    
    
  }  
}
nodes<-rbind(cbind(nodes_gene_out,"Genes"),cbind(fname,"Features"),cbind(nodes_function_out,"Functions"))

a<-unique(data2[,1])
edges<-as.data.frame(edges)
nodes<-as.data.frame(nodes)
write.table(edges,"D:\\HER2\\fig\\4\\edges.txt", sep = "\t", quote = F, row.names = F)
write.table(nodes,"D:\\HER2\\fig\\4\\nodes.txt", sep = "\t", quote = F, row.names = F)
mer<-merge(dataout1,data2,by.x = "to",by.y = "from")
mer1<-unique(mer[,c(2,4)])
#2)导入数据后，要转化成图数据才能用R作图，不同数据格式用不同方式
#graph_from_literal 通过文字创建，graph_from_edgelist通过边列表（矩阵）创建，graph_from_adj_list通过邻接列表（矩阵）创建
#graph_from_adjacency_matrix 通过邻接矩阵（所有点的纵横矩阵）创建，graph_from_incidence_matrix通过关联矩阵（两组内部独立点的矩阵）创建
#graph_from_data_frame通过数据框创建，详细介绍可参见《igraph manual》http://igraph.org/r/doc/igraph.pdf
#这里数据格式是数据框，所以用graph_from_data_frame
graph <- graph_from_data_frame(edges, directed = TRUE) #directed = TRUE表示有方向,如果不需要点数据，可以设置vertices=NULL

#3)转换完成后，有两种生成方式，一是直接plot,参数放里面；二是通过修改图的方式设置参数，然后plot
source("http://michael.hahsler.net/SMU/ScientificCompR/code/map.R")

#生成方式1：
plot(graph,  
     layout=layout.reingold.tilford,  #layout.fruchterman.reingold表示弹簧式发散的布局，
     #其他还有环形布局layout.circle，分层布局layout.reingold.tilford，中心向外发散layout.reingold.tilford(graph,circular=T) ，核心布局layout_as_star，大型网络可视化layout_with_drl
     vertex.size=5,     #节点大小  
     vertex.shape='circle',    #节点不带边框none,,圆形边框circle,方块形rectangle  
     vertex.color=map(degree(graph), c(1,20)),#设置颜色，其他如red,blue,cyan,yellow等
     vertex.label=NULL, #NULL表示不设置，为默认状态，即默认显示数据中点的名称，可以是中文。如果是NA则表示不显示任何点信息	 
     vertex.label.cex=0.1,    #节点字体大小  
     vertex.label.color='black',  #节点字体颜色,red  
     vertex.label.dist=0.2,   #标签和节点位置错开
     edge.arrow.size=0.01,#连线的箭头的大小,若为0即为无向图，当然有些数据格式不支持有向图  
     edge.width = 0.1, #连接线宽度
     edge.label=NA, #不显示连接线标签，默认为频次
     edge.color="gray50")  #连线颜色 
#生成方式2set.seed(50) #生成随机数，这样图的布局就会可重复，而不是每次生成的时候都变
l<-layout.circle(graph) #设置图的布局方式为弹簧式发散的布局
#具体修改过程
V(graph)$size <- degree(graph)  #节点大小与点中心度成正比，中心度即与该点相连的点的总数
colrs <- c("gray50", "tomato", "orange", "red", "yellow")
V(graph)$color <- colrs[V(graph)$V2] #根据类型设置颜色,按照类型分组
V(graph)$label.color <- 'white' #设置节点标记的颜色
E(graph)$width <- E(graph)$fre #根据频次列设置边宽度
E(graph)$label <- E(graph)$fre #根据频次列设置边标签
E(graph)$arrow.size=0.03 #设置箭头大小
#生成图
plot(graph, layout=l)

