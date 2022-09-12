setwd("")
ctdpath=read.table("KEGG PATHWAY GENE.txt",sep = '\t',header = T,quote="")


load("brest.MCF.rdata")
load("liver.HCC.rdata")
load("lung.A549.rdata")

brest_result=read.csv("lung score.csv", sep=",")
head(brest_result)
brest_order=brest_result[order(-brest_result$Scores),]
top_gene=brest_order$Gene[1:500]
cluster=read.table("cluster.txt",sep="\t",header=T, check.names = F)
dim(cluster)
genename=intersect(cluster$Gene,rownames(brest_ds_data))
cluster=cluster[cluster$Gene %in% genename,]
dim(cluster)
path_gene=genename

toppath_gene=top_gene[top_gene %in% path_gene]
topnew_gene=setdiff(top_gene,path_gene)

cluster1=cluster$Gene[cluster$Cluster==1]
cl1=top_gene[top_gene %in% cluster1]
cluster2=cluster$Gene[cluster$Cluster==2]
cl2=top_gene[top_gene %in% cluster2]
cluster3=cluster$Gene[cluster$Cluster==3]
cl3=top_gene[top_gene %in% cluster3]
cluster4=cluster$Gene[cluster$Cluster==4]
cl4=top_gene[top_gene %in% cluster4]
cluster5=cluster$Gene[cluster$Cluster==5]
cl5=top_gene[top_gene %in% cluster5]
cluster1=cbind(cl1,rep("1",length(cl1)))
cluster2=cbind(cl2,rep("2",length(cl2)))
cluster3=cbind(cl3,rep("3",length(cl3)))
cluster4=cbind(cl4,rep("4",length(cl4)))
cluster5=cbind(cl5,rep("5",length(cl5)))

topcluster=rbind(cluster1,cluster2,cluster3,cluster4,cluster5)

write.table(topcluster,file="liver_top500cluster.txt",sep="\t",row.names = F)
write.table(topcluster,file="brest_top1000cluster.txt",sep="\t",row.names = F)
write.table(topcluster,file="lung_top500cluster.txt",sep="\t",row.names = F)


####diff stat####

fdrFilter=0.05                                                            #fdr¡ŸΩÁ÷µ
logFCfilter=1      
pvaluefilter=0.05

cluGroup=levels(factor(cluster$Cluster))

statNames=c("id")
diffStat=data.frame(id=row.names(rt))
for(group in cluGroup){
  sample1=as.character(cluster$Gene[cluster$Cluster==group])
  sample2=as.character(cluster$Gene[cluster$Cluster!=group])
  conNum=length(sample2)
  treatNum=length(sample1)
  grade=c(rep(1,conNum),rep(2,treatNum))
  data=cbind(rt[,sample2],rt[,sample1])
  
  outTab=data.frame()
  
  #diff analysis
  for(i in row.names(data)){
    rt1=rbind(expression=data[i,],grade=grade)
    rt1=as.matrix(t(rt1))
    wilcoxTest<-wilcox.test(expression ~ grade, data=rt1)
    # conGeneMeans=mean(data[i,1:conNum])
    # treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
    conGeneMeans=abs(mean(data[i,1:conNum]))
    treatGeneMeans=abs(mean(data[i,(conNum+1):ncol(data)]))
    logFC=log2(treatGeneMeans)-log2(conGeneMeans)
    pvalue=wilcoxTest$p.value
    conMed=median(data[i,1:conNum])
    treatMed=median(data[i,(conNum+1):ncol(data)])
    diffMed=treatMed-conMed
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
  }
  
  
  pValue=outTab[,"pValue"]
  fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
  outTab=cbind(outTab,fdr=fdr)
  
  write.table(outTab,file=paste0("cluster",group,".all.xls"),sep="\t",row.names=F,quote=F)
  
  outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
  write.table(outDiff,file=paste0("cluster",group,".diff.xls"),sep="\t",row.names=F,quote=F)
  
  diffExp=rt[as.vector(outDiff[,1]),]
  dim(diffExp)
  diffExp=cbind(id=row.names(diffExp),diffExp)
  write.table(diffExp,file=paste0("cluster",group,".diffExp.txt"),sep="\t",row.names=F,quote=F)
  
  # stat=ifelse(abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter,1,0)
  stat=ifelse(as.numeric(as.vector(outTab$pValue))<pvaluefilter,1,0)
  diffStat=cbind(diffStat,stat)
  statNames=c(statNames,paste0("C",group))
}

colnames(diffStat)=statNames
write.table(diffStat,file="diffStat.txt",sep="\t",row.names=F,quote=F)




