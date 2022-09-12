setwd("") 
install.packages("Metrics",repos='http://cran.us.r-project.org')

library(glmnet)
library(Metrics)
library(cmapR)

col_meta <- read_gctx_meta("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx", dim="col")
row_meta <- read_gctx_meta("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx", dim="row")
info=read.table("GSE70138_Broad_LINCS_inst_info.txt",sep = '\t',header = T,quote="")

cellname=c("MCF7")  
cellname=c("HCC515")               
cellname=c("A549","A549.11")     

#index=which(info$cell_id %in% "MCF7")
index=which(info$cell_id %in% cellname)
brest_info=info[index,]
insname=info[index,]$inst_id
group1=sapply(strsplit(insname,"\\_"),"[",-4)
group=group1[4,]
group=sapply(strsplit(group,"\\:"),"[",-1)

insgctname=matrix(data=NA, nrow=ncol(group1))
for (i in 1:ncol(group1)){
  insgctname1=paste(group1[1,i],group1[2,i],group1[3,i],sep="_")
  insgctname[i]=paste(insgctname1,group[i],sep = ":")
}

head(insgctname)
brestinfonew=cbind(insgctname,brest_info[,2:12])
write.table(brestinfonew,file="brest.MCF_info.txt",sep='\t',row.names = T)
write.table(brestinfonew,file="liver.HCC_info.txt",sep='\t',row.names = T)
write.table(brestinfonew,file="lung.A549_info.txt",sep='\t',row.names = T)

brest_ds=parse_gctx("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328.gctx", cid = unique(insgctname))


