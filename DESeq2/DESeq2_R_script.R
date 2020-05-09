library(DESeq2)

count<-read.csv("readcounts_data.txt",sep="\t",row.names=1) # read counts file
count<-as.matrix(count)

header<-read.csv("DESeq2_header.txt",header=F)   # header file
header=as.vector(as.character(header[,1]))

group <- data.frame(con = factor(header))

dds <- DESeqDataSetFromMatrix(countData = count, colData = group, design = ~ con)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

#save(dds,file="dds.rda")


### Compare samples at 10:00 ###
res <- results(dds, contrast = c("con", "Bd21.10", "col.10"))
write.table(res,"Bd_col.10.txt")
res <- results(dds, contrast = c("con", "Bd21.10", "sta.10"))
write.table(res,"Bd_Bs.10.txt")
res <- results(dds, contrast = c("con", "Bd21.10", "Bh.10"))
write.table(res,"Bd_Bh.10.txt")
res <- results(dds, contrast = c("con", "col.10", "sta.10"))
write.table(res,"col_Bs.10.txt")
res <- results(dds, contrast = c("con", "col.10", "Bh.10"))
write.table(res,"col_Bh.10.txt")
res <- results(dds, contrast = c("con", "sta.10", "Bh.10"))
write.table(res,"Bs_Bh.10.txt")
res <- results(dds, contrast = c("con", "BhBd.10", "BhBs.10"))
write.table(res,"BhBd_BhBs.10.txt")



### Compare samples at 14:00 ###
res <- results(dds, contrast = c("con", "Bd21.14", "col.14"))
write.table(res,"Bd_col.14.txt")
res <- results(dds, contrast = c("con", "Bd21.14", "sta.14"))
write.table(res,"Bd_Bs.14.txt")
res <- results(dds, contrast = c("con", "Bd21.14", "Bh.14"))
write.table(res,"Bd_Bh.14.txt")
res <- results(dds, contrast = c("con", "col.14", "sta.14"))
write.table(res,"col_Bs.14.txt")
res <- results(dds, contrast = c("con", "col.14", "Bh.14"))
write.table(res,"col_Bh.14.txt")
res <- results(dds, contrast = c("con", "sta.14", "Bh.14"))
write.table(res,"Bs_Bh.14.txt")
res <- results(dds, contrast = c("con", "BhBd.14", "BhBs.14"))
write.table(res,"BhBd_BhBs.14.txt")



### Compare samples at 18:00 ###
res <- results(dds, contrast = c("con", "Bd21.18", "col.18"))
write.table(res,"Bd_col.18.txt")
res <- results(dds, contrast = c("con", "Bd21.18", "sta.18"))
write.table(res,"Bd_Bs.18.txt")
res <- results(dds, contrast = c("con", "Bd21.18", "Bh.18"))
write.table(res,"Bd_Bh.18.txt")
res <- results(dds, contrast = c("con", "col.18", "sta.18"))
write.table(res,"col_Bs.18.txt")
res <- results(dds, contrast = c("con", "col.18", "Bh.18"))
write.table(res,"col_Bh.18.txt")
res <- results(dds, contrast = c("con", "sta.18", "Bh.18"))
write.table(res,"Bs_Bh.18.txt")
res <- results(dds, contrast = c("con", "BhBd.18", "BhBs.18"))
write.table(res,"BhBd_BhBs.18.txt")



### Compare samples at 22:00 ###
res <- results(dds, contrast = c("con", "Bd21.22", "col.22"))
write.table(res,"Bd_col.22.txt")
res <- results(dds, contrast = c("con", "Bd21.22", "sta.22"))
write.table(res,"Bd_Bs.22.txt")
res <- results(dds, contrast = c("con", "Bd21.22", "Bh.22"))
write.table(res,"Bd_Bh.22.txt")
res <- results(dds, contrast = c("con", "col.22", "sta.22"))
write.table(res,"col_Bs.22.txt")
res <- results(dds, contrast = c("con", "col.22", "Bh.22"))
write.table(res,"col_Bh.22.txt")
res <- results(dds, contrast = c("con", "sta.22", "Bh.22"))
write.table(res,"Bs_Bh.22.txt")
res <- results(dds, contrast = c("con", "BhBd.22", "BhBs.22"))
write.table(res,"BhBd_BhBs.22.txt")



### Compare samples at 2:00 ###
res <- results(dds, contrast = c("con", "Bd21.2", "col.2"))
write.table(res,"Bd_col.2.txt")
res <- results(dds, contrast = c("con", "Bd21.2", "sta.2"))
write.table(res,"Bd_Bs.2.txt")
res <- results(dds, contrast = c("con", "Bd21.2", "Bh.2"))
write.table(res,"Bd_Bh.2.txt")
res <- results(dds, contrast = c("con", "col.2", "sta.2"))
write.table(res,"col_Bs.2.txt")
res <- results(dds, contrast = c("con", "col.2", "Bh.2"))
write.table(res,"col_Bh.2.txt")
res <- results(dds, contrast = c("con", "sta.2", "Bh.2"))
write.table(res,"Bs_Bh.2.txt")
res <- results(dds, contrast = c("con", "BhBd.2", "BhBs.2"))
write.table(res,"BhBd_BhBs.2.txt")



### Compare samples at 6:00 ###
res <- results(dds, contrast = c("con", "Bd21.6", "col.6"))
write.table(res,"Bd_col.6.txt")
res <- results(dds, contrast = c("con", "Bd21.6", "sta.6"))
write.table(res,"Bd_Bs.6.txt")
res <- results(dds, contrast = c("con", "Bd21.6", "Bh.6"))
write.table(res,"Bd_Bh.6.txt")
res <- results(dds, contrast = c("con", "col.6", "sta.6"))
write.table(res,"col_Bs.6.txt")
res <- results(dds, contrast = c("con", "col.6", "Bh.6"))
write.table(res,"col_Bh.6.txt")
res <- results(dds, contrast = c("con", "sta.6", "Bh.6"))
write.table(res,"Bs_Bh.6.txt")
res <- results(dds, contrast = c("con", "BhBd.6", "BhBs.6"))
write.table(res,"BhBd_BhBs.6.txt")



