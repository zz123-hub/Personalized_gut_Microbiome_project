all <- read.table("abundance.txt",header=T,row.names=1,sep="\t")
g<-factor(rep(1:3,c(6,6,6)),labels=c("BT1","GT1","TT1"))
result<-data.frame(as.character(rownames(all)))
for (i in 1:nrow(all)){
  result[i,2] <- kruskal.test(as.numeric(all[i,]), g)$p.value
}
write.table(result,"p_value.txt",sep="\t",row.names=F,col.names=F)
