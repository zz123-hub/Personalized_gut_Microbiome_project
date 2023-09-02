library(VennDiagram)
data= read.table("222.txt",header =TRUE,row=1,sep="\t")
option ="BT1-GT1-TT1"
ID =rownames(data)
sets =colnames(data)
newgroup =unlist(strsplit(option,"-"))
x =list()
for(i in newgroup){
  x[[i]]=ID[as.numeric(data[,i])>0]
}
color_list = c("#A9B88A","#B49592","#98B3C0")
fill_colors = color_list[1:length(x)]
P <- venn.diagram(x,col=fill_colors,fill=fill_colors,lwd=3,filename=NULL,fontfamily = "sans",alpha=0,cex=2,
                  cat.cex=2,cat.fontfamily = "sans",width=1200,height=1200)
pdf("venn-pdf.pdf")
grid.draw(P)
dev.off()
#"#7C8A6C","#95514F","#7A8DA6"
#https://www.jianshu.com/p/f858521828a5
#https://zhuanlan.zhihu.com/p/370916031