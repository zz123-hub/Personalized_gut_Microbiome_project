library(reshape2)
library(ggplot2) 
df <- read.table("222.txt",header=T, sep="\t")
p<- ggplot(df, aes(x = group, y = genus, size = value2, color=value2)) + geom_point()+ theme_bw()+ 
  scale_colour_gradientn(colors=c("#5E789E","#DDDDDD","#95514F"))
p+x11(width=10.5,height=8.5)