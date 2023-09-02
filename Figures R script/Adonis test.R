all <- read.table("distance.txt",header=T,row.names=1,sep="\t")
a <- read.table("group.txt",header=T,row.names=1,sep="\t")
library(permute)
library(lattice)
library(vegan)
adonis2(all~ group,data=a,permutations = 9999)

