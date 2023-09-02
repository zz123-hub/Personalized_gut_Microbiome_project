library(treemapify)
library(ggplot2)
df<- read.csv("33.csv", check.names = F)
p <-ggplot(df, aes(area = gdp_mil_usd, fill = country,
                   label = region,
                   subgroup = country)) +
  geom_treemap() +
  geom_treemap_subgroup_border() +
  geom_treemap_subgroup_text(place = "centre", alpha = 0.2, 
                             colour = "black", fontface = "italic",start="topleft",grow=TRUE) +
  geom_treemap_text(colour = "white", place = "topleft", reflow = T,min.size="10")+
  scale_fill_manual(values = c("#7C8A6C","#7A8DA6","#95514F"))
p+x11(width=6.5,height=5.5)