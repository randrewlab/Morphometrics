clus <- cutree(clust, 4)
g <- split(names(clus), clus)

p <- ggtree(clust, linetype='dashed')
clades <- sapply(g, function(n) MRCA(p, n))

p <- groupClade(p, clades, group_name='subtree') + aes(color=subtree)

# d <- data.frame(label = names(clus), 
#                 cyl = datasetisla[names(clus), "cyl"])

## Dummy grouping 
d <- data.frame(label = names(clus), 
                pop = samplemeta$group1[match(names(clus),samplemeta$long.id)])

p %<+% d + 
  layout_dendrogram() + 
  geom_tippoint(aes(fill=factor(pop), x=x+.5),
                size=5, shape=21, color='black') +
  geom_tiplab(aes(label=pop), size=3, hjust=.5, color='black') +
  # geom_tiplab(angle=90, hjust=1, offset=-10, show.legend=FALSE) + 
  scale_color_brewer(palette='Set1', breaks=1:4) +
  theme_dendrogram(plot.margin=margin(6,6,80,6)) +
  theme(legend.position=c(.9, .6))
