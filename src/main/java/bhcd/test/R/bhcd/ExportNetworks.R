##############################################################################
# Export real networks to CSV files
#
# 
##############################################################################
library(igraph)
library(igraphdata)

setwd("~/Google_Drive/blang/workspace/blangBHCD")

data(package="igraphdata")

data("yeast")
vcount(yeast)
ecount(yeast)
is_directed(yeast)
is_weighted(yeast)
edge_density(yeast)

write.graph(yeast, file="data/yeast/yeast.txt", format="edgelist")


data("karate")
vcount(karate)
ecount(karate)
is_directed(karate)
is_weighted(karate)
plot(karate)
V(karate)
edge_density(karate)

as_edgelist(karate)
write.graph(karate, file="data/karate/karate.txt", format="edgelist")
