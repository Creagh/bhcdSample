##############################################################################
# Read in a simulated network stored as an adjacency matrix in CSV format
#
# See TestForwardSampling.java for the simulation
##############################################################################
library(igraph)
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
#library(hrbrthemes)

setwd("~/Google_Drive/blang/workspace/blangBHCD")


# dataframe of simulation parameters
size <- rep(500, 9)
seed <- as.character(seq(from=2020, to=2028, by=1))
alpha <- rep("0.1", 9)
beta <- rep("1.0", 9)
params <- data.frame(size, seed, alpha, beta)

# dataframe of vertex degrees
deg_df <- data.frame(seed=integer(), deg=integer())
# dataframe of internal node probs
probs_df <- data.frame(seed=integer(), p=double())

for(i in 1:nrow(params)) {
	
	file <- paste("n", params$size[i], "seed", params$seed[i], "alpha", params$alpha[i], "beta", params$beta[i], ".csv",
											sep="")
	file_p <- paste("probs_n", params$size[i], "seed", params$seed[i], "alpha", params$alpha[i], "beta", params$beta[i], ".csv",
									sep="")
	
	# save adjacency matrix as data frame
	adj <- read.csv(paste("data/simulated/", file, sep=""), header = FALSE)
	# convert data frame to matrix
	adj <- data.matrix(adj)
	colnames(adj) <- NULL
	# convert the matrix to an igraph object
	graph <- graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)
	
	# get the vertex degrees
	deg <- degree(graph, mode='all')
	# append to data frame
	deg_df <- rbind(deg_df, data.frame(seed=params$seed[i], deg))
	
	# save the internal node probabilties as data frame
	p <- read.csv(paste("data/simulated/", file_p, sep=""), header = FALSE)
	# append to the data frame
	probs_df <- rbind(probs_df, data.frame(seed=params$seed[i], p=p$V1))
	
}

# plot all degree distributions together
ggplot(deg_df, aes(x=deg, color=seed, fill=seed)) + geom_histogram(alpha=0.6, binwidth = 1) +
	ggtitle("Histogram of Vertex Degree: n=500, {alpha, beta} = {0.1, 1.0}; internal node probs. fixed") + xlab("degree") + facet_grid(seed~.) 

#+ scale_fill_viridis(discrete=TRUE) + scale_color_viridis(discrete=TRUE)

# plot all of the internal node prob histograms together
ggplot(probs_df, aes(x=p, color=seed, fill=seed)) + geom_histogram(alpha=0.6, bins=100, boundary=0) +
	ggtitle('Histogram of Internal Node Probabilities') + xlab('internal node probability') + facet_grid(seed~.)


# Network summary stats:

edge_density(graph) # proportion of present edges from all possible edges in the network

transitivity(graph, type="global") # ratio of triangles (direction disregarded) to connected triples
# Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
# This is sometimes also called the clustering coefficient.

diameter(graph, directed=F) # the longest geodesic distance (length of the shortest path between two nodes) in the network

mean_distance(graph, directed=F) # the mean of the shortest distance between each pair of nodes in the network

mean(degree(graph, mode="all")) # mean vertex degree

# string of simulation parameters
#a <- params$alpha[i] # alpha
#b <- params$beta[i] # beta
#title <- expression(paste("n=200; (", alpha, ", ", beta, ")=(1.0, 0.1); seed 2020"))


# plot the graph using a force-directed layout
#plot(graph, vertex.size=5, vertex.label=NA, main=title)


# histogram of the degrees
#hist(deg, breaks=1:vcount(graph)-1, main="Histogram of node degree")


p1 <- ggplot(deg_df, aes(x=deg)) + geom_histogram(aes(y=..density..), color="darkblue", fill="lightblue", binwidth = 1) +
	ggtitle("Histogram of Vertex Degree") + geom_density(alpha=.2, fill="#FF6666") + ylab("fraction of vertices") +
	xlab("degree")

# culmulative degree distribution
deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
#plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")


# save the internal node probabilties as data frame
#p <- read.csv(paste("data/simulated/", filep3, sep=""), header = FALSE)

# histogram of the internal node probabilties
#ggplot(p, aes(x=V1)) + geom_histogram(color='darkred', fill='pink', bins=50) +
#	ggtitle('Histogram of Internal Node Probabilities') + xlab('internal node probability')

# histogram with overlay of Beta distribution
x <- seq(0, 1, 0.01)
y <- dbeta(x, shape1=a, shape2=b)
density <- data.frame(x=x, y=y)
p2 <- ggplot(p, aes(x=V1)) + geom_histogram(color='black', fill='lightgrey', bins=50, boundary=0) +
	geom_area(data=density, aes(x=x, y=dbeta(x, a, b)), fill='#FF6666', alpha=.4) +
	ggtitle('Histogram of Internal Node Probabilities') + xlab('internal node probability')


# create beta distribution
#ggplot(density, aes(x, y)) + geom_area(fill="lightblue")


# save plots to single pdf
pdf("saved_results/simulation/plots.pdf", onefile = TRUE)
plot(graph, vertex.size=5, vertex.label=NA, main=title)
grid.arrange(p2, p1, ncol = 1)
dev.off()





# Community detection based on edge betweenness (Newman-Girvan)
# High-betweenness edges are removed sequentially (recalculating at each step) 
# and the best partitioning of the network is selected.
#
# Note: The edge betweenness centrality is defined as the number of the shortest paths 
# that go through an edge in a graph or network (Girvan and Newman 2002). Each edge in the 
# network can be associated with an edge betweenness centrality value. An edge with a high 
# edge betweenness centrality score represents a bridge-like connector between two parts of a 
# network, and the removal of which may affect the communication between many pairs of nodes 
# through the shortest paths between them.
ceb <- cluster_edge_betweenness(graph) 
dendPlot(ceb, mode="hclust")
