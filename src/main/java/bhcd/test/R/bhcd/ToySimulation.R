##############################################################################
# Import Toy Network and Export as Edgelist
#
# 
##############################################################################
library(igraph)
library(ggplot2)


setwd("~/Google_Drive/blang/workspace/blangBHCD")

file <- "n100seed2021.csv"

# save adjacency matrix as data frame
adj <- read.csv(paste("data/simulated/toy/", file, sep=""), header = FALSE)
# convert data frame to matrix
adj <- data.matrix(adj)
colnames(adj) <- NULL
# convert the matrix to an igraph object
graph <- graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)

# plot graph
plot(graph, vertex.size=5)

# Export the graph as an edgelist
write.graph(graph, file="data/simulated/toy/edgelist_n100seed2021alpha15beta5.txt", format="edgelist")



# edge count
ecount(graph) # 96 # 28 # 35 # 3786 # 330 # 3229
vcount(graph) # 20 # 10 # 10 # 100 # 100 # 100

# Network summary stats:

# create a dataframe to store stats
stats <- data.frame(edg_dens=edge_density(graph), cc=transitivity(graph, type="global"), 
											diam=diameter(graph, directed=F), avg_dist=mean_distance(graph, directed=F), 
											avg_deg=mean(degree(graph, mode="all")))


# plot the degree dist'n

# get the vertex degrees
deg <- data.frame(deg=degree(graph, mode='all'))

(p1 <- ggplot(deg, aes(x=deg)) + geom_histogram(aes(y=..density..), color="darkblue", fill="lightblue", binwidth = 1) +
		ggtitle("Histogram of Vertex Degree") + ylab("fraction of vertices") +
		xlab("degree"))

# culmulative degree distribution
deg.dist <- degree_distribution(graph, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")


#######################################
# Simulation from inferred params
#
#######################################


# dataframe of simulation parameters
size <- 62
alpha <- 1.29644374922349E-06 # posterior mean
beta <- 0.115393178016197 # posterior mean

# plot the beta prior
x <-  seq(0,1, length=100000)
plot(x, dbeta(x, alpha, beta), ylab="density", type ="l", col=2, main="Beta Prior with Posterior Mean Parameters")

# quick sim of internal node probs
unique(rbeta(10000, alpha, beta))
max(rbeta(10000, alpha, beta))
(exp <- alpha / (alpha + beta))# expected value
(size/2)^2 * exp # expected num of edges from even root split
1 - dbinom(0, size=961, prob=exp) # prob that num of edges is >= 1 from even root split

num_sims <- 10000
seed <- as.character(seq(2020, length.out=num_sims, by=1))

params <- data.frame(seed)

# add columns for clustering coefficient, avg distance, edge density, diameter, avg degree
params$cc <- NA
params$avg_dist <- NA
params$edg_dens <- NA
params$diam <- NA
params$avg_deg <- NA
params$ecount <- NA # edge count : check if graph is empty
params$max_pr <- NA
params$avg_pr <- NA


for(i in 1:nrow(params)) {
	
	file <- paste("n", size, "seed", params$seed[i], ".csv", sep="")
	file_p <- paste("probs_n", size, "seed", params$seed[i], ".csv", sep="")
	
	# save adjacency matrix as data frame
	adj <- read.csv(paste("data/simulated/terrorist/", file, sep=""), header = FALSE)
	# convert data frame to matrix
	adj <- data.matrix(adj)
	colnames(adj) <- NULL
	# convert the matrix to an igraph object
	graph <- graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)
	
	
	# Network summary stats:
	params$edg_dens[i] <- edge_density(graph) # proportion of present edges from all possible edges in the network
	
	
	params$cc[i] <- transitivity(graph, type="global") # ratio of triangles (direction disregarded) to connected triples
	# Transitivity measures the probability that the adjacent vertices of a vertex are connected. 
	# This is sometimes also called the clustering coefficient.
	
	params$diam[i] <- diameter(graph, directed=F) # the longest geodesic distance (length of the shortest path between two nodes) in the network
	
	params$avg_dist[i] <- mean_distance(graph, directed=F) # the mean of the shortest distance between each pair of nodes in the network
	
	params$avg_deg[i] <- mean(degree(graph, mode="all")) # mean vertex degree
	
	params$ecount[i] <- ecount(graph) # number of edges
	
	# save the internal node probabilties as data frame
	p <- read.csv(paste("data/simulated/terrorist/", file_p, sep=""), header = FALSE)
	
	params$max_pr <- max(p$V1) # maximum internal node prob.
	params$avg_pr <- mean(p$V1) # mean internal node prob.
	
	
}


which(params$ecount > 0)
length(which(params$ecount > 0)) # of non-empty networks
unique(params$max_pr)

mean(params$ecount) # vs. 152
max(params$ecount)

mean(params$edg_dens)
mean(params$diam)
mean(params$avg_deg)


#params$ab <- interaction(params$alpha, params$beta)
#params$ab <- factor(params$ab , levels=c("1.0.1.0", "0.1.0.1", "0.1.1.0", "1.0.0.1"))

# PLOTS
ggplot(params, aes(y=ecount)) + geom_boxplot()


p1 <- ggplot(params, aes(y=cc, x=ab, fill=ab)) + geom_boxplot() + ggtitle("Clustering Coefficient by Hyperparameters") +
	xlab("alpha, beta") + ylab("clustering coefficient") +
	scale_fill_discrete(name="Prior\nShape",
											breaks=c("1.0.1.0", "0.1.0.1", "0.1.1.0", "1.0.0.1"),
											labels=c("Uniform", "U-shape", "L-shape", "Reverse-L"))


p2 <- ggplot(params, aes(y=avg_dist, x=ab, fill=ab)) + geom_boxplot() + ggtitle("Average Vertex-Vertex Distance by Hyperparameters") +
	xlab("alpha, beta") + ylab("avg. distance") +
	scale_fill_discrete(name="Prior\nShape",
											breaks=c("1.0.1.0", "0.1.0.1", "0.1.1.0", "1.0.0.1"),
											labels=c("Uniform", "U-shape", "L-shape", "Reverse-L"))


p3 <- ggplot(params, aes(y=edg_dens , x=ab, fill=ab)) + geom_boxplot() + ggtitle("Edge Density by Hyperparameters") +
	xlab("alpha, beta") + ylab("edge density") +
	scale_fill_discrete(name="Prior\nShape",
											breaks=c("1.0.1.0", "0.1.0.1", "0.1.1.0", "1.0.0.1"),
											labels=c("Uniform", "U-shape", "L-shape", "Reverse-L"))


p4 <- ggplot(params, aes(y=avg_deg , x=ab, fill=ab)) + geom_boxplot() + ggtitle("Average Degree by Hyperparameters") +
	xlab("alpha, beta") + ylab("avg. degree") +
	scale_fill_discrete(name="Prior\nShape",
											breaks=c("1.0.1.0", "0.1.0.1", "0.1.1.0", "1.0.0.1"),
											labels=c("Uniform", "U-shape", "L-shape", "Reverse-L"))


p5 <- ggplot(params, aes(y=diam , x=ab, fill=ab)) + geom_boxplot() + ggtitle("Graph Diameter by Hyperparameters") +
	xlab("alpha, beta") + ylab("diameter") +
	scale_fill_discrete(name="Prior\nShape",
											breaks=c("1.0.1.0", "0.1.0.1", "0.1.1.0", "1.0.0.1"),
											labels=c("Uniform", "U-shape", "L-shape", "Reverse-L"))






# save plots to single pdf
pdf("saved_results/simulation/summary_stats.pdf", onefile = TRUE)
p1
p2
p3
p4
p5
dev.off()





# plot all degree distributions together
ggplot(deg_df, aes(x=deg, color=seed, fill=seed)) + geom_histogram(alpha=0.6, binwidth = 1) +
	ggtitle("Histogram of Vertex Degree: n=500, {alpha, beta} = {0.1, 1.0}; internal node probs. fixed") + xlab("degree") + facet_grid(seed~.) 

#+ scale_fill_viridis(discrete=TRUE) + scale_color_viridis(discrete=TRUE)

# plot all of the internal node prob histograms together
ggplot(probs_df, aes(x=p, color=seed, fill=seed)) + geom_histogram(alpha=0.6, bins=100, boundary=0) +
	ggtitle('Histogram of Internal Node Probabilities') + xlab('internal node probability') + facet_grid(seed~.)




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
y <- dbeta(x, shape1=alpha, shape2=beta)
density <- data.frame(x=x, y=y)
ggplot(p, aes(x=V1)) + geom_histogram(color='black', fill='lightgrey', bins=50, boundary=0) +
	geom_area(data=density, aes(x=x, y=dbeta(x, alpha, beta)), fill='#FF6666', alpha=.4) +
	ggtitle('Histogram of Internal Node Probabilities') + xlab('internal node probability')


# create beta distribution
#ggplot(density, aes(x, y)) + geom_area(fill="lightblue")


# save plots to single pdf
pdf("saved_results/simulation/plots.pdf", onefile = TRUE)
plot(graph, vertex.size=5, vertex.label=NA, main=title)
grid.arrange(p2, p1, ncol = 1)
dev.off()



######
# Val checking
alpha <- 0.1
beta <- 1
a <- 0.34
as <- 0.07
b <- 2.22
bs <- 0.91

a - 3*as
b - 1*bs

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
