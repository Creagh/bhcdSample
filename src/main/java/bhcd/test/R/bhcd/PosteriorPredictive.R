##############################################################################
# Posterior Predictive Checks
#
# 
##############################################################################
library(igraph)
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
library(dplyr)
#library(hrbrthemes)

setwd("~/Google_Drive/blang/workspace/blangBHCD")

# Read in Terrorist data
tadj <- read.csv("data/terrorists/terrorist_pairs.csv", header = FALSE)
tadj <- read.csv("data/grass_web/grass_web_reduced.csv", header = FALSE)

# convert data frame to matrix
tadj <- data.matrix(tadj)
colnames(tadj) <- NULL
# convert the matrix to an igraph object

# FOR GRASS WEB
tadj <- tadj+1 # increase all vertex ID's by 1, so they don't start at zero

tgraph <- graph_from_edgelist(tadj, directed=FALSE)

# remove duplicate edges
tgraph <- simplify(tgraph, remove.multiple = TRUE, remove.loops = TRUE)
is_simple(tgraph)

# plot the graph using a force-directed layout
plot(tgraph, vertex.size=10, main="Terrorist Network")

# edge count
ecount(tgraph) # 152 t # 113 g
vcount(tgraph) # 62 t # 75 g

# Network summary stats:

# create a dataframe to store stats
stats <- data.frame(edg_dens=edge_density(tgraph), cc=transitivity(tgraph, type="global"), 
											diam=diameter(tgraph, directed=F), mean_dist=mean_distance(tgraph, directed=F), 
											mean_deg=mean(degree(tgraph, mode="all")))


##############################################################################
# Read in monitored statistics
##############################################################################
setwd("~/Desktop/computeCanada")

mon <- read.csv("terrorist/2020-10-31-01-38-53-f00p1BP8.exec/samples/monitor.csv", header = TRUE)

mon <- read.csv("grassWeb/2020-11-06-15-11-33-pYSkxguw.exec/samples/monitor.csv", header = TRUE)
mon <- read.csv("grassWeb/2020-11-06-15-23-14-8FqK7S4a.exec/samples/monitor.csv", header = TRUE)

thining <- 10
iters <- seq(0,max(mon$sample), by=thining)
maxV <- 75 # max vertex number
sim_stats <- data.frame(iters)

# Read in the forward-sampled edge lists
for(i in 1:length(iters)) {

	file <- paste("grassWeb/2020-11-06-15-23-14-8FqK7S4a.exec/sampled-edges/", iters[i]/thining, ".csv", sep='')
	edges <- tryCatch( 
		read.csv(file, header = FALSE),
		error = function(e) { # error occurs when the edge list is empty (no edges)
			return(data.frame(V1=integer(), V2=integer())) # create an empty data frame
		})
	
	edges <- data.matrix(edges)
	colnames(edges) <- NULL
	dframe <- as.data.frame(edges)
	graph <- graph.data.frame(dframe, vertices = 1:maxV , directed=F)
	#graph <- graph_from_edgelist(edges, directed=FALSE)
	
	# Calculate the summary stats
	sim_stats$edg_dens[i] <- edge_density(graph) 
	sim_stats$mean_deg[i] <- mean(degree(graph, mode="all"))
	sim_stats$cc[i] <- transitivity(graph, type="global")
	sim_stats$diam[i] <- diameter(graph, directed=F)
	sim_stats$mean_dist[i] <- mean_distance(graph, directed=F)
}


# OPTIONAL - Drop first x-values from burn-in
burn <- 25000
sim_stats <- tail(sim_stats, -burn/thining)

# Create new data frame for plotting purposes
df <- sim_stats
df$iters <- NULL
df <- stack(df)
df$real <- ifelse(df$ind == "edg_dens", stats$edg_dens,
									ifelse(df$ind == "mean_deg",stats$mean_deg,
										ifelse(df$ind == "cc",stats$cc,
											ifelse(df$ind == "diam", stats$diam,
												stats$mean_dist))))
df$indicator <- ifelse(df$values >= df$real, 1, 0) # indicator that test stat >= real value

# Function to calculate two more percentiles
add.percentiles = function(x) {
	r = quantile(x, probs=c(0.025, 0.975))
	names(r) = c("lwr", "upr")
	r
}

ggplot(df, aes(x = factor(ind, levels = names(sim_stats)), y = values)) + geom_boxplot() + 
	theme(legend.position = "none") + facet_wrap(~ind, scales="free") + xlab("") +
	geom_point(aes(y=real, x=ind, colour='red')) +
	stat_summary(fun.y=add.percentiles, geom="point", pch="_", colour="grey40", size=6) +
	ggtitle("Reference Distributions of Posterior Predictive Network Statistics - Zero Inflated Model")


# Posterior predictive p-values from Eqn 5 of Gelman 
# "POSTERIOR PREDICTIVE ASSESSMENT OF MODEL FITNESS VIA REALIZED DISCREPANCIES" (see also pg 741 Sec. 2.3)
# -- tail area probability of test stat under its reference distribution
mean(subset(df, ind == "edg_dens")$indicator)
mean(subset(df, ind == "mean_deg")$indicator)
mean(subset(df, ind == "cc")$indicator)
mean(subset(df, ind == "diam")$indicator)
mean(subset(df, ind == "mean_dist")$indicator)


# Combine boxplots from ZI and standard models
dfS <- df
dfS$model <- "Standard"
dfZI <- df
dfZI$model <- "ZI"

dfC <- rbind(dfS, dfZI)
levels(dfC$ind) <- c('Edge Density','Mean Degree','Clustering Coefficient', 'Diameter', 'Mean Vertex-Vertex Distance')

ggplot(dfC, aes(x = factor(model), y = values, fill=model)) + geom_boxplot() + 
	facet_wrap(~ind, scales="free") + theme(legend.position = "none") + xlab("") + ylab("") + 
	ggtitle("Reference Distributions of Posterior Predictive Network Statistics by Model") +
	geom_point(aes(y=real), colour='purple') +
	stat_summary(fun.y=add.percentiles, geom="point", pch="_", colour="grey40", size=6) 







dens <- subset(mon, statistic=='density')
deg <- subset(mon, statistic=='MeanDegree')

which(sim_stats$mean_deg != deg$value)
which(sim_stats$edg_dens != dens$value)


mean(sim_stats$mean_deg)
sd(sim_stats$mean_deg)
mean(sim_stats$edg_dens)
sd(sim_stats$edg_dens)




ggplot(mon, aes(y=value, x=statistic)) + geom_boxplot() + theme(legend.position = "none") +
	geom_point(aes(y=stats$edg_dens, x=factor("density"), colour='red')) + 
	geom_point(aes(y=stats$mean_deg, x=factor("MeanDegree"), colour='red')) +
	facet_wrap(~statistic)

ggplot(dens, aes(y=value, x=statistic)) + geom_boxplot() + theme(legend.position = "none") +
	geom_point(aes(y=stats$edg_dens, x=factor("density"), colour='red')) 

ggplot(deg, aes(y=value, x=statistic)) + geom_boxplot() + theme(legend.position = "none") +
	geom_point(aes(y=stats$mean_deg, x=factor("MeanDegree"), colour='red'))



ggplot(mon, aes(y=value, x=statistic)) + geom_point() + theme(legend.position = "none")



# Calculate CIs  ! WRONG DONT USE

n <- nrow(dens)
m1 <- mean(dens$value)
s1 <- sd(dens$value)
error1 <- qt(0.975,df=n-1) * s1/sqrt(n) # 95% CI using t-dist
lower1 <- m1 - error1
upper1 <- m1 + error1


m2 <- mean(deg$value)
s2 <- sd(deg$value)
error2 <- qt(0.975,df=n-1) * s2/sqrt(n) # 95% CI using t-dist
lower2 <- m2 - error2
upper2 <- m2 + error2

ci <- data.frame(statistic = levels(mon$statistic), mean = c(m1, m2), sd = c(s1, s2), 
								 error = c(error1, error2), lower = c(lower1, lower2), upper = c(upper1, upper2))

pd <- position_dodge(0.78)
ggplot(ci, aes(x=statistic, y = mean, group = statistic)) +
	#draws the means
	geom_point(position=pd) +
	#draws the CI error bars
	geom_errorbar(data=ci, aes(ymin=lower, ymax=upper), width=.1, position=pd) +
	theme(legend.position = "none") + ggtitle("95% Confidence Intervals") +
	geom_point(aes(y=stats$edg_dens, x=factor("density"), colour='red')) + 
	geom_point(aes(y=stats$avg_deg, x=factor("MeanDegree"), colour='red'))


ci_dens <- subset(ci, statistic=='density')
ggplot(ci_dens, aes(x=statistic, y = mean, group = statistic)) +
	#draws the means
	geom_point(position=pd) +
	#draws the CI error bars
	geom_errorbar(data=ci_dens, aes(ymin=lower, ymax=upper), width=.1, position=pd) +
	theme(legend.position = "none") + ggtitle("95% Confidence Interval") +
	geom_point(aes(y=stats$edg_dens, x=factor("density"), colour='red')) 


ci_deg <- subset(ci, statistic=='MeanDegree')
ggplot(ci_deg, aes(x=statistic, y = mean, group = statistic)) +
	#draws the means
	geom_point(position=pd) +
	#draws the CI error bars
	geom_errorbar(data=ci_deg, aes(ymin=lower, ymax=upper), width=.1, position=pd) +
	theme(legend.position = "none") + ggtitle("95% Confidence Interval") +
	geom_point(aes(y=stats$avg_deg, x=factor("MeanDegree"), colour='red'))
