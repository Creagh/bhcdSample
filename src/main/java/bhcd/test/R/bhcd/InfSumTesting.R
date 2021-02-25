library(ggplot2)

# Input variables
Er <- 600
Lr <- 30
Rr <- 20
alpha <- 0.1
beta <- 0.00000001

maxIter <- 30


# Produce plot
df <- produceParetoPlot(Er, Lr, Rr, alpha, beta, maxIter)




########################## FUNCTIONS

# Function to produce pareto style plot of sum and cumulative sum
produceParetoPlot <- function(Er, Lr, Rr, alpha, beta, maxIter) { 
	
	df <- data.frame(k = 1:maxIter)
	
	for(k in 1:maxIter) {
		# first, compute the summand on the log-scale
		logSummand <- lgamma(Er + alpha + (k * (Lr * Rr + 1)) ) - lgamma(Er + alpha + (k * (Lr * Rr + 1)) + beta + 1)
		df$summand[k] <- exp(logSummand)
	}
	
	# replace NaNs with zero
	df[is.nan(df)] <- 0
	
	# calculate cumulative sum
	df$cuSum <- cumsum(df$summand)
	
	
	# plot sum (pareto chart style)
	print(ggplot(df, aes(x=k)) + geom_bar(aes(y=summand), stat='identity') + geom_point(aes(y=cuSum), color='red') + 
		geom_path(aes(y=cuSum, group=1), colour="slateblue1", lty=3, size=0.9) +
		labs(title='Infinite Sum within Marginalized Likelihood', y='Summand Value', x='Index k') +
		scale_x_continuous(breaks = df$k) )
	
	# return df of values
	return(df)
	
}

# Helper function for finding NaN values
is.nan.data.frame <- function(x) {
	do.call(cbind, lapply(x, is.nan))
}