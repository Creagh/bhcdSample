##############################################################################
# Posterior Plots of sampled parameters & hyperparameters
#
# 
##############################################################################
library(ggplot2)
library(hexbin)
library(ggpubr)


setwd("~/Desktop/computeCanada")

### TERRORIST ###
alpha<- read.csv("terrorist/2020-10-31-01-12-58-WaySFoKT.exec/samples/alpha.csv", header = TRUE)
beta <- read.csv("terrorist/2020-10-31-01-12-58-WaySFoKT.exec/samples/beta.csv", header = TRUE)

aZI<- read.csv("terrorist/2020-10-31-01-38-53-f00p1BP8.exec/samples/alpha.csv", header = TRUE)
bZI <- read.csv("terrorist/2020-10-31-01-38-53-f00p1BP8.exec/samples/beta.csv", header = TRUE)
theta<- read.csv("terrorist/2020-10-31-01-38-53-f00p1BP8.exec/samples/theta.csv", header = TRUE)

### GRASS WEB ###
alpha<- read.csv("grassWeb/2020-11-06-15-11-33-pYSkxguw.exec/samples/alpha.csv", header = TRUE)
beta <- read.csv("grassWeb/2020-11-06-15-11-33-pYSkxguw.exec/samples/beta.csv", header = TRUE)

aZI<- read.csv("grassWeb/2020-11-06-15-23-14-8FqK7S4a.exec/samples/alpha.csv", header = TRUE)
bZI <- read.csv("grassWeb/2020-11-06-15-23-14-8FqK7S4a.exec/samples/beta.csv", header = TRUE)
theta<- read.csv("grassWeb/2020-11-06-15-23-14-8FqK7S4a.exec/samples/theta.csv", header = TRUE)

# Combine samples into data frame

data <- data.frame(alpha=c(alpha$value, aZI$value), beta=c(beta$value, bZI$value), 
									 model=c(rep("Standard", nrow(alpha)), rep("ZI", nrow(aZI))))

dataZI <- data.frame(alpha=aZI$value, beta=bZI$value, theta=theta$value)

# 2D Density Plots

# Bin size control + color palette
ggplot(data, aes(x=alpha, y=beta) ) +	geom_bin2d(bins = 200) +	scale_fill_continuous(type = "viridis") + theme_bw() +
	ggtitle("Joint Posterior Density Plot by Model") + facet_wrap(~model)

# Hex bins
ggplot(data, aes(x=alpha, y=beta) ) +	geom_hex(bins = 150) +	scale_fill_continuous(type = "viridis") + theme_bw() +
	ggtitle("Joint Posterior Density Plot by Model") + facet_wrap(~model)

gA <- ggplot(dataZI, aes(x=alpha, y=theta) ) +	geom_hex(bins = 150) +	scale_fill_continuous(type = "viridis") + theme_bw()
gB <- ggplot(dataZI, aes(x=beta, y=theta) ) +	geom_hex(bins = 150) +	scale_fill_continuous(type = "viridis") + theme_bw()
ggarrange(gA, gB, ncol=2)

# Contours + density
ggplot(subset(data, model=="Standard"), aes(x=alpha, y=beta) ) + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +	
	scale_fill_continuous(type = "viridis") + ggtitle("Joint Posterior Density Plot - Standard Model")

ggplot(subset(data, model=="ZI"), aes(x=alpha, y=beta) ) + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +	
	scale_fill_continuous(type = "viridis") + ggtitle("Joint Posterior Density Plot - ZI Model")

ggplot(data, aes(x=alpha, y=beta) ) + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +	
	scale_fill_continuous(type = "viridis") + ggtitle("Joint Posterior Density Plot by Model") +
	facet_wrap(~model, scales="free")


gAC <- ggplot(dataZI, aes(x=alpha, y=theta) ) + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +	
	scale_fill_continuous(type = "viridis")
gBC <- ggplot(dataZI, aes(x=beta, y=theta) ) + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +	
	scale_fill_continuous(type = "viridis")
ggarrange(gAC, gBC, ncol=2)



############ test 3D plots

library(rayshader)

data <- data.frame(alpha = alpha$value, beta = beta$value)

ggCont <- ggplot(data, aes(x=alpha, y=beta) ) + stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
	scale_fill_continuous(type = "viridis")

(ggHex <- ggplot(data, aes(x=alpha, y=beta) ) +	geom_hex(bins = 20, size=0.25, color="black") +	scale_fill_viridis_c(option = "C") + theme_bw()+
		scale_x_continuous(limits = c(0.1, 0.6)) + scale_y_continuous(limits=c(0.1, 1.5)))

plot_gg(ggCont,multicore=TRUE,width=5,height=5,scale=250,windowsize=c(1400,866),
				zoom = 0.55, phi = 30)



par(mfrow = c(1, 1))
plot_gg(ggHex, width = 5, height = 4, scale = 300, multicore = TRUE, windowsize = c(1200, 960),
				fov = 70, zoom = 0.4, theta = 330, phi = 40)
Sys.sleep(0.2)
render_snapshot(clear = TRUE)
