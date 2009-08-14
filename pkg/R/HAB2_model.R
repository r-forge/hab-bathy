###  HAB 2 model

### 	R implementation of the remote sensing based Hydraulically Assisted Bathymetry models (HAB) developed by Fonstad, M. & Marcus, W. (2005). Functions are provided for modelling with a range of custom variables. Summary functions are provided as well.

###   Fonstad, M. & Marcus, W. (2005), ‘Remote sensing of stream depths with hydraulically assisted bathymetry (HAB) models’, Geomorphology 72, 320–339.



library(car)
#library(caret) # includes postResample
#library(NCStats) # includes fit.plot

printGG <- 
function(filename="R_Graph")
{
	dev.print(pdf, file=paste(filename,"_7x4.pdf"), width=7, height=4, pointsize=10)
	dev.print(pdf, file=paste(filename,"_7x7.pdf"), width=7, height=7, pointsize=10)
	dev.print(pdf, file=paste(filename,"55x55.pdf"), width=5.5, height=5.5, pointsize=10)
	return("Printed!")
	}
	
	par(las=1) # axis lables horzontal

# Q = discharge (m^3/s)
# V = average velocity (m/s)
# A = cross section area (m^2)     = W * Da
# W = width (m)
# Da = average depth (m)
# Dmin = minimum depth (cm) = 5 cm
# Dmax = maximum depth (cm) = 2 * Da
# R = Hydraulic radius (m)     = Da
# P = wetted perimeter (m)
# S = longitudinal energy gradient of flow (m/m)
# n = hydraulic resistance = 0.033 (\ref{NVE})
#     HAB2 variables: 
# I = intensity of light at some depth = DN
## I0 = intensity of light immediately prior to entering water = DN0
## 		DN0 = DN of shallowest depth or shore
# 		B10 = max(Lines$max.B1)
# 		B20 = max(Lines$max.B2)
# 		B30 = max(Lines$max.B3)
# beta = diffuse attenuation coefficient
# D = distance of light travelled through water = Depth

#  Q = A * V = W * Da * V 
#
#  V = R^(2/3) * S^(0.5) / n
#
#  assuming: R = Da
#
#  Q = W * (Da^1.83) * (S^0.12) / 0.32
#
#  Da = (Q / (3.125 * W * S^0.12))^0.55

# Wetted Perimeter:
	# Aw <- sqrt((D2 - D1)^2 + W^2)
	#where:
	#Aw = wetted area between two points
	#D1 = depth at first point
	#D2 = depth at second point
	#W = distance between measured points

# The tab-delimited data.frame "Lines" should be of the following structure:
# min.B1 avg.B1 max.B1 min.B2 avg.B2 max.B2 min.B3 avg.B3 max.b3 Length LineNo
# Length (lenght of the cross-section) should be in meters

"lines.data" <- read.table("Lines.csv", header=TRUE, sep="\t", na.strings="NA", dec=".", , strip.white=TRUE)

	
#######################################################

#Q <- 194.5
#n <- 0.033
#S <- 5/2700
#
#beta = 0.5 # seed value
#point.dist = 1 # distance between points in meters (defaults to 1)
#
#data <- lines.data
#band = 1
#CrsSec = 1


HAB2.depth <- function(data, CrsSec, band, DN.0, Q, n, S, beta, point.dist, auto = FALSE, verbose = FALSE)
{
	
	depth.estimate <- function(DN, DN0, b )
	{
		De <- log(DN/DN0)/-b
		return(De)
		}
	
	wetted.perim<- function(D = De, W = width, p.dist = point.dist)
	{
		Aw <- W - length(D) # adjust for sub-meter accuracy
		# compute first step (shore - D[i])
		Awi <- sqrt((D[1] - 0)^2 + p.dist^2)
		Aw <- Aw + Awi
		for (i in 2:length(D)-1)
		{
			Awi <- sqrt((D[i+1] - D[i])^2 + p.dist^2)
			Aw <- Aw + Awi
			}
		#compute last step
		Awi <- sqrt((0 - D[length(D)])^2 + p.dist^2)
		Aw <- Aw + Awi
		
		P <- Aw
		return(P)
		}

	velocity.estimate <- function(R = hydr.radius, S ,n )
	{
		V <- R^(2/3) * S^(0.5) / n
		return(V)
		}
	
	discharge <- function(W = width, Da = depth.avg, V = V.est)
	{
		Q <- W * Da * V
		Q.round <- round(Q, digits=1)
		return(Q.round)
		}
	
	line <- subset(data, data$LineNo == CrsSec )
	
	# get DN0 values automatically if none given
	if (length(DN.0) < 1)  
		DN0 <- max(line[2 + band])
	else
		DN0 <- DN.0[band]
	
	width = line[1,2] # in meters
	
#	estimate beta until Q is correct	
#  crude target finding method, but it works
	Q.est = 0
	iterations = 1
	while (Q.est != Q)
	{
		if (Q.est > Q) 
			beta <- beta * 1.2 
		else 
			beta <- beta / 1.1

		D.est <- mapply( depth.estimate, list(line[2 + band]), 					MoreArgs=list(DN0 = DN0, b = beta))
		De <- D.est[[1]]
	
		depth.avg <- mean(De)  # in meters
		A <- width * depth.avg  # cross sectional area in m^2

		P <- wetted.perim(D = De, W = width, p.dist = point.dist)
		hydr.radius <- A / P	 # Hydraulic radius in meters

		V.est <- velocity.estimate(R=hydr.radius, S=S, n=n)

		Q.est <- discharge(W = width, Da = depth.avg, V = V.est)

		iterations <- iterations + 1

		}
	De.round <- round(De,digits=2)
	De.round.mean <- round(mean(De),digits=2)
	beta.round <- round(beta, digits=3)
	DN0.mean <- round(mean(DN0),digits=0)
	if (verbose == TRUE)
	{
		print(paste("Q estimate:                ", Q.est), quote=F)
		print(paste("beta:                      ", beta.round),quote=F)
		print(paste("mean depth estimate:       ", De.round.mean," m"), quote=F)
		print(paste("max depth estimate:        ", max(De.round)," m"), quote=F)
		print(paste("Iterations:                ", iterations), quote=F)
		print(paste("DN0:                       ", DN0), quote=F)
		}
	return(list(De.round, beta.round, De.round.mean, max(De.round), 					DN0.mean))

	}
	
HAB2.model <- function(data = lines.data, DN.0 = c(), CrsSec = 1, band = 1, Q, n, S, beta = 0.5, point.dist = 1, auto = FALSE, verbose = FALSE)
{
	if (auto == TRUE)
	{
		if (verbose == FALSE)
			print(paste("thinking..."), quote=F)
		# matrix length = number of cross sections + mean De + max De
		#times number of bands -2 (-lineNo and LineLength columns)

		CS.results <- matrix(nrow = 4, ncol = length(data[1,])-2) 
		# 3 rows: avg beta, avg mean De, avg max De, mean DN0
		
		beta.bands <- matrix(nrow = max(data$LineNo), ncol = length(data[1,])-2) 	
		D.est.means <- matrix(nrow = max(data$LineNo), ncol = length(data[1,])-2)
		D.est.max <- matrix(nrow = max(data$LineNo), ncol = length(data[1,])-2)
		DN0.means <- matrix(nrow = max(data$LineNo), ncol = length(data[1,])-2)
	
		for (j in 1:max(data$LineNo))    # cross sections
			{
				if (verbose == TRUE)
					print(paste("********* CrossSection: ",j), quote=F)
				for (k in 1:(length(data[1,])-2))    # bands
					{
						if (verbose == TRUE)
							print(paste("band:  ",k), quote=F)
						model <- HAB2.depth(data, CrsSec = j, band = k, 											DN.0 = DN.0, Q=Q, n=n, S=S, beta=beta, 											point.dist=point.dist, verbose=verbose)
						beta.bands[j,k] <- model[[2]]
						D.est.means[j,k] <- model[[3]] 
						D.est.max[j,k] <- model[[4]] 
						DN0.means[j,k] <- model[[5]]
						}
				}
			# include mean De in second-to-last column
			# Include max De in last column
			for (a in 1:length(CS.results[1,]))
			{
				CS.results[1,a] <- mean(beta.bands[,a], na.rm=T)
				CS.results[2,a] <- mean(D.est.means[,a], na.rm=T)
				CS.results[3,a] <- max(D.est.max[,a], na.rm=T)
				CS.results[4,a] <- mean(DN0.means[,a], na.rm=T)
				}

		return(CS.results)		
		}
	else
	{
		model <- HAB2.depth(data=data, CrsSec=CrsSec, band= band, DN.0 = DN.0 								, Q=Q, n=n, S=S, beta=beta, point.dist=point.dist, 								verbose=verbose)
		return(model)
		}
	}

#### Example:	
#Call either as:
#	HAB2 <-HAB2.model(lines.data, CrsSec = 1, band = 1, Q = 194.5, n = 0.033, S=5/2700, DN.0 = c())
#	
#or with set to automatically step through all bands and cross sections: 	HAB2 <-HAB2.model(lines.data, Q = 194.5, n = 0.033, S=5/2700, DN.0 = c(), auto=TRUE)

#Get the means of beta for each band:
#				A0033 <- apply(HAB2.0033, 2, mean)	


######## Summary function

HAB2.summary <- function(model)
{
	if (is.list(model))
	{
		print(paste("beta:    DN0:   mean depth estimate:   max depth estimate:"), quote=F)
		print(paste(model[[2]],"  ",round(model[[5]],digits=0),"  ",model[[3]]," m               ", model[[4]]	," m"), quote=F)
		}
	else
		
	if (is.matrix(model))
	{
		Tm <- t(model)
		TmF <- format(Tm, digits=3)
		print(paste("Summary of HAB2 matrix (means for all cross sections)"), quote=F)
		print(paste("band:   beta:    DN0:   mean depth estimate:   max depth estimate:"), quote=F)
		for (i in 1:length(Tm[,1]))
		{
			print(paste(i,"   ",TmF[i,1],"  ",round(Tm[i,4],digits=0),"",TmF[i,2]," m            ",TmF[i,3]," m"), quote=F)
			
			}
		}
	else
		print(paste("ERROR: model to summarise is neither a list nor a matrix"), quote=F)	
	}	


#### Example:
# Call either one band at a time explicitly:
#	HAB2 <-HAB2.model(lines.data, CrsSec = 1, band = 1, Q = 194.5, n = 0.033, S=5/2700)
#	
# or set to automatically step through all bands and cross sections: 	HAB2 <-HAB2.model(lines.data, Q = 194.5, n = 0.033, S=5/2700, auto=TRUE)

#Get the means of beta for each band:
#				A0033 <- apply(HAB2.003, 2, mean)

# or call manually
#		CS.bands <- c(1:15) 
#		 dim(CS.bands) <- c(5,3)
#	
#for (j in 1:5)
#	{
#		
#		print(paste("\r CrossSection: ",j))
#		for (k in 1:3)
#		{
#			print(paste("band:     ",k))
#			
#			HAB2 <-HAB2.depth(lines.data, CrsSec = j, band = k, n = 0.033, S=5/2700)
#			CS.bands[j,k] <- HAB2[[2]]
#			}
#
#		}	







#Example:
# First run the model:
#HAB2.0033 <-HAB2.model(lines.data, Q = 194.5, n = 0.033, S=5/2700, DN.0 = c(), auto=T, verbose=F)
# then get the summary:
#HAB2.summary(HAB2.0033)


