###  HAB 1 model

### 	R implementation of the remote sensing based Hydraulically Assisted Bathymetry models (HAB) developed by Fonstad, M. & Marcus, W. (2005). Functions are provided for modelling with a range of custom variables. Summary functions are provided as well.

###   Fonstad, M. & Marcus, W. (2005), ‘Remote sensing of stream depths with hydraulically assisted bathymetry (HAB) models’, Geomorphology 72, 320–339.


library(car)
library(caret) # includes postResample
library(NCStats) # includes fit.plot

printG <- 
function(filename="R_Graph")
{
	dev.print(pdf, file=paste(filename,"_7x4.pdf"), width=7, height=4, pointsize=10)
	dev.print(pdf, file=paste(filename,"_7x7.pdf"), width=7, height=7, pointsize=14)
	dev.print(pdf, file=paste(filename,"_55x55.pdf"), width=5.5, height=5.5, pointsize=10)
	return("Printed!")
	}
	
par(las=1) # axis lables horizontal

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

#  Q = A * V = W * Da * V 
#
#  V = (R^(2/3) * S^(0.5)) / n
#
#  assuming: R = Da
#
#  Q = W * (Da^1.83) * (S^0.12) / 0.32
#
#  Da = (Q / (3.125 * W * S^0.12))^0.55

# The data.frame "Lines" should be of the following structure:
# min.B1 avg.B1 max.B1 min.B2 avg.B2 max.B2 min.B3 avg.B3 max.b3 Length LineNo
# Length (lenght of the cross-section) should be in meters

"Lines" <- read.table("Lines.csv", header=TRUE, sep="\t", na.strings="NA", dec=".", , strip.white=TRUE)

# S <- height/distance in metres
# Dmin <- 5 in cm

HAB1.model <- function(CrossSecs=Lines, Q, n, S, Dmin, writeTables = FALSE)
{

	Depth.avg <- function( Q, n, W, S)
	{
		Da <- (W * S^(1/2) / (n * Q))^(-3/5)
		return(Da)
		}		

	HAB1.avg.max.depth <- function(CrossSecs, Q, n, S) 
	{
		for(i in 1:length(CrossSecs$LineNo))
		{
			Daa <- Depth.avg(Q,n,CrossSecs$Length[i],S)
			CrossSecs$Davg[i] <- (Daa) * 100   # m to cm
			CrossSecs$Dmax[i] <- CrossSecs$Davg[i] * 2
			}
		return(CrossSecs)
		}

	Lines <- HAB1.avg.max.depth(CrossSecs=CrossSecs, Q=Q, n=n, S=S)

	Depth <- c(Lines$Dmax, Lines$Davg, rep(Dmin,length(Lines$Line)))
	B1 <- c(Lines$min.B1, Lines$avg.B1, Lines$max.B1)
	B2 <- c(Lines$min.B2, Lines$avg.B2, Lines$max.B2)
	B3 <- c(Lines$min.B3, Lines$avg.B3, Lines$max.B3)

	HAB1.DB123 <- cbind(Depth,B1,B2,B3)

if (writeTables == TRUE)
	{
		write.table(Lines, paste("CrsSec",length(CrossSecs[,1]),"Lines model n",n,".csv"), col.names=TRUE, sep="\t", dec=".")	
		write.table(HAB1.DB123, paste("HAB1 DB123",length(CrossSecs[,1]),"n",n,".csv"), col.names=TRUE, sep="\t", dec=".")
		}

	plot(B1,Depth)
	par(new=T)
	plot(smooth.spline(B1,Depth,spar=0.8), type="l", ann=F, xaxt="n", yaxt="n")
	
	return(HAB1.DB123)
	
	}

# Example:
# HAB1.0033 <- HAB1.model(Q=194.5, n=0.03, S=5/2700, Dmin=5, writeTables=T)





HAB1.auto.fit <- function(model = HAB1.model, printGraphs = FALSE, legend = FALSE, n, A = 10, B = -0.001) 
# alternatively A=50, B=-0.04 has worked well in some cases
{
	for (i in 1:(length(model[1,])-1))
	{
		X <- model[,i+1]
		nls.fit <- nls(model[,1] ~ A * exp( B * X ), start=list(A=A,B=B))
		print(paste("*********** Band",i,"***********"),quote=F)
		print(summary(nls.fit))
		ab <- coef(nls.fit)
			round.a <- round(ab[1],digits=4)
			round.b <- round(ab[2],digits=4)
		print(paste("y =",round.a,"* exp(",round.b,"* B",i,")"),quote=F)
		plot(model[,i+1],model[,1], xlab=paste("Band",i), ylab="Depth")
		par(new=T)
		curve(ab[1] * exp( ab[2] * x ),ann=F, xaxt="n", yaxt="n", xlim=c(0,255))
		if (legend == TRUE)
		{
			legend(legend=c(paste("y =",round.a,"exp(",round.b,"x)")), 			x="topright", box.lty=0, inset=0.025)
			}
		
		if (printGraphs == TRUE)
			printG(paste("nls.D-B",i,"y=a*exp(b*x) n",n))
		}
	}
	
	
HAB1.fit <- function(model = HAB1.model, band, A = 50, B = -0.04, printGraphs = FALSE, legend= FALSE, n )
{

	BN <- model[,band+1]
	nls.fit <- nls(model[,1] ~ A * exp( B * BN ), start=list(A=A,B=B))
	print(paste("*********** Band",band,"***********"),quote=F)
	print(summary(nls.fit))
	ab <- coef(nls.fit)
		round.a <- round(ab[1],digits=4)
		round.b <- round(ab[2],digits=4)
	print(paste("y =",round.a,"* exp(",round.b,"* B",band,")"),quote=F)
	plot(model[,band+1],model[,1], xlab=paste("Band",band), ylab="Depth")
	par(new=T)
	curve(ab[1] * exp( ab[2] * x ),ann=F, xaxt="n", yaxt="n", xlim=c(0,255))
	if (legend == TRUE)
	{
		legend(legend=c(paste("y =",round.a,"exp(",round.b,"x)")),		x="topright", box.lty=0, inset=0.025)
		}
		
	if (printGraphs == TRUE)
		printG(paste("nls.D-B",band,"y=a*exp(b*x) n",n))
	}


#
################# Band 1 ################################
#
#
#nls.B1.exp <- nls(Depth ~ A * exp( B * B1 ), start=list(A=50,B=-0.04))
#summary(nls.B1.exp)
#ab.b1.exp <- coef(nls.B1.exp)
#plot(B1,Depth)
#par(new=T)
#curve(ab.b1.exp[1] * exp( ab.b1.exp[2] * x ),ann=F, xaxt="n", yaxt="n", xlim=c(0,255))
#printGG("nls.D-B1.exp_y=a*exp(b*x)")
#
#
#
################# Band 2 ########################################
#
#
#nls.b2.exp <- nls(Depth ~ A * exp( B * B2 ), start=list(A=2,B=0.05))
#summary(nls.b2.exp)
#ab.b2.exp <- coef(nls.b2.exp)
#plot(B2,Depth)
#par(new=T)
#curve(ab.b2.exp[1] * exp( ab.b2.exp[2] * x ),ann=F, xaxt="n", yaxt="n")
#printGG("nls.D-b2.exp_y=a*exp(b*x)")
#
#
################# Band 3 ######################################
#
#
#nls.B3.exp <- nls(Depth ~ A * exp( B * B3 ), start=list(A=2,B=0.05))
#summary(nls.B3.exp)
#ab.b3.exp <- coef(nls.B3.exp)
#plot(B3,Depth)
#par(new=T)
#curve(ab.b3.exp[1] * exp( ab.b3.exp[2] * x ),ann=F, xaxt="n", yaxt="n")
#printGG("nls.D-B3.exp_y=a*exp(b*x)")













