
################################ PROJECT ADV DATA ANALYSIS ################################

wind <- read.csv("/Adv Data Analysis/Project/HornsRevData.csv")
turb <- read.csv("/Adv Data Analysis/Project/hr2_turbines_98.csv")
pred <- read.csv("/predictionData.csv")
wind$Density <- wind$Nhat/wind$area
pred$day <- 19
pred$month <- 11

# EDA
par(mfrow= c(2,2))
plot(wind$x.pos, wind$Nhat, xlab= "x.pos", ylab = "Counts")
plot(wind$y.pos, wind$Nhat, xlab= "y.pos", ylab = "Counts")
plot(wind$month, wind$Nhat, xlab= "month", ylab = "Counts")
plot(wind$Depth, wind$Nhat, xlab= "depth", ylab = "Counts")
plot(as.factor(wind$impact), wind$Nhat, xlab= "impact", ylab = "Counts")
plot(wind$day, wind$Nhat, xlab= "day", ylab = "Counts")
par(mfrow= c(1,1))
library(fields)
quilt.plot(wind$x.pos[wind$impact==0], wind$y.pos[wind$impact==0], wind$Density[wind$impact==0])
quilt.plot(wind$x.pos[wind$impact==1], wind$y.pos[wind$impact==1], wind$Density[wind$impact==1])

# -------------------------------------------------------------------------------------------- #

# Fitting a CReSS using SALSA2D: see notes on paper
library(MRSea)
wind$foldid <- getCVids(wind, folds = 2, block = NULL)
wind$response <- wind$Nhat
bind <- cbind(wind$x.pos, wind$y.pos)
knotgrid <- getKnotgrid(coordData = bind, plot = F)
#OneD model
splineParams <- makesplineParams(data= wind, varlist = c('Depth', 'day'), predictionData = pred)
initialModel <- glm(response~impact+as.factor(month)+offset(log(area)), 
                    data = wind, family = "quasipoisson")

salsa1dlist <- list(fitnessMeasure= "QAIC", minKnots_1d= c(1, 1), maxKnots_1d= c(10, 10), 
                    startKnots_1d= c(4, 8), degree= c(2, 2), maxIterations= 10, gaps= c(0, 0))
salsa1dOutput <- runSALSA1D_withremoval(initialModel, salsa1dlist, varlist=c('Depth', 'day'),
                                        factorlist = c('impact', 'month'), 
                                        splineParams = splineParams, datain = "wind")

distMats <- makeDists(cbind(wind$x.pos, wind$y.pos), na.omit(knotgrid))
r_seq <- getRadiiChoices(numberofradii = 10, distMatrix = distMats$dataDist)

salsa2dlist <- list(fitnessMeasure= "QAIC", knotgrid= knotgrid, knotdim= c(100,100),
                   minKnots= 4, maxKnots= 12, r_seq= r_seq, gap= 0, 
                   interactionTerm= "impact", startKnots= 10)
splineParams <- salsa1dOutput$splineParams
updatedModel <- runSALSA2D(salsa1dOutput$bestModel, salsa2dlist, d2k = distMats$dataDist,
                           k2k= distMats$knotDist, splineParams = salsa1dOutput$splineParams)

library(car)
pvals <- Anova(updatedModel$bestModel, test= "F")
rownames(pvals) <- c("impact", "month", "s(Depth)", "s(day)", "s(x.pos, y.pos)", "s(x.pos, y.pos):impact",
                     "Residuals")

# Should try CV to see if we got stuck in a local minimum (1D)
cv1 <- getCV_CReSS(wind, salsa1dOutput$bestModel, salsa1dOutput$splineParams)


# -------------------------------------------------------------------------------------------- #
# Model Diagnostic

# Mean-variance relationship
cut.fit <- cut(fitted(updatedModel$bestModel), breaks= quantile(fitted(updatedModel$bestModel), 
                                                               probs = c(seq(0,1, length=40))))
means <- tapply(fitted(updatedModel$bestModel), cut.fit, mean)
vars <- tapply(residuals(updatedModel$bestModel, type = "response"), cut.fit, var)
plot(means, vars, xlab = "Fitted Values", ylab= "Variance of Residuals")
lines(fitted(updatedModel$bestModel), fitted(updatedModel$bestModel)*606.97, col= "blue")
# The increasing and overdispersed mean-var relationship seems respected for fitted values
 


# Check p-values of covariates
# See text file pvalues_gee.csv


# Collinearity
coll <- vif(updatedModel$bestModel)
rownames(coll) <- c("impact", "month", "s(Depth)", "s(day)", "s(x.pos, y.pos)", "s(x.pos, y.pos):impact")
coll # Huge collinearity issues. May solve with scaling the x.pos/y.pos ---> useless: VIFs increased
# But I should center the spline not the individual x.pos and y.pos!

library(lawstat)
splineParams <- updatedModel$splineParams
radiusIndices <- splineParams[[1]]$radiusIndices
dists <- splineParams[[1]]$dist
radii <- splineParams[[1]]$radii
aR <- splineParams[[1]]$invInd[splineParams[[1]]$knotPos]



inter <- update(updatedModel$bestModel, .~. -LocalRadialFunction(radiusIndices, dists, 
                                            radii, aR)-impact:LocalRadialFunction(radiusIndices, 
                                            dists, radii, aR), data= wind)
inter <- update(inter, .~. +scale(LocalRadialFunction(radiusIndices, dists, 
                          radii, aR), scale= F)+impact:scale(LocalRadialFunction(radiusIndices, 
                          dists, radii, aR), scale= F), data= wind, x= T, y=T)
# Rechecking collinearity
coll1 <- vif(inter)
rownames(coll1) <- c("impact", "month", "s(Depth)", "s(day)", "scaled(s(x.pos, y.pos))", 
                     "scaled(s(x.pos, y.pos)):impact")
coll1


# Inependence of residuals
plotRunsProfile(inter, varlist= c("Depth", "day"))
# The residuals are not independent. Should model this dependence through a GEE
# ACF to chose the blocking structure
acf(residuals(inter)) 
runACF(wind$TNO, inter, store = F) # The correlation declines, hence TNO blocking structure works

# Building a GEE
baseModel <- update(inter, .~.,  data= wind)
library(geepack)
gee <- geeglm(formula(inter), data = wind, family = poisson,
              id = wind$TNO, corstr = "ar1", x= T) 

# Calculating QICu
library(MuMIn)
QIC(gee)

# Rechecking autocorrelation: 
plotRunsProfile(gee, varlist= c("Depth", "day")) 
acf(residuals(gee)) 
runACF(wind$TNO, gee, store = F) # The model seems still highly correlated: # AR(1) not enough


# Linearity on the link-scale: covariate
runPartialPlots(gee, wind, factorlist.in = c("month", "impact"), varlist.in = c("Depth", "day"))
# Huge CIs for Depth and day ---> try to include them as linear terms
gee.DepthDay.lin <- update(gee, .~. -bs(Depth, knots = splineParams[[2]]$knots, 
                        degree = splineParams[[2]]$degree, Boundary.knots = splineParams[[2]]$bd)+Depth)
# Recheck linearity on the link scale for this updated model
runPartialPlots(gee.DepthDay.lin, wind, factorlist.in = c("month", "impact"), varlist.in = c("Depth", "day"))
QIC(gee.DepthDay.lin) # The GEE with Depth as linear has much lower QIC.


# Influential points
plot(cooks.distance(inter), ylab= "Cook's Distance")
# Highly influential points between observations N° 15000-20000

# -------------------------------------------------------------------------------------------- #
# Checking fitted values

library(fields)
quilt.plot(wind$x.pos[wind$impact==0], wind$y.pos[wind$impact==0], 
           fitted(gee, type= "response")[wind$impact==0], nrow=25, ncol=60, 
           zlim=range(fitted(gee, type= "response")))

quilt.plot(wind$x.pos[wind$impact==1], wind$y.pos[wind$impact==1], 
           fitted(gee, type= "response")[wind$impact==1], nrow=25, ncol=60, 
           zlim=range(fitted(gee, type= "response")))
points(turb$X, turb$Y, type = "p")

# Observed vs Fitted
plot(wind$Nhat, predict(gee, type= "response"))
abline(0, 1, col= "red")
# Awful fit


# CI for impact parameter
exp(c(1.3504 - qt(0.975, 27853)*0.2799, 1.3504 +  qt(0.975, 27853)*0.2799))

# P-values
getPvalues(gee, varlist = c("Depth", "day"), factorlist = c("month","impact"))


# Cumulative Residuals
plotCumRes(gee, varlist=c("Depth"), splineParams=splineParams, d2k=dists)





