source('initialise.R')
source('config.R')
source('runInitialiseData.R')
cat(sprintf('Reading from %s\n',location.output))
calDat <- readRDS(file.path(location.output,'calDat.RDS'))$calDat
resSigma <- readRDS(file.path(location.output,'sigma-indepParms.RDS'))

# run files ####
# for input
location.runFiles <- file.path(location.output,'detectedParmSpace')
runFilesList <- list.files(location.runFiles,pattern = 'workUnit-[0-9]+\\.RDS')

if(length(runFilesList)==0){
	stop('no run files to process\nHave you run runMLEandParmSpace?\n')
}

# read ####
cat('reading sampleParms...')
sampleParms <- readRDS(file.path(location.output,'sampleParmsParscaleRanged.RDS'))
cat('done\nreading samplePoints...')
samplePoints <- as.data.frame(readRDS(file.path(location.output,'samplePoints.RDS')))
if(nrow(samplePoints)!=numSample){
	stop('Number of sample points and numSample do not match\n')
}
if(ncol(samplePoints)!=nrow(sampleParms)){
	integerParms <- read.csv('frida_integer_parms.csv')
	if(sum(integerParms$Variable%in%colnames(samplePoints) &&
									!integerParms$Variable%in%sampleParms$Variable)>0){
		sampleParms <- prepareSampleParms(sampleParms=sampleParms,integerParms = integerParms)
	}
}
cat('done\n')

# log like ####
cat('reading log likelihoods...\n')
logLike <- rep(NA,numSample)
completeRunsSoFar <- 0
for(f.i in 1:length(runFilesList)){
	cat(sprintf('\r chunk %i of %i',f.i,length(runFilesList)))
	parOutput <- readRDS(file.path(location.output,'detectedParmSpace',paste0('workUnit-',f.i,'.RDS')))
	for(l in 1:length(parOutput)){
		logLike[parOutput[[l]]$parmsIndex] <- parOutput[[l]]$logLike
		if(!is.na(parOutput[[l]]$runDat[[1]][length(parOutput[[l]]$runDat[[1]])])){
			completeRunsSoFar <- completeRunsSoFar + 1
		}
	}
	rm(parOutput)
}
cat(sprintf('\r read %i files, collected %i sample log likes, %i runs in data where complete\n',
						length(runFilesList),numSample,completeRunsSoFar))


samplePoints$logLike <- logLike

parVect <- sampleParms$Value
names(parVect) <- sampleParms$Variable
baseLogLike <- -negLLike(parVect)

logLike.badRM <- logLike
logLike.badRM[logLike.badRM <= -.Machine$double.xmax+.Machine$double.eps*(1+nrow(calDat))] <- NA
cat(sprintf('%i bad log likelihoods\n',sum(is.na(logLike.badRM))))
png(file.path(location.output,'logLikeHist.png'),
		width=10,height=10,res=150,units='cm')
hist(logLike.badRM,xlim=c(min(logLike.badRM,na.rm=T),baseLogLike))
abline(v=baseLogLike,col='red')

# plot log like accross all parms ####
sqrtNplots <- sqrt(nrow(sampleParms))
plotCols <- round(sqrtNplots)
plotRows <- ceiling(sqrtNplots)
png(file.path(location.output,'parameterLogLikelihoods.png'),
		width=plotCols*10,height=plotRows*10,res=150,units='cm')
funPlotParRangesLikelihoods(sampleParms,samplePoints=samplePoints,like=logLike.badRM,
														baseLike = baseLogLike,minLike = baseLogLike-abs(0.999*baseLogLike))#,ylim=range(logLike))
dev.off()


# experiments ####
if(F){ # experiments
	
	# log like ceteris paribusish ####
	dim(samplePoints)
	dim(sampleParms)
	par.i <- 1
	parTol <- (sampleParms$Max-sampleParms$Min)*0.1
	samplePointsOnlyVaryPari.idc <- which(colSums(t(samplePoints[,-par.i])>=(sampleParms$Value[-par.i]-parTol[-par.i]))==(nrow(samplePoints)-1) & 
																					t(samplePoints[1:10,-par.i])<(sampleParms$Value[-par.i])+parTol[-par.i])
	samplePointsOnlyVarPari <- samplePoints[samplePointsOnlyVaryPari.idc,]
	funPlotParRangesLikelihoodsI(par.i,sampleParms,sampleParms,samplePointsOnlyVarPari,
															 logLike[samplePointsOnlyVaryPari.idc])
	
	# calculate prob ###
	range(logLike)
	like <- (exp(logLike-max(logLike)+log(.Machine$double.xmax)))
	funPlotParRangesLikelihoodsI(1,sampleParms,sampleParms,
															 samplePoints,logLike,yaxPad=0.04,
															 outputFilePath=NULL,
															 outputFileName='parameterLogLikelihoods.png',
															 NULL,
															 includeZeroYVal=FALSE,logY=F,ylim=NULL,
															 like.range=NULL,
															 minLike=NULL, parallelPlot = F,
															 plotWidth = 10,
															 plotHight = 10,
															 plotRes = 150,
															 plotUnits = 'cm')
	
	
	rangeLogLike <- range(logLike,na.rm=T)
	logLikeCompressed <- (logLike-rangeLogLike[1])/(rangeLogLike[2]-rangeLogLike[1])
	like <- exp(logLikeCompressed)^(1/(rangeLogLike[2]-rangeLogLike[1]))
	hist(logLikeShifted)
	logLikeShifted <- logLike/max(logLike,na.rm=T)
	like <- exp(logLikeShifted)
	likeSum <- sum(like)
	prob <- like/likeSum
	# assuming each samplePoint Likelihood has the same weight in the distribution
	prob <- as.double(like/likeSum)
	
	samplePoints$like <- like
	samplePoints$prob <- prob
	
	# 
	# plot prob ####
	png(file.path(location.output,'parameterLogLikelihoods.png'),
			width=plotCols*10,height=plotRows*10,res=150,units='cm')
	funPlotParRangesLikelihoods(sampleParms,samplePoints=samplePoints,like=prob)
	dev.off()
	
	
	# par.i <- 1
	# samplePoints.sorted <- sort_by(samplePoints,samplePoints[[par.i]])
	# ecdf.byPar <- list()
	# ecdf.byPar[[colnames(samplePoints.sorted)[par.i]]] 
	# ecdf.par.i <- data.frame(parValue=samplePoints.sorted[[par.i]])
	# colnames(ecdf.par.i)[1] <- colnames(samplePoints.sorted)[par.i]
	# ecdf.par.i$cdf <- NA
	# ecdf.par.i$cdf[1] <- samplePoints.sorted$prob[1]
	# for(r.i in 2:nrow(samplePoints.sorted)){
	# 	ecdf.par.i$cdf[r.i] <- 
	# 		samplePoints.sorted$prob[r.i] + 
	# 		ecdf.par.i$cdf[r.i-1]
	# }
}							
