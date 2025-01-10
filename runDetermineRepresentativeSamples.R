# The likelihood that a given run is "close" to the desired quantile of all runs in all variable
# dimensions is the likelihood that the difference between this run and the desired quantile is
# the same as the difference between the desired quantile and other runs in gerneral.
# 
# How to account for quantiles different from 0.5 where the differences are assymmetric?
# The likelihood for the quantiles is the likelihood that the errors are from 
# a distribution with the variance at median and a mean of the difference between the
# desired quantile and the median

source('config.R')
source('runInitialiseData.R')
CIsToPlot <- c(0,0.5,0.95)
varsToRead <- colnames(calDat)

cat('reading sampleParms...')
sampleParms <- readRDS(file.path(location.output,'sampleParmsParscaleRanged.RDS'))
cat('done\nreading samplePoints...')
samplePoints <- as.data.frame(readRDS(file.path(location.output,'samplePoints.RDS')))
if(nrow(samplePoints)!=numSample){
	stop('Number of sample points and numSample do not match\n')
}
cat('done\n')
# for input
location.runFiles <- file.path(location.output,'detectedParmSpace')
runFilesList <- list.files(location.runFiles,pattern = 'workUnit-[0-9]+\\.RDS')
if(length(runFilesList)==0){
	stop('no run files to process\nHave you run runMLEandParmSpace?\n')
}

# collect time series ####
defRun <- runFridaDefaultParms()
yearsToRead <- rownames(defRun)
varsToRead.lst <- list()
workUnitBoundaries <- seq(1,ncol(calDat)+1,100000000)
workUnitBoundaries <- c(workUnitBoundaries,ncol(calDat)+1)
for(i in 1:(length(workUnitBoundaries)-1)){
	varsToRead.lst[[i]] <- colnames(calDat)[workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)]
}
vars.i <- 1
cat(sprintf('Reading %i vars: %i to %i of %i total vars...\n   ',
						length(varsToRead),workUnitBoundaries[vars.i],workUnitBoundaries[vars.i+1]-1,ncol(calDat)))
cat(paste0(varsToRead,collapse='\n   '))
cat('\n')
# dimensions time in rows, run IDs in columns,  variables to read
# read all the years, selectively plot later
runsData <- array(NA,dim=c(nrow(defRun),nrow(samplePoints),length(varsToRead)))
dimnames(runsData) <- list(rownames(defRun),1:nrow(samplePoints),varsToRead)
for(f.i in 1:length(runFilesList)){
	cat(sprintf('\r  reading chunk %i of %i',f.i,length(runFilesList)))
	parOutput <- readRDS(file.path(location.output,'detectedParmSpace',paste0('workUnit-',f.i,'.RDS')))
	for(l in 1:length(parOutput)){
		runsData[,parOutput[[l]]$parmsIndex,] <- unlist(parOutput[[l]]$runDat[,varsToRead])
	}
	rm(parOutput)
}
cat('\r  reading done                                                            \n')

# prep quantiles ####
# determine the values of the desired quantiles in all of the variables
plotWeightType <- 'equaly'
if(!exists('logLike')&&plotWeightType %in% c('likelihood','logCutoff','linearly')){
	# log like ####
	cat(' reading log likelihoods...\n')
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
	logLike.ecdf <- ecdf(logLike)	
}
if(plotWeightType=='likelihood'){
	samplePoints$plotWeight <- exp(logLike)
} else if(plotWeightType == 'logCutoff'){
	# somewhat wrong likelihood weighting
	samplePoints$plotWeight <- logLike.ecdf(logLike)
} else if(plotWeightType == 'equaly'){
	# equal weighting
	samplePoints$plotWeight <- rep(1,nrow(samplePoints))
} else if(plotWeightType == 'linearly'){
	samplePoints$plotWeight <- order(logLike)/nrow(samplePoints)
} else {
	stop('unknown plotWeightType\n'	)
}
# median ####
medians <- array(NA,dim=c(nrow(defRun),length(varsToRead)))
for(var.i in 1:length(varsToRead)){
	for(year.i in 1:nrow(defRun)){
		medians[year.i,var.i] <- Quantile(runsData[year.i,,var.i],
																			 weights = samplePoints$plotWeight,
																			 probs = 0.5,
																			 na.rm = T)
	}
}
# variances ####
if(!treatVarsAsIndep){
	stop('need to think about this for dependent vars\n')
}
stdDevs <- rep(NA,length(varsToRead))
names(stdDevs) <- varsToRead
for(var.i in 1:length(varsToRead)){
	errors <- array(NA,dim=c(nrow(defRun),nrow(samplePoints)))
	for(year.i in 1:length(varsToRead)){
		errors[year.i,] <- runsData[year.i,,var.i] - medians[year.i,var.i]
	}
	stdDevs[var.i] <- sd(as.vector(errors),na.rm=T)
}
rm(errors)

# quantile likelihood ###
# per sample log likelihood for each quantile
ciBoundQs <- unique(c(rev((1-CIsToPlot)/2),1-(1-CIsToPlot)/2))
ciLogLikes <- array(0,dim=c(length(ciBoundQs),length(varsToPlot),nrow(samplePoints)))
for(q.i in 1:length(ciBoundQs)){
	for(var.i in 1:length(varsToRead)){
		ciBounds <- c()
		errors <- array(NA,dim=c(nrow(defRun),nrow(samplePoints)))
		for(year.i in 1:nrow(defRun)){
			ciBounds[year.i] <- Quantile(runsData[year.i,,var.i],
																		weights = samplePoints$plotWeight,
																		probs = ciBoundQs[q.i],
																		na.rm = T)
			errors[year.i,] <- runsData[year.i,,var.i] - medians[,var.i]
		}
		for(sample.i in 1:nrow(samplePoints)){
			# ciLogLikes[q.i,var.i,sample.i] <- sum(pnorm(errors[,sample.i],
			# 								 mean = medians[,var.i]-ciBounds,sd = stdDevs[var.i],log=T))
			ciLogLikes[q.i,var.i,sample.i] <- sum(samplePoints$plotWeight[sample.i] * (pnorm(errors[,sample.i],
											 mean = medians[,var.i]-ciBounds,sd = stdDevs[var.i])-ciBoundQs[q.i])^2)
		}
	}
}
cat('\n    ')
means.store <- means
ciBounds.store <- ciBounds


# the likelihood for the quantiles is the likelihood that the errors are from 
# a distribution with the variance at median and a mean of the difference between the
# desired quantile and the median
if(F){
	# for testing, ylim values correspond to testing with var.i 132 (STA)
	plot(rownames(defRun),runsData[,1,var.i],ylim=c(0,10),type='l')
	for(i in 1:100){
		lines(rownames(defRun),runsData[,i,var.i])
	}
	lines(rownames(defRun),medians[,var.i],col='red',lwd=3)
	lines(rownames(defRun),ciBounds,col='blue',lwd=3)
	lines(rownames(defRun),runsData[,which.min(ciLogLikes[q.i,var.i,]),var.i],col='purple',lwd=3)
	plot(rownames(defRun),rownames(defRun),type='n',ylim=c(-3,3))
	for(i in 1:100){
		lines(rownames(defRun),errors[,i])
	}
	lines(rownames(defRun),medians[,var.i]-ciBounds,col='red',lwd=3)
	idc <- which.min(ciLogLikes[q.i,var.i,1:100])
	lines(rownames(defRun),runsData[,idc,var.i]-ciBounds,col='purple',lwd=3)
	
	lines(rownames(defRun),runsData[,sample.i,var.i]-ciBounds,col='green',lwd=4)
}

