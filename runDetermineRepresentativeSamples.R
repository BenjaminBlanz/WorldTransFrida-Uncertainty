# The likelihood that a given run is "close" to the desired quantile of all runs in all variable
# dimensions is the likelihood that the difference between this run and the desired quantile is
# the same as the difference between the desired quantile and other runs in gerneral.
# 
# How to account for quantiles different from 0.5 where the differences are assymmetric?
# The likelihood for the quantiles is the likelihood that the errors are from 
# a distribution with the variance at median and a mean of the difference between the
# desired quantile and the median

cat('\nRunning runDetermineRepresentativeSample.R\n\n')

library('parallel')
source('initialise.R')
source('config.R')

if(exists('preexsistingBaselineFolder')&&!is.na(preexsistingBaselineFolder)){
	location.output <- preexsistingBaselineFolder
	cat(sprintf('Overriding output location:\n%s\n',location.output))
}

source('runInitialiseData.R')
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
location.runFiles <- file.path(location.output,'detectedParmSpace',paste0('PerVarFiles-',perVarOutputTypes[1]))

# collect time series ####
cat('reading vars from run data...\n')
defRun <- runFridaDefaultParms()
yearsToRead <- rownames(defRun)
# dimensions time in rows, run IDs in columns,  variables to read
varsToRead <- subSample.TargetVars
runsData <- array(NA,dim=c(nrow(defRun),nrow(samplePoints),length(varsToRead)))
dimnames(runsData) <- list(rownames(defRun),1:nrow(samplePoints),varsToRead)
for(f.i in 1:length(varsToRead)){
	cat(sprintf('reading file %i of %i: %s...',f.i,length(varsToRead),varsToRead[f.i]))
	perVarData <- readPerVarFile(file.path(location.runFiles,varsToRead[f.i]),outputType = perVarOutputTypes[1])
	if(nrow(perVarData)!=nrow(samplePoints)||(ncol(perVarData)-1)!=nrow(defRun)){
		stop(sprintf(paste0('%s holds %i runs over %i years, expected %i over %i.\n',
												'Rerun the ensemble for this variable.\n'),
								 varsToRead[f.i],nrow(perVarData),ncol(perVarData)-1,
								 nrow(samplePoints),nrow(defRun)))
	}
	runsData[,,f.i] <- t(perVarData[,-1])
	cat('done\n')
}

# prep quantiles ####
# determine the values of the desired quantiles in all of the variables
for(plotWeightType in plotWeightTypes){
	cat('determining weights...\n')
	if(plotWeightType %in% c('likelihood','logCutoff','linearly','completeEqually')){
		# log like ####
		cat(' reading log likelihoods...\n')
		logLike.perVar <- readPerVarFile(file.path(location.runFiles,'logLike'),outputType = perVarOutputTypes[1])
		logLike <- rep(NA,nrow(samplePoints))
		logLike[logLike.perVar$id] <- logLike.perVar$logLike
		rm(logLike.perVar)
		# whether a run completed comes from the run status rather than from the log
		# likelihood markers, indexed by id so that a run missing from the file does
		# not shift everything after it
		runStatus <- funReadRunStatus(location.output,outputType = perVarOutputTypes[1],
																	numSample = nrow(samplePoints))
		completed <- rep(NA,nrow(samplePoints))
		completed[runStatus$id] <- runStatus$completed
		likelihoodOK <- rep(NA,nrow(samplePoints))
		likelihoodOK[runStatus$id] <- runStatus$likelihoodOK
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
	} else if(plotWeightType == 'completeEqually'){
		# equal weighting of completed runs. A run counts when it reached the final
		# year and its log likelihood is a real value. likelihoodOK is NA where there
		# is no calibration likelihood at all, and that must not zero every weight.
		samplePoints$plotWeight <- 0
		samplePoints$plotWeight[completed%in%1 & !(likelihoodOK%in%0)] <- 1
	} else if(plotWeightType == 'linearly'){
		samplePoints$plotWeight <- order(logLike)/nrow(samplePoints)
	} else {
		stop('unknown plotWeightType\n'	)
	}
	# With no weight anywhere every weighted quantile below is undefined, and the
	# search for the sample points closest to them dies inside the mini cluster
	# on a which.min of nothing, reported as 'replacement has length zero' with
	# no hint of what went wrong. Say what went wrong here instead.
	if(!any(samplePoints$plotWeight>0,na.rm=TRUE)){
		stop(sprintf(paste0('No sample point carries any weight for plotWeightType %s, ',
												'so there is no ensemble to draw a representative sample from.\n',
												'%i of %i runs completed. Check the parameter ranges the\n',
												'sampling drew from and the run completion summary in runMetadata.txt.\n'),
								 plotWeightType,
								 if(exists('completed')){sum(completed%in%1)}else{NA},
								 nrow(samplePoints)))
	}
	# median ####
	medians <- array(NA,dim=c(nrow(defRun),length(varsToRead)))
	for(var.i in 1:length(varsToRead)){
		for(year.i in 1:nrow(defRun)){
			medians[year.i,var.i] <- weighted.quantile(runsData[year.i,,var.i],
																								 w = samplePoints$plotWeight,
																								 probs = 0.5,
																								 na.rm = T)
		}
	}
	cat('...weights determined\n')

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

	# quantile min SSEs ###
	cat('determining samples that closest match the desired quantiles...')
	parMinSSEFun <- function(p.i){
		minSSEidc <- rep(NA,length(subSample.TargetVars))
		for(var.i in 1:length(varsToRead)){
			ciBounds <- rep(NA,nrow(defRun))
			errors <- array(NA,dim=c(nrow(defRun),nrow(samplePoints)))
			SSEs <- rep(NA,nrow(samplePoints))
			for(year.i in 1:nrow(defRun)){
				ciBounds[year.i] <- weighted.quantile(runsData[year.i,,var.i],
																							w = samplePoints$plotWeight,
																							probs = subSample.Ps[p.i],
																							na.rm = T)
				errors[year.i,] <- runsData[year.i,,var.i] - ciBounds[year.i]
			}
			for(sample.i in 1:nrow(samplePoints)){
				SSEs[sample.i] <- sum(samplePoints$plotWeight[sample.i] *
																(errors[,sample.i])^2)
			}
			minSSEidc[var.i] <- which.min(SSEs)
		}
		return(minSSEidc)
	}
	minicl <- makeForkCluster(3)
	minSSEidc.lst <- parLapply(minicl,1:subSample.NumSamplePerVar,parMinSSEFun)
	stopCluster(minicl)
	minSSEidc <- array(NA,dim=c(subSample.NumSamplePerVar,length(subSample.TargetVars)))
	for(p.i in 1:subSample.NumSamplePerVar){
		minSSEidc[p.i,] <- minSSEidc.lst[[p.i]]
	}
	cat('\n    ')

	repSample <- samplePoints[as.vector(minSSEidc),]
	colnames(repSample) <- gsub('\\[1\\]','',colnames(repSample))
	colnames(repSample) <- gsub('\\[1,','[*,',colnames(repSample))
	# plotWeight and logLike are bookkeeping this script put on samplePoints, not
	# model parameters. Both have to go before the file is written: FRIDA imports it
	# as a parameter file (see policy_update_frida.sh), so a column left in here
	# reaches the model looking like a parameter. Selecting by name rather than by
	# -which() also avoids dropping every column when neither is present.
	repSample <- repSample[,!colnames(repSample)%in%c('plotWeight','logLike')]
	cat('done\n')
	cat('writing out to subSampleParameterValues.csv ...')
	location.output.repSample <- file.path(location.output,'repSample',plotWeightType)
	dir.create(location.output.repSample,F,T)
	write.table(repSample,file.path(location.output.repSample,'subSampleParameterValues.csv'),
							append = F,sep = ',',row.names = F)
	write.table(data.frame(id=minSSEidc),file.path(location.output.repSample,'subSampleParameterIndices.csv'),
							append = F,sep = ',',row.names = F)
	saveRDS(repSample,file.path(location.output.repSample,'subSampleParameterValues.RDS'))
	saveRDS(data.frame(id=minSSEidc),file.path(location.output.repSample,'subSampleParameterIndices.RDS'))
	cat('done\n')
	
	# plot ####
	cat('plotting selected samples...')
	sqrtNumPlots <- sqrt(length(subSample.TargetVars))
	plotCols <- round(sqrtNumPlots)
	plotRows <- ceiling(sqrtNumPlots)
	par(mfrow=c(plotRows,plotCols))
	subsampleFigDir <- file.path(location.output,'figures','subSample.TargetVars')
	dir.create(file.path(subsampleFigDir),showWarnings = F,recursive = T)
	for(var.i in 1:length(subSample.TargetVars)){
		png(file.path(subsampleFigDir,paste0(subSample.TargetVars[var.i],'.png')),
				width = plotWidth, height = plotHeight, units = plotUnit,res = plotRes)
		plot(rownames(defRun),defRun[,subSample.TargetVars[var.i]],ylim=range(runsData[,minSSEidc[,var.i],var.i]),
				 type='n',xlab='year',ylab=subSample.TargetVars[var.i])
		for(var.ii in 1:length(subSample.TargetVars)){
			for(p.i in 1:subSample.NumSamplePerVar){
				lines(rownames(defRun),runsData[,minSSEidc[p.i,var.ii],var.i],
							col=adjustcolor(hsv(var.ii/length(subSample.TargetVars),1,0.25+p.i/(subSample.NumSamplePerVar*3)),
															alpha.f=0.5),
							lwd=3)
			}
		}
		dev.off()
	}
	cat('done\n')
}

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
	lines(rownames(defRun),runsData[,minSSEidc[p.i,var.i],var.i],col='purple',lwd=3)
	
	plot(rownames(defRun),rownames(defRun),type='n',ylim=c(-3,3))
	for(i in 1:100){
		lines(rownames(defRun),errors[,i])
	}
	abline(h=0,col='red',lwd=3)
	# lines(rownames(defRun),medians[,var.i]-ciBounds,col='red',lwd=3)
	lines(rownames(defRun),runsData[,minSSEidc[p.i,var.i],var.i]-ciBounds,col='purple',lwd=3)
}
