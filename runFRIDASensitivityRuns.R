#
# Manager script for the uncertainty analysis of FRIDA
# 
#
# 2024 Benjamin Blanz
# 

require(SobolSequence)
require(tictoc)
require(parallel)
require(scales)
source('funRunFRIDA.R')

# config ####
cat('Config...')
source('config.R')
dir.create(location.output,showWarnings=F,recursive=T)
file.copy('config.R',location.output)
cat('done\n')

# read input files for likelihood ####
if(file.exists(file.path(location.output,'sigma.RDS'))){
	resSigma <- readRDS(file.path(location.output,'sigma.RDS'))
} else {
	stop('Missing covariance matrix file. Run runFRIDASensitivityRunLikelihood first.\n')
}
if(file.exists(file.path(location.output,'calDat.RDS'))){
	calDat <- readRDS(file.path(location.output,'calDat.RDS'))
} else {
	stop('Missing calDat file. Run runFRIDASensitivityRunLikelihood first.\n')
}
# specify sampling parameters ####
cat('Specify sampling parameter...')
# read in the parameters in frida that have ranges defined
frida_info <- read.csv("frida_info.csv")
columnsThatAreFlags <- c(2,3,4,5,6,7,8,9,10)
# select the parameters to be sampled
sampleParms <- frida_info[rowSums(frida_info[,columnsThatAreFlags])>0 &
														frida_info$No.Sensi==0 & 
														frida_info$Policy==0,
													-columnsThatAreFlags]
cat('done\n')

# generate sobol sequence ####
cat('Generate sobol sequence...')
# sobolSequence.points generates points on the unit interval for each var
# transformed, so vars are in rows samples in cols, makes the next steps easier
samplePoints <- sobolSequence.points(nrow(sampleParms),31,numSample)
if(sum(duplicated(samplePoints))>0){
	stop('Not enough unique sample points. Check the sobol generation\n')
}
samplePoints <- t(samplePoints)

if(!restretchSamplePoints){
	# Substract the min and multiply by max-min to strecth the unit interval to the 
	# actual sampling range.
	samplePointsStretched <- samplePoints*(sampleParms$Max-sampleParms$Min) + sampleParms$Min
	# plot(samplePointsStretched[1,],samplePointsStretched[2,])
	# abline(v=sampleParms$Value[1],h=sampleParms$Value[2],col='red')
	samplePoints <- samplePointsStretched
	rm(samplePointsStretched)
} else {
	# stretch the sample points to be left and right of the mean centre value of the 
	# description file
	lowIdc <- samplePoints<0.5
	highIdc <- samplePoints>=0.5
	samplePointsLow <- samplePointsHigh <- samplePoints
	samplePointsLow[highIdc] <- NA
	samplePointsLow <- samplePointsLow*2*(sampleParms$Value-sampleParms$Min) + sampleParms$Min
	samplePointsHigh[lowIdc] <- NA
	samplePointsHigh <- (samplePoints-0.5)*2*(sampleParms$Max-sampleParms$Value) + sampleParms$Value
	samplePointsReStretched <- samplePoints
	samplePointsReStretched[lowIdc] <- samplePointsLow[lowIdc]
	samplePointsReStretched[highIdc] <- samplePointsHigh[highIdc]
	# plot(samplePointsReStretched[1,],samplePointsReStretched[2,])
	# abline(v=sampleParms$Value[1],h=sampleParms$Value[2],col='red')
	samplePoints <- samplePointsReStretched
	rm(samplePointsHigh,samplePointsLow,samplePointsReStretched,lowIdc,highIdc)
}
cat('done\n')

# run FRIDA with the samples ####
negLogLikes <- data.frame(parmsIndex = 1:numSample,
													negLogLike=rep(0,numSample))
## cluster ####
### cluster setup ####
cat('cluster setup...')
baseWD <- getwd()
workDirBasename <- 'workDir_'
# start cluster
cl <- makeForkCluster(numWorkers,renice=15)
workers <- 1:length(cl)
# make working directories
gobble <- clusterApply(cl,workers,function(i){
	dir.create(file.path('workerDirs',paste0(workDirBasename,i)),showWarnings = F,recursive = T)}) 
gobble <- clusterApply(cl,workers,function(i){
	setwd(file.path('workerDirs',paste0(workDirBasename,i)))})
#clusterEvalQ(cl,getwd())
# copy over the model and simulator
gobble <- clusterApply(cl,workers,function(i){
	system(paste('cp -r',file.path(baseWD,location.frida),getwd()))})
gobble <- clusterApply(cl,workers,function(i){
	system(paste('cp -r',file.path(baseWD,location.stella),getwd()))})
cat('done\n')

### cluster run ####
cat('cluster run...')
if(plotWhileRunning){
	defDat <- runFridaDefaultParms()
	par(mfrow=c(1,1))
	par(mar=c(5.1,5.1,4.1,2.1))
	ylims <- range(defDat[[whatToPlot]])
	plot(rownames(defDat),defDat[[whatToPlot]],type='l',
			 ylim=c(ylims[1]-diff(ylims)*0.2,ylims[2]+diff(ylims)*0.2),xlim=c(1980,2130),
			 xlab='year',ylab='Real GDP in 2021 bn intl$',
			 xaxt='n',lwd=4)
	axis(1,at=seq(1980,2130,10))
}
workUnitBoundaries <- seq(1,numSample,chunkSizePerWorker*length(cl))
# in case the chunkSize is not a perfect divisor of the numSample, add numSample as the 
# final boundary
if(workUnitBoundaries[length(workUnitBoundaries)]!=numSample){
	workUnitBoundaries <- c(workUnitBoundaries,numSample)
}
# add one to the last work unit boundary, as during running we always deduct one from the next boundary
workUnitBoundaries[length(workUnitBoundaries)] <- numSample+1
cat(sprintf('  Run of %i runs split up into %i work units.\n',
						numSample,length(workUnitBoundaries)-1))
chunkTimes <- c()
for(i in 1:(length(workUnitBoundaries)-1)){
	if(file.exists(file.path(location.output,paste0('workUnit-',i,'.RDS')))){
		cat(sprintf('\r   Reading existing unit %i',i))
		tryCatch({parOutput <- readRDS(file.path(location.output,paste0('workUnit-',i,'.RDS')))},
						 error = function(e){},warning=function(w){})
	} 
	if(!exists('parOutput')){
		cat(sprintf('\r   Running unit %i: samples %i to %i',
								i, workUnitBoundaries[i],workUnitBoundaries[i+1]-1))
		if(length(chunkTimes>1)){
			cat(sprintf(', average duration per unit so far %i sec, expect completion in %i sec',
									round(mean(chunkTimes,na.rm=T)),round(mean(chunkTimes,na.rm=T))*(length(workUnitBoundaries)-i)))
		}
		tic()
		workUnit <- workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)
		parOutput <- parLapplyLB(cl=cl,X=workUnit,fun = runFridaParmsByIndex)
		timing <- toc(quiet=T)
		chunkTimes[i] <- timing$toc-timing$tic
		saveRDS(parOutput,file.path(location.output,paste0('workUnit-',i,'.RDS')))
	}
	for(l in 1:length(parOutput)){
		negLogLikes$negLogLike[parOutput[[l]]$parmsIndex] <- parOutput[[l]]$negLogLike
	}
	if(plotWhileRunning){
		# readline(prompt="Press [enter] to continue")
		for(l in 1:length(parOutput)){
			lines(rownames(parOutput[[l]]$runDat),parOutput[[l]]$runDat[[whatToPlot]],
						col=alpha(i,min(0.25,max(0.05,parOutput[[l]]$negLogLike/1e7))))
		}
	}
	rm(parOutput)
}
if(length(chunkTimes)==0){
	cat('\r    all runs read, no calculation necessary.                                \n')
} else {
	cat(sprintf('\r    runs completed average chunk time %i sec, over all run time %i sec                                       \n',
							round(mean(chunkTimes,na.rm=T)),round(sum(chunkTimes,na.rm=T))))
}
cat('done\n')

### save Likelihoods and sample points ####
saveData <- list(sampleParms=sampleParms,samplePoints = samplePoints,negLogLikes=negLogLikes)
saveRDS(saveData,file.path(location.output,paste0('sensiParmsAndLikes.RDS')))

### save figure ####
if(plotWhileRunning){
	dev.print(pdf,
						file.path(location.output,
											paste0(whatToPlot,'.pdf')))
}

### cluster cleanup ####
cat('cluster cleanup...')
# stop cluster
stopCluster(cl)
# clean up working directories
system('rm -r workerDirs')
cat('done\n')


