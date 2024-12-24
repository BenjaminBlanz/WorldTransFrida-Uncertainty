#
# Manager script for the uncertainty analysis of FRIDA
# 
#
# 2024 Benjamin Blanz
# 
gc()
require(SobolSequence)
require(tictoc)
require(parallel)
require(scales)
source('funRunFRIDA.R')
source('funPlot.R')

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
	calDat.lst <- readRDS(file.path(location.output,'calDat.RDS'))
	calDat <- calDat.lst$calDat
	calDat.impExtrValue <- calDat.lst$calDat.impExtrValue
} else {
	stop('Missing calDat file. Run runFRIDASensitivityRunLikelihood first.\n')
}
# specify sampling parameters ####
# reads frida_info.csv and outputs the SampleParms
# also removes parms we will not sample
# and complains about invalid lines in frida_info.csv
if(file.exists(file.path(location.output,'sampleParms.RDS'))){
	cat('Reading sampling parameters...')
	sampleParms <- readRDS(file.path(location.output,'sampleParms.RDS'))
	cat('done\n')
} else {
	cat('Specify sampling parameters...')
	# read in the parameters in frida that have ranges defined
	frida_info <- read.csv("frida_info.csv")
	columnsThatAreFlags <- c(2,3,4,5,6,7,8,9,10)
	# select the parameters to be sampled
	sampleParms <- frida_info[rowSums(frida_info[,columnsThatAreFlags])>0 &
															frida_info$No.Sensi==0 &
															frida_info$Policy==0,
														-columnsThatAreFlags]
	invalidLines <- which(!((sampleParms$Max-sampleParms$Min)>0 &
														sampleParms$Min <= sampleParms$Value &
														sampleParms$Value <= sampleParms$Max))
	# cat('invalid parm specs min val max\n')
	# for(i in invalidLines){
	# 	cat(sprintf('%5i %-100s %10.f %10.f %10.f\n',
	# 							i,
	# 							sampleParms$Variable[i],
	# 							sampleParms$Min[i],
	# 							sampleParms$Value[i],
	# 							sampleParms$Max[i]))
	# }
	suppressWarnings(file.remove('frida_info_errorCases.csv'))
	if(length(invalidLines)>0){
		cat('invalid lines detected, see frida_info_errorCases.csv...')
	}
	write.csv(sampleParms[invalidLines,],'frida_info_errorCases.csv')
	sampleParms <- sampleParms[-invalidLines,]
	cat('done\n')
}
sampleParms.orig <- sampleParms
saveRDS(sampleParms,file.path(location.output,'sampleParms.RDS'))

# generate sobol sequence ####
if(file.exists(file.path(location.output,'samplePoints.RDS'))){
	cat('Reading sampling points...')
	samplePoints <- readRDS(file.path(location.output,'samplePoints.RDS'))
	samplePoints.base <- readRDS(file.path(location.output,'samplePointsBase.RDS'))
	cat('done\n')
} else {
	cat('Generate sampling points using sobol sequence...')
	# sobolSequence.points generates points on the unit interval for each var
	# transformed, so vars are in rows samples in cols, makes the next steps easier
	samplePoints.base <- sobolSequence.points(nrow(sampleParms),31,numSample) 
	if(sum(duplicated(samplePoints.base))>0){
		stop('Not enough unique sample points. Check the sobol generation\n')
	}
	samplePoints <- funStretchSamplePoints(samplePoints.base,sampleParms,restretchSamplePoints)
	# samplePoints <- rbind(samplePoints, t(sampleParms$Value))
	cat('done\nWriting sample points to file...')
	rownames(samplePoints) <- 1:nrow(samplePoints)
	colnames(samplePoints) <- sampleParms[,1]
	saveRDS(samplePoints,file.path(location.output,'samplePoints.RDS'))
	saveRDS(samplePoints.base,file.path(location.output,'samplePointsBase.RDS'))
	cat('done\n')
}
samplePoints.orig <- samplePoints

# run FRIDA with the samples ####
logLike <- rep(NA,numSample)
like <- rep(NA,numSample)
names(logLike) <- 1:numSample
names(like) <- 1:numSample

## cluster setup ####
cat('cluster setup...')
baseWD <- getwd()
workDirBasename <- 'workDir_'
# start cluster
cl <- makeForkCluster(numWorkers,renice=15)
workers <- 1:length(cl)
# make working directories
gobble <- clusterApply(cl,workers,function(i){
	workerID <- i
	dir.create(file.path('workerDirs',paste0(workDirBasename,i)),showWarnings = F,recursive = T)
	setwd(file.path('workerDirs',paste0(workDirBasename,i)))
	})
# clusterEvalQ(cl,getwd())
# copy over the model and simulator
gobble <- clusterApply(cl,workers,function(i){
	system(paste('cp -r',file.path(baseWD,location.frida),getwd()))})
gobble <- clusterApply(cl,workers,function(i){
	system(paste('cp -r',file.path(baseWD,location.stella),getwd()))})
cat('done\n')


## tighten parms ####
location.output <- file.path(location.output,'BaseParmRange')
dir.create(location.output,showWarnings = F,recursive = T)
doneChangingParms <- F
tight.i <- 0
while(!doneChangingParms){
	## plot setup ####
	if(plotWhileRunning&&plotDatWhileRunning){
		subPlotLocations <- funPlotDat(calDat,calDat.impExtrValue,yaxPad = yaxPad)
	}
	## cluster run ####
	cat('cluster run...\n')
	workUnitBoundaries <- seq(1,numSample,chunkSizePerWorker*numWorkers)
	# in case the chunkSize is not a perfect divisor of the numSample, add numSample as the 
	# final boundary
	if(workUnitBoundaries[length(workUnitBoundaries)]!=numSample){
		workUnitBoundaries <- c(workUnitBoundaries,numSample)
	}
	# add one to the last work unit boundary, as during running we always deduct one from the next boundary
	workUnitBoundaries[length(workUnitBoundaries)] <- numSample+1
	
	### initialise  cluster ####
	cat('  initialising cluster global env...')
	baseWD <- getwd()
	clusterExport(cl,list('location.output','baseWD','sampleParms',
												'chunkSizePerWorker','runFridaParmsBySamplePoints',
												'calDat','resSigma',
												'runFridaParmsByIndex'))
	cat('done\n')
	### running ####
	cat(sprintf('  Run of %i runs split up into %i work units.\n',
							numSample,length(workUnitBoundaries)-1))
	chunkTimes <- c()
	for(i in 1:(length(workUnitBoundaries)-1)){
		if(file.exists(file.path(location.output,paste0('workUnit-',i,'.RDS')))){
			cat(sprintf('\r    Reading existing unit %i',i))
			tryCatch({parOutput <- readRDS(file.path(location.output,paste0('workUnit-',i,'.RDS')))},
							 error = function(e){},warning=function(w){})
		} 
		if(!exists('parOutput')){
			cat(sprintf('\r(r) Running unit %i: samples %i to %i',
									i, workUnitBoundaries[i],workUnitBoundaries[i+1]-1))
			if(length(chunkTimes>1)){
				cat(sprintf(', average duration per unit so far %i sec (%.2f r/s, %.2f r/s/thread), expect completion in %i sec',
										round(mean(chunkTimes,na.rm=T)),
										length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
										chunkSizePerWorker/mean(chunkTimes,na.rm=T),
										round(mean(chunkTimes,na.rm=T))*(length(workUnitBoundaries)-i)))
			}
			tic()
			workUnit <- workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)
			workerWorkUnits <- chunk(workUnit,numWorkers)
			# write the samplePoints of the work units to the worker directories
			for(w.i in workers){
				if(!is.null(workerWorkUnits[[w.i]])){
					saveRDS(samplePoints[workerWorkUnits[[w.i]],],
									file.path('workerDirs',paste0(workDirBasename,w.i),'samplePoints.RDS'))
				} else {
					saveRDS(samplePoints[c(),],
									file.path('workerDirs',paste0(workDirBasename,w.i),'samplePoints.RDS'))
				}
			}
			gobble <- clusterEvalQ(cl,{
				samplePoints <- readRDS('samplePoints.RDS')
			})
			parOutput <- unlist(
				clusterEvalQ(cl,runFridaParmsBySamplePoints()),
				recursive = F)
			timing <- toc(quiet=T)
			chunkTimes[i] <- timing$toc-timing$tic
			cat('\r(s)')
			saveRDS(parOutput,file.path(location.output,paste0('workUnit-',i,'.RDS')))
			cat('\r   ')
		}
		cat('\r(l)')
		for(l in 1:length(parOutput)){
			logLike[parOutput[[l]]$parmsIndex] <- parOutput[[l]]$logLike
			like[parOutput[[l]]$parmsIndex] <- parOutput[[l]]$like
		}
		cat('\r   ')
		if(plotWhileRunning&&plotDatWhileRunning){
			cat('\r(p)')
			for(dat.i in 1:ncol(calDat)){
				par(mfg = which(subPlotLocations==dat.i,arr.ind = T))
				xlims <- c(min(as.numeric(rownames(calDat)))-0.5,
									 max(as.numeric(rownames(calDat)))+0.5)
				yrange <- range(calDat[[dat.i]],na.rm=T)
				ylims <- c(yrange[1]-abs(diff(yrange))*yaxPad,
									 yrange[2]+abs(diff(yrange))*yaxPad)
				plot(rownames(calDat),calDat[[dat.i]],type='n',
						 xaxt='n',yaxt='n',
						 xaxs='i',yaxs='i',
						 xlim=xlims,
						 ylim=ylims)
				for(l in 1:length(parOutput)){
					lines(rownames(parOutput[[l]]$runDat),parOutput[[l]]$runDat[[dat.i]],
								col=i)
								# col=alpha(i,min(0.01,max(1,0.01*parOutput[[l]]$like/.Machine$double.eps))))
				}
				cat('\r   ')
			}
		}
		rm(parOutput)
	}
	if(length(chunkTimes)==0){
		cat('\r    all runs read, no calculation necessary.                                \n')
	} else {
		cat(sprintf('\r    runs completed average chunk time %i sec (%.2f r/s, %.2f r/s/thread), over all run time %i sec %s\n',
								round(mean(chunkTimes,na.rm=T)),
								length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								round(sum(chunkTimes,na.rm=T)),
								'                                                                             '))
	}
	
	### save Likelihoods and sample points ####
	cat(sprintf('  Saving run data...'))
	saveData <- list(sampleParms=sampleParms,samplePoints = samplePoints,logLike=logLike)
	saveRDS(saveData,file.path(location.output,paste0('sensiParmsAndLikes.RDS')))
	cat('done\n')
	
	### save figure ####
	if(plotWhileRunning&&plotDatWhileRunning){
		cat(sprintf('  Saving figure...'))
		dev.print(png,width=5*ncol(subPlotLocations),
							height=5*(nrow(subPlotLocations)-1)+5/4,
							unit='cm',res=150,
							file.path(location.output,'likelihoodWeightedModelRuns.png'))
		cat('done\n')
	}
	
	### calculate probability ####
	# cat('    Calculating sample probabilities...')
	# TODO: weight by spacing!
	# TODO: check likelihoods for infinities and use Rmpfr in that case
	
	# suppressPackageStartupMessages(require(Rmpfr))
	# logLike.mpfr <-  mpfr(logLike,32)
	# like.mpfr <- exp(-logLike.mpfr)
	# if(sum(is.infinite(as.double(like.mpfr)))==0){
	# 	like <- as.double(like.mpfr)
	# }
	# like <- exp(logLike)
	# prob <- as.double(prob.mpfr/sum(prob.mpfr))
	# names(prob) <- names(negLogLike)
	# rm(negLogLike.mpfr,prob.mpfr)
	# cat('done\n')
	
	### tighten parms ####
	cat('checking parm bound tightness...\n')
	fullTermLike <- nrow(calDat)*.Machine$double.eps
	likeThreshold <- max(fullTermLike,-fullTermLike+max(like)/likeThresholdRatio)
	if(sum(like>=likeThreshold)==0){
		doneChangingParms <- T
		cat(' seem tight enough\n')
	} else {
		relTightening <- c()
		for(i in 1:nrow(sampleParms)){
			minmax <- range(samplePoints[like>=likeThreshold,i])
			relTightening[i] <- 1-diff(minmax)/(sampleParms[i,c('Max')]-sampleParms[i,c('Min')])
			sampleParms$Min[i] <- minmax[1]
			sampleParms$Max[i] <- minmax[2]
		}
		if(max(relTightening)>0){
			doneChangingParms <- F
			tight.i <- tight.i + 1
			cat(sprintf(' Tightening parm bounds: %i\n',tight.i))
			if(plotWhileRunning){
				funPlotParRangesLikelihoods(sampleParms,sampleParms.orig,samplePoints,like,
																		savePlotFilePath = file.path(location.output,paste0('parmLikelihoods-tightening-',tight.i,'.png')))
			}
			samplePoints <- funStretchSamplePoints(samplePoints.base,sampleParms,restretchSamplePoints)
			location.output <- file.path(location.output,'..',paste0('tightening-',tight.i))
			dir.create(location.output,recursive = T,showWarnings = F)
		} else {
			doneChangingParms <- T
			cat(' seem tight enough\n')
		}
		if(doneChangingParms){
			### plot parm Range ####
			if(plotWhileRunning){
				cat('  Plotting parm range and likelihoods...')
				funPlotParRangesLikelihoods(sampleParms,sampleParms.orig,
																		samplePoints,like,
																		savePlotFilePath = file.path(location.output,'parmLikelihoods-baseParmRange.png'))
				cat('done\n')
			}
		}
	}
	
	### expand parms ####
	#TODO: expand parms
	# if(like does not decrease below maxLike/likeThresholdRatio at the edges
}
## cluster cleanup ####
cat('cluster cleanup...')
# stop cluster
stopCluster(cl)
# clean up working directories
system('rm -r workerDirs')
cat('done\n')


