#
# Manager script for the uncertainty analysis of FRIDA
# 
#
# 2024 Benjamin Blanz
# 
source('initialise.R')

# config ####
cat('Config...')
source('config.R')
location.output.base <- location.output
cat('done\n')

# read input files for likelihood ####
if(file.exists(file.path(location.output,'sigma.RDS'))){
	resSigma <- readRDS(file.path(location.output,'sigma.RDS'))
} else {
	stop('Missing covariance matrix file. Run runInitialiseData.R first.\n')
}
if(file.exists(file.path(location.output,'calDat.RDS'))){
	calDat.lst <- readRDS(file.path(location.output,'calDat.RDS'))
	calDat <- calDat.lst$calDat
	calDat.impExtrValue <- calDat.lst$calDat.impExtrValue
} else {
	stop('Missing calDat file. Run runInitialiseData.R first.\n')
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
	frida_info <- read.csv(file.path(location.frida.info,name.frida_info))
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
samplePoints <- generateSobolSequenceForSampleParms(sampleParms,numSample,
																										restretchSamplePoints,
																										redoAllCalc)
samplePoints.orig <- samplePoints

# run FRIDA with the samples ####


## cluster setup ####
source('clusterHelp.R')

## tighten parms ####
location.output <- file.path(location.output.base,'BaseParmRange')
dir.create(location.output,showWarnings = F,recursive = T)
doneChangingParms <- F
tight.i <- 0
while(!doneChangingParms){
	## cluster run ####
	clusterRunRetList <- clusterRunFridaForSamplePoints(samplePoints,chunkSizePerWorker,
																											calDat=calDat,
																											resSigma=resSigma,
																											location.output,
																											redoAllCalc,
																											plotDatWhileRunning)
	logLike <- clusterRunRetList$logLike
	like <- clusterRunRetList$like
	
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
	
	stop()
	
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
	if(!tightenParms){
		doneChangingParms <- T
	} else {
		cat('checking parm bound tightness...\n')
		fullTermLike <- nrow(calDat)*.Machine$double.eps
		likeThreshold <- max(fullTermLike,-fullTermLike+max(like)/likeThresholdRatio)
		if(sum(like>=likeThreshold)==0){
			doneChangingParms <- T
			cat(' seem tight enough\n') 
		}	else if (sum(like>=likeThreshold)<=nrow(sampleParms)){
			cat(sprintf('Not enough samples above threshold.\n %i above threshold, require %i.\n',
					sum(like>=likeThreshold),
					nrow(sampleParms)))
			doneChangingParms <- T
		} else {
			cat(sprintf(' Tightening parm bounds: %i\n',tight.i))
			relTightening <- c()
			for(i in 1:nrow(sampleParms)){
				minmax <- range(samplePoints[like>=likeThreshold,i])
				relTightening[i] <- 1-diff(minmax)/(sampleParms[i,c('Max')]-sampleParms[i,c('Min')])
				sampleParms$Min[i] <- minmax[1]
				sampleParms$Max[i] <- minmax[2]
			}
			if(max(relTightening)>0 && max(relTightening)<1){
				doneChangingParms <- F
				tight.i <- tight.i + 1
				cat(' ...done\n')
				if(plotWhileRunning){
					cat('   plotting...')
					funPlotParRangesLikelihoods(sampleParms,sampleParms.orig,samplePoints,like,
																			savePlotFilePath = file.path(location.output,paste0('parmLikelihoods-tightening-',tight.i,'.png')))
					cat('done\n')
				}
				cat('   stretching sample points accross new parm space...')
				samplePoints <- funStretchSamplePoints(samplePoints.base,sampleParms,restretchSamplePoints)
				cat('done\n')
				location.output <- file.path(location.output.base,paste0('tightening-',tight.i))
				dir.create(location.output,recursive = T,showWarnings = F)
			} else if (sum(relTightening==1)<=nrow(sampleParms)){
				cat('Tightening singularity. Can not continue\n')
				doneChangingParms <- T
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
	}
	
	### expand parms ####
	#TODO: expand parms
	# if(like does not decrease below maxLike/likeThresholdRatio at the edges
}


