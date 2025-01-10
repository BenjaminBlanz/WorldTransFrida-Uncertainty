# change the varExclusion list before running
# you can use the varExclusionListExcludeEverythig... file

# system(paste('rm -r',location.output))
# dir.create(location.output,F,T)
source('initialise.R')
source('config.R')
location.output <- 'testing'
source('runInitialiseData.R')
resSigma <- resDat.cv
# 
# samplePoints.sample <- samplePoints.logLikeSorted[1:1000,]
# rownames(samplePoints.sample) <- 1:nrow(samplePoints.sample)
# if('logLike'%in%colnames(samplePoints)){
# 	samplePoints.sample <- samplePoints.sample[,-which(colnames(samplePoints.sample)=='logLike')]
# }

sampleParms <- readRDS(file.path(location.output,'sampleParmsParscaleRanged.RDS'))

sampleParms.test <- sampleParms
# sampleParms.test$parscale[246] <- sampleParms$parscale[246]/1e5
sampleParms.test$parscale[246] <- 0.01222698
for(i in 1:nrow(sampleParms.test)){
	sampleParms.test$Min <- sampleParms.test$Value-1e-04*abs(sampleParms.test$parscale)
	sampleParms.test$Max <- sampleParms.test$Value+1e-04*abs(sampleParms.test$parscale)
}

samplePoints.sample <- generateSobolSequenceForSampleParms(sampleParms.test,1400,ignoreExistingResults = T)

billysFile <- read.csv('for_ben_symetric.csv')

for(var.i in 1:length(billysFile$Variable)){
	sampleParms.idx <- which(sampleParms$Variable==billysFile$Variable[var.i])
	varmin <- sampleParms$Min[sampleParms.idx]
	varmax <- sampleParms$Max[sampleParms.idx]
	varvalue <- sampleParms$Value[sampleParms.idx]
	maxDistance <- max(varvalue-varmin,varmax-varvalue)
	sampleParms$Max[sampleParms.idx] <- sampleParms$Value[sampleParms.idx]+maxDistance
	sampleParms$Min[sampleParms.idx] <- sampleParms$Value[sampleParms.idx]-maxDistance
}
samplePoints <- generateSobolSequenceForSampleParms(sampleParms,1000,ignoreExistingResults = T)
write.csv(samplePoints,'samples1400ForBilly.csv')
numSample <- nrow(samplePoints.sample)
# write.csv(samplePoints.sample,file.path(location.output,'teensyRangeExperimentSamplePoints.csv'))

# runFRIDASpecParms(samplePoints.sample[1,])
# runFridaDefaultParms()

baseNegLL <- -defLogLike

# frida_info <- read.csv(file.path(location.frida.info,name.frida_info))
# parVect <- frida_info$Value
# names(parVect) <- frida_info$Variable
# negLLike(parVect)
# 
# parVect <- sampleParms$Value
# names(parVect) <- sampleParms$Value
# negLLike(parVect)

# length(samplePoints.sample[1,]) - length(parVect)
# samplePoints.sample[1,] - parVect

source('cleanup.R')
source('clusterHelp.R')
system(paste('rm',file.path(location.output,'workUnit_*')))
redoAllCalc <- T
clusterExport(cl,'redoAllCalc')
logLike.subsample <- clusterRunFridaForSamplePoints(samplePoints.sample,chunkSizePerWorker,
																										 calDat=calDat,
																										 resSigma=resSigma,
																										 location.output=file.path(location.output,'detectedParmSpace'),
																										 redoAllCalc=redoAllCalc,
																										 plotDatWhileRunning=F,
																										 plotDatPerChunWhileRunning=T,
																										 baseLL=-baseNegLL)
hist(logLike.subsample[logLike.subsample>-1e5],xlim=range(logLike.subsample[logLike.subsample>-1e5],-baseNegLL),
		 breaks = 150)
abline(v=defLogLike,col='red',lwd=3)
sqrtNplots <- sqrt(nrow(sampleParms))
plotCols <- round(sqrtNplots)
plotRows <- ceiling(sqrtNplots)
funPlotParRangesLikelihoods(sampleParms.test,samplePoints=samplePoints.sample,like=logLike.subsample,
														baseLike = defLogLike,minLike=defLogLike-abs(0.99*defLogLike),
														outputFilePath = file.path(location.output,'figures'))#,ylim=range(logLike))
if(FALSE){
	# experiments with shifting the likelihood into the double range
range(exp(mpfr(logLike.subsample,32)) * mpfr(10,32)^0)
like.subsample <- exp(mpfr(logLike.subsample,32)) * mpfr(10,32)^0
like.subsample <- as.double(like.subsample)
range(like.subsample)
plot(samplePoints.sample[,1],like.subsample)

png(file.path(location.output,'parameterLikelihoods.png'),
		width=plotCols*10,height=plotRows*10,res=150,units='cm')
funPlotParRangesLikelihoods(sampleParms.test,samplePoints=samplePoints.sample,like=like.subsample,
														baseLike = exp(defLogLike))#,ylim=range(logLike))
dev.off()


max(logLike.subsample)

samplePoints.sample$logLike <- logLike.subsample
}
