source('initialise.R')
source('config.R')
location.output <- file.path('workOutput','testCases')
source('runInitialiseData.R')
calDat <- readRDS(file.path(location.output,'calDat.RDS'))$calDat
resSigma <- readRDS(file.path(location.output,'sigma-indepParms.RDS'))
location.output.base <- location.output
dir.create(file.path(location.output),recursive = T,showWarnings = F)
sampleParms <- prepareSampleParms()
parVect.default <- sampleParms$Value
names(parVect.default) <- sampleParms$Variable

# test 1 ####
thingToMod <- "Finance.sensitivity of effect of sta on failure rate[1]"
idxToMod <- which(sampleParms$Variable==thingToMod)
parVect.highFailure <- parVect.default
parVect.highFailure[thingToMod] <- 0.85
parVect.lowFailure <- parVect.default
parVect.lowFailure[thingToMod] <- 0.3

samplePoints <- rbind(parVect.default,parVect.highFailure,parVect.lowFailure)
# colnames(samplePoints) <-
runDats <- runFridaParmsBySamplePoints()

runDat.default <- runDats[[1]]$runDat
runDat.highFailure <- runDats[[2]]$runDat
runDat.lowFailure <- runDats[[3]]$runDat

plot(rownames(calDat),calDat$gdp_nominal_gdp,
		 xlab='year',ylab='gdp_nominal_gdp',
		 ylim=c(0,205000),pch='.',cex=3)
lines(rownames(calDat),runDat.default$gdp_nominal_gdp,col='red',lwd=2)
lines(rownames(calDat),runDat.lowFailure$gdp_nominal_gdp,col='pink',lwd=3)
lines(rownames(calDat),runDat.highFailure$gdp_nominal_gdp,col='darkgreen',lwd=2)

# runDat.noParm <- runFRIDASpecParms(c())

# test 2 ####
# test case 2 is doThis.csv
# run as cluster with runsPerWorkerperChung 10
# put parameters in the export spec to verify 

# do this
doThis <- read.csv(file.path(location.output,'Do This.csv'),row.names = 1,header = F)
parNames <- unlist(doThis[1,])
doThis <- read.csv(file.path(location.output,'Do This.csv'),row.names = 1)
colnames(doThis) <- parNames
samplePoints <- doThis


source('clusterHelp.R')

# write firda export vars
# doing this after clusterHelp, i.e. after tmpfs means this is not persistent after cleanup
varsForExport.fridaNames <- c(varsForExport.fridaNames,parNames)
writeFRIDAExportSpec(varsForExport.fridaNames,location.frida)

## run single ####
tic()
runDat.singleThread <- runFridaParmsBySamplePoints()
timing <- toc(quiet = T)
cat(sprintf('%.2f r/s, %i runs in %s',
						(timing$toc-timing$tic)/nrow(samplePoints),
						nrow(samplePoints),dseconds(round(timing$toc-timing$tic,digits=1))))

### plot ####
baseLL <- runDats[[1]]$logLike
dat.is <- c(1,132)
par(mfrow=c(1,2))
for(dat.i in dat.is){
	yrange <- c(NA,NA)
	for(l in 1:nrow(samplePoints)){
		yrange[1] <- min(yrange[1],runDat.singleThread[[l]]$runDat[[dat.i]],na.rm=T)
		yrange[2] <- max(yrange[2],runDat.singleThread[[l]]$runDat[[dat.i]],na.rm=T)
	}
	plot(rownames(calDat),calDat[[dat.i]],
			 xlab='year',ylab=colnames(calDat)[dat.i],
			 ylim=yrange,pch='.',cex=3)
	for(l in 1:nrow(samplePoints)){
		runLL <- runDat.singleThread[[l]]$logLike
		lines(rownames(runDat.singleThread[[l]]$runDat),runDat.singleThread[[l]]$runDat[[dat.i]],
					col=1
					# col=adjustcolor(1,min(1,max(0.01,
					# 														100/(abs(runLL-baseLL)+1))))
					)
	}
	points(rownames(calDat),calDat[[dat.i]],col='red',pch=20)
}

### write to csv ####
test2outputPath <- file.path(location.output,'test2','singleThreaded')
dir.create(test2outputPath,F,T)
for(l in 1:nrow(samplePoints)){
	write.csv(runDat.singleThread[[l]]$runDat,file.path(test2outputPath,paste0(rownames(samplePoints)[l],'.csv')))
}

## run cluster ####
for(worker in 1:length(cl)){
	writeFRIDAExportSpec(varsForExport.fridaNames,file.path('workerDirs',paste0('workDir_',worker),'FRIDAforUncertaintyAnalysis'))
}
chunkSizePerWorker <- 5
test2outputPathMT <- file.path(location.output,'test2','multiThreaded')
dir.create(test2outputPathMT,F,T)
logLike.multiThreaded <- clusterRunFridaForSamplePoints(samplePoints,chunkSizePerWorker,
																												location.output=test2outputPathMT,
																												calDat,resSigma,redoAllCalc = T,
																												plotDatWhileRunning = T,
																												plotDatPerChunWhileRunning = T,
																												plotPerChunk = T)
runDat.multiThreaded <- loadClusterRuns(test2outputPathMT)

### plot ####
baseLL <- runDats[[1]]$logLike
dat.is <- c(1,132)
par(mfrow=c(1,2))
for(dat.i in dat.is){
	yrange <- c(NA,NA)
	for(l in 1:nrow(samplePoints)){
		yrange[1] <- min(yrange[1],runDat.multiThreaded[[l]]$runDat[[dat.i]],na.rm=T)
		yrange[2] <- max(yrange[2],runDat.multiThreaded[[l]]$runDat[[dat.i]],na.rm=T)
	}
	plot(rownames(calDat),calDat[[dat.i]],
			 xlab='year',ylab=colnames(calDat)[dat.i],
			 ylim=yrange,pch='.',cex=3)
	for(l in 1:nrow(samplePoints)){
		runLL <- runDat.multiThreaded[[l]]$logLike
		lines(rownames(runDat.multiThreaded[[l]]$runDat),runDat.multiThreaded[[l]]$runDat[[dat.i]],
					# col=1
					col=adjustcolor(1,min(1,max(0.01,
																			100/(abs(runLL-baseLL)+1))))
		)
	}
	points(rownames(calDat),calDat[[dat.i]],col='red',pch=20)
}

### write to csv ####
dir.create(test2outputPathMT,F,T)
for(l in 1:nrow(samplePoints)){
	write.csv(runDat.multiThreaded[[l]]$runDat,file.path(test2outputPathMT,paste0(rownames(samplePoints)[l],'.csv')))
}

















