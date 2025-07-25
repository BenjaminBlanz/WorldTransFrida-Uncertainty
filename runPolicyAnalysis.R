# initialise ####
cat('Initialise...')
source('initialise.R')

# config ####
cat('done\nConfig...')
baseWD <- getwd()
source('configPolicyAnalysis.R')
# hardcoding this here, as the below code relies on it.
# this should become standard also in uncertainty in the future
writePerWorkerFiles <- TRUE
doNotReturnRunDataSavePerWorkerOnly <- TRUE
extraVarsToExport <- unique(read.csv(file.path(location.frida.info,name.frida_extra_variables_to_export_list))$FRIDA.FQN)
extraVarsToExport <- extraVarsToExport[nchar(extraVarsToExport)>4]
writeFRIDAExportSpec(extraVarsToExport,location.frida)
cat('done\n')
sink(file.path(location.output,'log.txt'),append=T)
cat(paste0(
	'\n###############################################################\n',
	format(Sys.Date(), "%c"),
	location.output,
	'\n###############################################################\n'))
sink()
sink(file.path(location.output,'log.txt'),append=T,split=T)

# read domain policies ####
cat('Read policies...')
pdp.lst <- list()
for(singleDomainPolicyFile in list.files(location.singleDomainPolicyFiles)){
	pdp.lst[[tools::file_path_sans_ext(singleDomainPolicyFile)]] <-
		read.csv(file.path(location.singleDomainPolicyFiles,singleDomainPolicyFile))
}
cat('done\nGenerate metadata...')
subDomainNames <- names(pdp.lst)
domainNames <-c()
for(i in 1:length(subDomainNames)){
	nameParts <- strsplit(subDomainNames[i],'\\+')[[1]]
	domainNames[i] <- nameParts[1]
}
domainNames <- unique(domainNames)
numDom <- length(domainNames)
nameParts <- strsplit(subDomainNames[1],'\\+')[[1]]
pdpMeta <- data.frame(pfID=1,
											domID=1,
											sdmID=1,
											domain=nameParts[1],
											subdomain=nameParts[2],
											numSDP=max(pdp.lst[[1]]$polID),
											numDP=0)
for(i in 2:length(subDomainNames)){
	nameParts <- strsplit(subDomainNames[i],'\\+')[[1]]
	domID <- pdpMeta$domID[which(pdpMeta$domain==nameParts[1])[1]]
	if(is.na(domID)){
		domID <- which.max(pdpMeta$domID)+1
	}
	sdmID <- max(c(-1,pdpMeta$sdmID[pdpMeta$domID==domID]))
	if(sdmID==-1){
		sdmID <- 1
	} else {
		sdmID <- sdmID + 1
	}
	pdpMeta[nrow(pdpMeta)+1,] <- 
		data.frame(pfID=i,
							 domID=domID,
							 sdmID=sdmID,
							 domain=nameParts[1],subdomain=nameParts[2],
							 numSDP=max(pdp.lst[[i]]$polID),
							 numDP=0)
}
for(domID in unique(pdpMeta$domID)){
	pdpMeta$numDP[pdpMeta$domID==domID] <- sum(pdpMeta$numSDP[pdpMeta$domID==domID])
}
write.csv(pdpMeta,file.path(location.output,'perDomainPolicyMetadata.csv'))
cat('done\n')

# setup policies ####
cat('setting up policies...')
jointPolicies <- data.frame(polID=numeric(),
														domID=numeric(),
														dplID=numeric(),
														sdmID=numeric(),
														sdpID=numeric())
maxPolID <- 0
for(pfID in pdpMeta$pfID){
	jointPolicies <- rbind(jointPolicies, 
												 data.frame(
												 	polID=(maxPolID+1):(maxPolID+pdpMeta$numSDP[pfID]),
												 	dplID=rep(0,pdpMeta$numSDP[pfID]),
												 	domID=rep(pdpMeta$domID[pfID],pdpMeta$numSDP[pfID]),
												 	sdmID=rep(pdpMeta$sdmID[pfID],pdpMeta$numSDP[pfID]),
												 	sdpID=1:pdpMeta$numSDP[pfID])
	)
	maxPolID <- maxPolID+pdpMeta$numSDP[pfID]
}
write.csv(jointPolicies,file.path(location.output,'policyIdentifierList.csv'))

sampleParms <- data.frame(Variable=pdpMeta$domID,
													Value=rep(0,nrow(pdpMeta)),
													Min=rep(1,nrow(pdpMeta)),
													Max=pdpMeta$numDP)
sampleParms <- sampleParms[!duplicated(sampleParms),]
write.csv(sampleParms, file.path(location.output,'policyIndexBounds.csv'))

singlePolicyArr <- array(NA,dim=c(sum(pdpMeta$numSDP),
																 numDom))
colnames(singlePolicyArr) <- as.character(sampleParms$Variable)
jPA.idx <- 0
for(domID in unique(pdpMeta$domID)){
	domCases <- which(jointPolicies$domID==domID)
	jointPolicies$dplID[domCases] <- 1:length(domCases)
	singlePolicyArr[(jPA.idx+1):(jPA.idx+length(domCases)),as.character(domID)] <- jointPolicies$dplID[domCases]
	jPA.idx <- jPA.idx+length(domCases)
}
cat('done\n')
# sample points generation is automatically verbose
samplePoints <- generateSobolSequenceForSampleParms(sampleParms,
																										numSample = numInitialJointPol,
																										restretchSamplePoints = F,
																										ignoreExistingResults = T,
																										integerParms = sampleParms,
																										nullProb = nullPolProb)
samplePoints <- rbind(singlePolicyArr,samplePoints)
numSample <- nrow(samplePoints)
rownames(samplePoints) <- 1:numSample
write.csv(samplePoints,file.path(location.output,'policySamplePoints.csv'))
# hist(rowSums(!is.na(samplePoints)))

# cluster setup ####
source('clusterHelp.R')
clusterExport(cl,list('pdpMeta',
											'pdp.lst',
											'jointPolicies',
											'writePerWorkerFiles',
											'perVarOutputTypes',
											'doNotReturnRunDataSavePerWorkerOnly'))

# work chunkification ####
cat('chunkification...')
# chunk list and save to workers
workUnitBoundaries <- seq(1,numSample,chunkSizePerWorker*numWorkers)
# in case the chunkSize is not a perfect divisor of the numSample, add numSample as the 
# final boundary
if(workUnitBoundaries[length(workUnitBoundaries)]!=numSample){
	workUnitBoundaries <- c(workUnitBoundaries,numSample)
}
# add one to the last work unit boundary, as during running we always deduct one from the next boundary
workUnitBoundaries[length(workUnitBoundaries)] <- numSample+1
cat('done\n')

# run ####
cat(sprintf('Run of %i runs split up into %i work units of size %i (%i per worker).\n',
						numSample,length(workUnitBoundaries)-1,chunkSizePerWorker*numWorkers,chunkSizePerWorker))
chunkTimes <- c()
completeRunsSoFar <- 0
i <- 0
while(i<(length(workUnitBoundaries)-1)){
	i <- i+1
	workUnit.i <- i
	clusterExport(cl,list('workUnit.i'),envir=environment())
	cat(sprintf('\r(r) Running unit %i: samples %i to %i. ',
							i, workUnitBoundaries[i],workUnitBoundaries[i+1]-1))
	if((workUnitBoundaries[i]-1)>0){
		cat(sprintf('So far: Complete runs %i (%.1f%%)',
								completeRunsSoFar,100*completeRunsSoFar/(workUnitBoundaries[i]-1)))
	}
	if(length(chunkTimes>1)){
		cat(sprintf(', time per unit %i s (%.2f r/s, %.2f r/s/thread), expect completion in %s',
								round(mean(chunkTimes,na.rm=T)),
								length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								if(exists('dseconds')){dseconds(round(mean(chunkTimes,na.rm=T))*
																									(length(workUnitBoundaries)-i))} 
								else{paste0(round(mean(chunkTimes,na.rm=T))*(length(workUnitBoundaries)-i),'s')}))
	}
	tic()
	workUnit <- workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)
	workerWorkUnits <- chunk(workUnit,numWorkers)
	# write the samplePoints of the work units to the worker directories
	for(w.i in workers){
		if(w.i <= length(workerWorkUnits) && !is.null(workerWorkUnits[[w.i]])){
			saveRDS(samplePoints[workerWorkUnits[[w.i]],],
							file.path(name.workDir,paste0(name.workerDirBasename,w.i),'samplePoints.RDS'))
		} else {
			saveRDS(samplePoints[c(),],
							file.path(name.workDir,paste0(name.workerDirBasename,w.i),'samplePoints.RDS'))
		}
	}
	gobble <- clusterEvalQ(cl,{
		samplePoints <- readRDS('samplePoints.RDS')
	})
	parOutput <- clusterEvalQ(cl,runFridaParmsBySamplePoints(policyMode=T))
	timing <- toc(quiet=T)
	cat('\r(d)') # diagnostics
	chunkTimes[i] <- timing$toc-timing$tic
	logLike.df <- data.frame(id=integer(),logLike=double())
	for(r.i in 1:length(parOutput)){
		logLike.df <- rbind(logLike.df,parOutput[[r.i]]$logLike.df)
	}
	logLike <- c()
	logLike[logLike.df$id] <- logLike.df$logLike
	completeRunsSoFar <- sum(logLike > -.Machine$double.xmax+(200*.Machine$double.eps))
	cat('\r   ')
}
cat(sprintf('\r    complete runs %i (%.2f%%), average chunk time %i sec (%.2f r/s, %.2f r/s/thread), over all run time %s %s\n',
						completeRunsSoFar,100*completeRunsSoFar/numSample,
						round(mean(chunkTimes,na.rm=T)),
						length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
						chunkSizePerWorker/mean(chunkTimes,na.rm=T),
						dseconds(round(sum(chunkTimes,na.rm=T))),
						'                                                                             '))
source('cleanup.R')

# merge files ####
numWorkers <- floor(parallel::detectCores()/3)
skipExtraVars <- T
source('clusterHelp.R')
mergePerVarFiles(verbosity = 1,compressCsv=compressCsv)
source('cleanup.R')




