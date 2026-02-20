cat('initialising\n')
source('initialise.R')

# workUnit.i has to be supplied in the environment at runtime
# it can be supplied in a scripted context by invoking an -e argument
# e.g.
# Rscript runPolicyAnalysisWorkUnit.R 1
args <- commandArgs(T)
workUnit.i <- as.numeric(args[1])
configFile <- as.character(args[2])
cat('Supplied workUnit.i: ')
cat(workUnit.i)
cat('\n')
if(is.na(workUnit.i)){
	stop('Supplied workunit.i is not a number\n')
}
cat('Supplied configFile: ')
cat(configFile)
cat('\n')
if(sum(grep('\\.R$',configFile))!=1){
	stop('Incorrect arg supplied as config file\n')
}
source(configFile)
location.workunit <- file.path(location.output,'workUnits',paste0('workUnit-',workUnit.i))
cat('workUnit directory:')
cat(location.workunit)
cat('\n')
if(!file.exists(location.workunit)){
	stop('Incorrect workUnit.i supplied as arg\n')
}
cat('location.output:\n')
cat(paste0(location.output,'\n'))
cat(sprintf('WorkUnit.i: %i\n',workUnit.i))
write('running',file.path(location.workunit,'status.txt'),append=F)
name.workerDirBasename <- paste0(origName.workerDirBasename,workUnit.i,'_')

# load calDat and resSigma for likelihood
if(file.exists(file.path(location.output,'calDat.RDS'))){
	calDat.lst <- readRDS(file.path(location.output,'calDat.RDS'))
	calDat <- calDat.lst$calDat
	calDat.impExtrValue <- calDat.lst$calDat.impExtrValue
	calDat.orig <- calDat.lst$calDat.orig
	calDat.withAllVars <- calDat.lst$calDat.withAllVars
} else {
	stop('Missing calDat file. Run runInitialiseData.R first.\n')
}
if(treatVarsAsIndep&&
	 file.exists(file.path(location.output,'sigma-indepParms.RDS'))){
	resSigma <- readRDS(file.path(location.output,'sigma-indepParms.RDS'))
} else if(file.exists(file.path(location.output,'sigma.RDS'))){
	resSigma <- readRDS(file.path(location.output,'sigma.RDS'))
	if(treatVarsAsIndep){
		# get the diagonal elements
		resSigma.var <- diag(resSigma)
		# make a diagonal matrix with those elements
		resSigma <- diag(resSigma.var)
	}
} else {
	stop('Missing covariance matrix file. Run runInitialiseData.R first.\n')
}

# setup frida working directory
source('setupTMPFS.R')
# write export spec
extraVarsToExport <- unique(read.csv(file.path(location.frida.info,name.frida_extra_variables_to_export_list))$FRIDA.FQN)
extraVarsToExport <- extraVarsToExport[nchar(extraVarsToExport)>4]
writeFRIDAExportSpec(extraVarsToExport,location.frida)
# start the cluster
cat('starting cluster\n')
source('clusterHelp.R')
clusterExport(cl,list('writePerWorkerFiles',
											'perVarOutputTypes',
											'doNotReturnRunDataSavePerWorkerOnly',
											'workUnit.i',
											'treatVarsAsIndep'))

# load the sample points to evaluate/the workUnit
cat('reading sample points\n')
samplePoints <- readRDS(file.path(location.workunit,'samplePoints.RDS'))
workUnit <- 1:nrow(samplePoints)
workerWorkUnits <- chunk(workUnit,numWorkers)
# list of workers is provided by clusterHelp
cat('writing per worker sample points\n')
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
cat('cluster runFridaParmsBySamplePoints\n')
parOutput <- clusterEvalQ(cl,runFridaParmsBySamplePoints(policyMode=F))
cat('cleanup\n')
source('cleanup.R')
write('completed',file.path(location.workunit,'status.txt'),append=F)
cat('done\n')
