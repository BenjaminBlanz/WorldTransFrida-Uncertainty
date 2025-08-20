cat('initialising\n')
source('initialise.R')
source('configPolicyAnalysis.R')
# workUnit.i has to be supplied in the environment at runtime
# it can be supplied in a scripted context by invoking an -e argument
# e.g.
# Rscript runPolicyAnalysisWorkUnit.R 1
args <- commandArgs(T)
if(length(args)!=1 && !is.numeric(args)){
	stop('Incorrect arg supplied as workUnit.i\n')
}
if(!file.exists(file.path(location.output,'workUnits',paste0('workUnit-',workUnit.i)))){
	stop('Incorrect workUnit.i supplied as arg\n')
}
workUnit.i <- commandArgs(T)
cat('location.output:\n')
cat(paste0(location.output,'\n'))
cat(sprintf('WorkUnit.i: %i\n',workUnit.i))
write('running',file.path(location.output,'workUnits',paste0('workUnit-',workUnit.i),'status.txt'),append=F)
name.workerDirBasename <- paste0(name.workerDirBasename,workUnit.i,'_')

# start the cluster
cat('starting cluster\n')
source('clusterHelp.R')
pdpMeta <- readRDS(file.path(location.output,'pdpMeta.RDS'))
pdp.lst <- readRDS(file.path(location.output,'pdp.lst.RDS'))
jointPolicies <- readRDS(file.path(location.output,'jointPolicies.RDS'))
clusterExport(cl,list('pdpMeta',
											'pdp.lst',
											'jointPolicies',
											'writePerWorkerFiles',
											'perVarOutputTypes',
											'doNotReturnRunDataSavePerWorkerOnly',
											'workUnit.i'))

# load the sample points to evaluate/the workUnit
cat('reading sample points\n')
samplePoints <- readRDS(file.path(location.output,'workUnits',paste0('workUnit-',workUnit.i),'samplePoints.RDS'))
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
parOutput <- clusterEvalQ(cl,runFridaParmsBySamplePoints(policyMode=T))
cat('cleanup\n')
source('cleanup.R')
write('completed',file.path(location.output,'workUnits',paste0('workUnit-',workUnit.i),'status.txt'),append=F)
cat('done\n')
