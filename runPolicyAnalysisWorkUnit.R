# unitID will be replaced by the main runscript
source('initialise.R')
source('configPolicyAnalysis.R')
workUnit.i <- unitID
name.workerDirBasename <- paste0(name.workerDirBasename,unitID,'_')
setwd(file.path(location.baseOutput,expID,'workunit-unitID'))

# start the cluster
source('clusterHelp-expID-unitID.R')
pdpMeta <- readRDS(file.path(location.baseOutput,expID,'pdpMeta.RDS'))
pdp.lst <- readRDS(file.path(location.baseOutput,expID,'pdp.lst.RDS'))
jointPolicies <- saveRDS(file.path(location.baseOutput,expID,'jointPolicies.RDS'))
clusterExport(cl,list('pdpMeta',
											'pdp.lst',
											'jointPolicies',
											'writePerWorkerFiles',
											'perVarOutputTypes',
											'doNotReturnRunDataSavePerWorkerOnly'))

# load the sample points to evaluate/the workunit
samplePoints <- readRDS('samplePoints-workunit-unitID.RDS')
workUnit <- 1:nrow(samplePoints)
workerWorkUnits <- chunk(workUnit,numWorkers)
# list of workers is provided by clusterHelp
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

