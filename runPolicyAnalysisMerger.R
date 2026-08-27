source('initialise.R')
source('configPolicyAnalysis.R')
numWorkersArg <- as.numeric(commandArgs(T))
origNumWorkers <- numWorkers
if(length(numWorkersArg)==1 && is.numeric(numWorkersArg)){
	numWorkers <- numWorkersArg
} else {
	numWorkers <- numWorkersFileMerge
}
# merge files ####
skipExtraVars <- T
source('clusterHelp.R')
mergePerVarFiles(verbosity = 1,parStrat=2,compressCsv=compressCsv)
# completion of the ensemble, alongside the merged files
runStatus <- tryCatch(funReadRunStatus(location.output,outputType = perVarOutputTypes[1],
																			 policyMode = T),
											error=function(e){
												warning(sprintf('could not read the merged run status: %s',
																				conditionMessage(e)),call.=FALSE,immediate.=TRUE)
												NULL
											})
if(!is.null(runStatus)){
	funWriteRunStatusPlainCsv(runStatus,location.output)
	funAppendRunCompletionSummary(runStatus=runStatus,location.output=location.output)
}
source('cleanup.R')
numWorkers <- origNumWorkers

