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
source('cleanup.R')
numWorkers <- origNumWorkers
write('completed',file.path(location.output,'mergeStatus.txt'))

