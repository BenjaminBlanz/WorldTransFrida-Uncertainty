source('initialise.R')
source('configPolicyAnalysis.R')
numWorkersArg <- as.numeric(commandArgs(T))
origNumWorkers <- numWorkers
if(length(numWorkersArg)!=1 && !is.numeric(numWorkersArg)){
	numWorkers <- numWorkersFileMerge
} else {
	numWorkers <- numWorkersArg
}
# merge files ####
skipExtraVars <- T
source('clusterHelp.R')
mergePerVarFiles(verbosity = 1,parStrat=1,compressCsv=compressCsv)
source('cleanup.R')
numWorkers <- origNumWorkers

