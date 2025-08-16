
# merge files ####
origNumWorkers <- numWorkers
numWorkers <- numWorkersFileMerge
skipExtraVars <- T
source('clusterHelp.R')
mergePerVarFiles(verbosity = 1,compressCsv=compressCsv)
source('cleanup.R')
numWorkers <- origNumWorkers

