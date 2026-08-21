# Merges the chunk files of a single variable in a process of its own.
# Invoked by workerMergePerVarFilesIndepProc (mergePerVarFiles parStrat 3).
args <- commandArgs(T)
if(length(args)<9){
	stop('usage: workerFileMergeScript.R v.i varNamesFile chunkFolder outputFolder outputTypes verbosity compressCsv rdsCompress fullPrecision [baseWD]\n')
}
v.i <- as.numeric(args[1])
varNamesFileName <- as.character(args[2])
chunkFolder <- as.character(args[3])
outputFolder <- as.character(args[4])
outputTypes <- strsplit(as.character(args[5]),',')[[1]]
verbosity <- as.numeric(args[6])
compressCsv <- as.logical(args[7])
rdsCompress <- as.logical(args[8])
fullPrecision <- as.logical(args[9])
baseWD <- if(length(args)>=10){as.character(args[10])} else {getwd()}

suppressPackageStartupMessages(library(data.table,quietly=T,warn.conflicts = F))
# one variable per process, no additional OpenMP threads please
data.table::setDTthreads(1)
source(file.path(baseWD,'funRunFRIDA.R'))
source(file.path(baseWD,'naturalsort.R'))

varNames <- readRDS(varNamesFileName)
workerMergePerVarFiles(v.i,varNames=varNames,
											 chunkFolder=chunkFolder,
											 outputFolder=outputFolder,
											 outputTypes=outputTypes,
											 verbosity=verbosity,
											 compressCsv=compressCsv,
											 rdsCompress=rdsCompress,
											 fullPrecision=fullPrecision)
