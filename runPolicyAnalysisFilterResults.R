source('initialise.R')
source('configPolicyAnalysis.R')
require(parallel)

args = commandArgs(trailingOnly=TRUE)
# args <- c('inflation_inflation_rate.RDS', 'TRUE', 'policy-workOutput/AllPolicies1e6-moreExports/detectedParmSpace/PerVarFiles-RDS')
varFile <- args[1]
useCluster <- args[2]
if(length(args)>=3){
	location.output <- args[3]
}
verbosity <- 9

outputFolder <- file.path(location.output,'detectedParmSpace','PerVarFiles-RDS')

polIDsToDrop <- NULL

varName <- tools::file_path_sans_ext(varFile)
if(varName %in% names(filterSpec)){
	if(verbosity>0){cat(sprintf('reading %s...',varName))}
	varDat <- readPerVarFile(file.path(outputFolder,varFile))
	if(verbosity>0){cat('determining filtered...')}
	# if the last column is entirely NA this is probably
	# a generated variable, drop that col to not mess with the complete cases
	# filter
	if(sum(is.na(varDat[,ncol(varDat)]))==nrow(varDat)){
		varDat <- varDat[,-ncol(varDat)]
	}
	# if we have ex ante known polIDs to drop, we do not have to filter them
	if(is.null(polIDsToDrop)){
		if(verbosity>0){cat('applying prefilter...')}
		varDat <- varDat[!varDat$polID %in% polIDsToDrop,]
	}
	cat('dropping incomplete and inf...')
	polIDsToDrop <- unique(varDat$polID[!complete.cases(varDat) | !is.finite(varDat[,ncol(varDat)])])
	# years
	years <- colnames(varDat)[-c(1,2)]
	filterFun <- function(year.i){
		year <- years[year.i]
		if(filterSpec[[varName]][1] == 'ltabs'){
			polIDsToDrop <- varDat$polID[abs(varDat[[year]])<filterSpec[[varName]][2]]
		} else if(filterSpec[[varName]][1] == 'gtabs'){
			polIDsToDrop <- varDat$polID[abs(varDat[[year]])>filterSpec[[varName]][2]]
		} else if (filterSpec[[varName]][1] == 'ltval'){
			polIDsToDrop <- varDat$polID[varDat[[year]]<filterSpec[[varName]][2]]
		} else if (filterSpec[[varName]][1] == 'gtval'){
			polIDsToDrop <- varDat$polID[varDat[[year]]>filterSpec[[varName]][2]]
		} else {
			stop('unkown filter spec\n')
		}
		if(verbosity>0){cat('.')}
	}
	if(verbosity>0){cat('processing years')}
	if(useCluster != T){
		if(verbosity>0){cat(' not using cluster...')}
		yearPolIDsToDrop <- lapply(1:length(years),filterFun)
	} else {
		if(verbosity>0){cat(' starting cluster...')}
		clFiltering <- makeForkCluster(min(171,detectCores()))
		if(verbosity>0){cat('running...')}
		yearPolIDsToDrop <- parLapply(clFiltering,1:length(years),filterFun)
		stopCluster(clFiltering)
	}
	polIDsToDrop <- unique(c(polIDsToDrop,unlist(yearPolIDsToDrop)))
	writeToFolder <- file.path(location.output,'filterResults')
	dir.create(writeToFolder,showWarnings = F,recursive = T)
	saveRDS(polIDsToDrop,file.path(writeToFolder,paste0(varName,'-filter.RDS')))
}
if(verbosity>0){cat('done\n')}
