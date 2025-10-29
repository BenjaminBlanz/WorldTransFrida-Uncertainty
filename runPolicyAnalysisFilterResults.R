source('initialise.R')
source('configPolicyAnalysis.R')

option_list = list(
	make_option(c("-f", "--varfile"), type="character", default=NULL, 
							help="file containing model output for variable to be filtered by", metavar="character"),
	make_option(c("-c", "--useCluster"), type="character", default="True", 
							help="should a cluster be started for filtering the years [default= %default]", metavar="character"),
	make_option(c("-d", "--desiredFilterFileName"), type="character", default=NULL, 
							help="name of file containing filters for desired filtering", metavar="character"),
	make_option(c("-p", "--droppedPolIDs"), type="character", default=NULL, 
							help="name of file containing dropped PolIDs for prefiltering", metavar="character"),
	make_option(c("-v", "--verbosity"), type="character", default=9, 
							help="verbosity 0 is silent", metavar="character"),
	make_option(c("-o", "--location.output"), type="character", default=NULL, 
							help="override the location.output parameter", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(is.null(opt$varfile)){
	stop('varfile has to be provided\n')
} else {
	varFile <- opt$varfile
}
useCluster <- opt$useCluster
if(!is.null(opt$location.output)){
	location.output <- opt$location.output
}
if(!is.null(opt$droppedPolIDs)){
	polIDsToDrop <- readRDS(file.path(location.output,opt$droppedPolIDs))
} else {
	polIDsToDrop <- NULL
}
if(!is.null(opt$desiredFilterFileName)){
	filterSpec <- readRDS(file.path(location.output,opt$desiredFilterFileName))
}
verbosity <- opt$verbosity

outputFolder <- file.path(location.output,'detectedParmSpace','PerVarFiles-RDS')


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
	if(!is.null(polIDsToDrop)){
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
		} else if (filterSpec[[varName]][1] == 'sltval'){
			polIDsToDrop <- varDat$polID[varDat$sowID==filterSpec[[varName]][3] & 
																	 	varDat[[year]]<filterSpec[[varName]][2]]
		} else if (filterSpec[[varName]][1] == 'sgtval'){
			polIDsToDrop <- varDat$polID[varDat$sowID==filterSpec[[varName]][3] &
																	 	varDat[[year]]>filterSpec[[varName]][2]]
		} else {
			stop('unkown filter spec\n')
		}
		if(verbosity>0){cat('.')}
		return(polIDsToDrop[!is.na(polIDsToDrop)])
	}
	if(verbosity>0){cat('processing years')}
	if(useCluster != T){
		if(verbosity>0){cat(' not using cluster...')}
		yearPolIDsToDrop <- lapply(1:length(years),filterFun)
	} else {
		if(verbosity>0){cat(' starting cluster...')}
		clFiltering <- makeForkCluster(min(floor(171/2),detectCores()))
		if(verbosity>0){cat('running...')}
		yearPolIDsToDrop <- parLapply(clFiltering,1:length(years),filterFun)
		stopCluster(clFiltering)
	}
	polIDsToDrop <- unique(c(polIDsToDrop,unlist(yearPolIDsToDrop)))
	writeToFolder <- file.path(location.output,'filterResults')
	dir.create(writeToFolder,showWarnings = F,recursive = T)
	saveRDS(polIDsToDrop,file.path(writeToFolder,paste0(varName,'-filter.RDS')))
} else {
	cat(sprintf('Var name %s not in filter spec.\n',varName))
}
if(verbosity>0){cat('done\n')}
