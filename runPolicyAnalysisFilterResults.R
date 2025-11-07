# Scriptable version of the filter function.
source('initialise.R')
source('configPolicyAnalysis.R')
source('funPolicyAnalysisFilter.R')
# config cotains the baseline filterSpec
# config also contains desiredFilterSpec

option_list = list(
	make_option(c("-f", "--varFile"), type="character", default=NULL, 
							help="file containing model output for variable to be filtered by", metavar="character"),
	make_option(c("-c", "--useCluster"), type="boolean", default=FALSE, 
							help="should a cluster be started for filtering the years [default= %default]", metavar="character"),
	make_option(c("-d", "--useDesiredFilterSpec"), type="boolean", default=FALSE, 
							help="TRUE to use desiredFilterSpec from config",
							metavar="boolean"),
	make_option(c("-p", "--droppedPolIDs"), type="character", default=NULL, 
							help="name of file containing dropped PolIDs for prefiltering", metavar="character"),
	make_option(c("-v", "--verbosity"), type="character", default=9, 
							help="verbosity 0 is silent", metavar="character"),
	make_option(c("-o", "--location.output"), type="character", default=NULL, 
							help="override the location.output parameter", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if(is.null(opt$varFile)){
	stop('varfile has to be provided\n')
} else {
	varFile <- opt$varFile
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
if(opt$useDesiredFilterSpec){
	useDesiredFilterSpec <- T
	filterSpec <- desiredFilterSpec
} else {
	useDesiredFilterSpec <- F
}
verbosity <- opt$verbosity

outputFolder <- file.path(location.output,'detectedParmSpace','PerVarFiles-RDS')

filterPolicyAnalysisResults(varFile,useCluster,useDesiredFilterSpec,
														droppedPolIDs,verbosity,location.output,
														returnPolIDsToDrop=FALSE)
