# plot policy analysis output
source('initialise.R')
source('configPolicyAnalysis.R')
source('funPolicyAnalysisPlots.R')
#override location.output
location.output <- 'policy-workOutput/AllPolicies1e6-moreExports'
# location.output <- 'policy-workOutput/'

figuresFolder <- file.path(location.output,'figures')
dir.create(figuresFolder,showWarnings = F,recursive = T)

outputFolder <- file.path(location.output,'detectedParmSpace','PerVarFiles-RDS')

if(!file.exists(file.path(outputFolder,'gdpgr.RDS'))){
	cat('gdpgr does not exist calculating...')
	if(file.exists(file.path(outputFolder,'gdp_real_gdp_in_2021c.RDS'))){
		gdp <- readPerVarFile(file = file.path(outputFolder,'gdp_real_gdp_in_2021c.RDS'))
	} else {
		stop('missing gdp in output files\n')
	}
	gdpgr <- gdp
	gdpgr[,3:(ncol(gdpgr)-1)] <- (gdpgr[,4:(ncol(gdpgr))] - gdpgr[,3:(ncol(gdpgr)-1)]) /  gdpgr[,3:(ncol(gdpgr)-1)]
	gdpgr[,ncol(gdpgr)] <- NA
	cat('saving...')
	saveRDS(gdpgr,file.path(outputFolder,'gdpgr.RDS'))
	cat('done\n')
}
if(!file.exists(file.path(outputFolder,'stagr.RDS'))){
	cat('stagr does not exist calculating...')
	if(file.exists(file.path(outputFolder,'energy_balance_model_surface_temperature_anomaly.RDS'))){
		sta <- readPerVarFile(file = file.path(outputFolder,'energy_balance_model_surface_temperature_anomaly.RDS'))
	} else {
		stop('missing sta in output files\n')
	}
	stagr <- sta
	stagr[,3:(ncol(stagr)-1)] <- (stagr[,4:(ncol(stagr))] - stagr[,3:(ncol(stagr)-1)]) /  stagr[,3:(ncol(stagr)-1)]
	stagr[,ncol(stagr)] <- NA
	cat('saving...')
	saveRDS(stagr,file.path(outputFolder,'stagr.RDS'))
	cat('done\n')
}

# read files list ####
varsFiles <- list.files(file.path(outputFolder),pattern='*.RDS')
varsFiles <- varsFiles[varsFiles!='logLike.RDS']
varsMeta <- read.csv(file.path(location.frida.info,name.frida_extra_variables_to_export_list))
varsMeta$cleanName <- cleanNames(varsMeta$FRIDA.FQN)

# polIDsToDrop.lst <- parLapplyLB(cl,1:length(varsFiles),parFilterResults,
# 														varsFiles = varsFiles,filterSpec=filterSpec,
# 														chunk.size = 1)

# sequential filtering could be faster
polIDsToDrop.lst <- list()
for(i in 1:length(filterSpec)){
	filteredFile <- paste0(names(filterSpec)[i],'.RDS')
	cat(sprintf('reading for filtering %s ',names(filterSpec)[i]))
	cat(sprintf('filtered file %i of %i\n ',i,length(filterSpec)))
	polIDsToDrop <- c()
	if(filteredFile %in% varsFiles){
		system(paste('Rscript --max-connections=1024 --no-site-file runPolicyAnalysisFilterResults.R -f',
								 filteredFile, '-c','TRUE', '-o',location.output))
		polIDsToDrop.lst[[i]] <- readRDS(file.path(location.output,'filterResults',
																					paste0(names(filterSpec)[i],'-filter.RDS')))
		polIDsToDrop.old <- polIDsToDrop
		polIDsToDrop <- sort(unique(unlist(polIDsToDrop.lst)))
	} else {
		cat('no such file\n')
	}
	cat(sprintf('PolIDs dropped so far: %i (%i new from this file)\n',
							length(polIDsToDrop),sum(!polIDsToDrop.lst[[i]] %in% polIDsToDrop.old)))
}
polIDsToDrop <- sort(unique(unlist(polIDsToDrop.lst)))
saveRDS(polIDsToDrop,file.path(location.output,'droppedPolIDs.RDS'))
polIDsToDrop <- readRDS(file.path(location.output,'droppedPolIDs.RDS'))
firstThingsToPlot <- c(69,112,106,107)
thingsToPlot <- c(firstThingsToPlot)#,seq(1:length(varsFiles))[-firstThingsToPlot])
clPlotting <- makeForkCluster(numPlotThreads)
for(plotType in plotTypes){
	cat(sprintf('plotting plot type %s\n',plotType))
	logmax <- log(numInitialJointPol*ifelse(plotType==3,11,1))
	colLevels <- exp(seq(0,logmax,length.out=plot.numColLevels))
	parRes <- parLapplyLB(clPlotting,thingsToPlot,parPlotPolResults,
												varsFiles=varsFiles,
												polIDsToDrop=polIDsToDrop,
												funFigFolder=NULL,
												plotType=plotType,
												colLevels=colLevels)
}
stopCluster(clPlotting)

# desired filtering ####
desiredFilterSpec <- list()
desiredFilterSpec$energy_balance_model_surface_temperature_anomaly <- c('sgtval',2,5)
saveRDS(desiredFilterSpec,file.path(location.output,'desiredFilterSpec.RDS'))
polIDsToDropDesired.lst <- list()
for(i in 1:length(desiredFilterSpec)){
	filteredFile <- paste0(names(desiredFilterSpec)[i],'.RDS')
	cat(sprintf('reading for desired filtering %s ',names(desiredFilterSpec)[i]))
	cat(sprintf('filtered file %i of %i\n ',i,length(desiredFilterSpec)))
	polIDsToDropDesired <- c()
	if(filteredFile %in% varsFiles){
		system(paste('Rscript --max-connections=1024 --no-site-file runPolicyAnalysisFilterResults.R -f',
								 filteredFile, '-c','TRUE', '-o',location.output,'-d','desiredFilterSpec.RDS'))
		polIDsToDropDesired.lst[[i]] <- readRDS(file.path(location.output,'filterResults',
																							 paste0(names(desiredFilterSpec)[i],'-filter.RDS')))
		polIDsToDrop.old <- polIDsToDrop
		polIDsToDrop <- sort(unique(unlist(polIDsToDrop.lst)))
	} else {
		cat('no such file\n')
	}
	cat(sprintf('PolIDs dropped so far: %i (%i new from this file)\n',
							length(polIDsToDrop),sum(!polIDsToDrop.lst[[i]] %in% polIDsToDrop.old)))
}
polIDsToDropDesired <- unique(polIDsToDropDesired,polIDsToDrop)
funFigFolder <- file.path(location.output,'figures',paste0('plotType',plotType,'-desiredFilters'))
clPlotting <- makeForkCluster(numPlotThreads)
for(plotType in plotTypes){
	cat(sprintf('plotting plot type %s\n',plotType))
	logmax <- log(numInitialJointPol*ifelse(plotType==3,11,1))
	colLevels <- exp(seq(0,logmax,length.out=plot.numColLevels))
	parRes <- parLapplyLB(clPlotting,thingsToPlot,parPlotPolResults,
												varsFiles=varsFiles,
												polIDsToDrop=polIDsToDropDesired,
												funFigFolder=funFigFolder,
												plotType=plotType,
												colLevels=colLevels)
}
stopCluster(clPlotting)
