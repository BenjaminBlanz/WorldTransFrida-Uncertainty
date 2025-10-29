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
		system(paste('Rscript --max-connections=1024 --no-site-file runPolicyAnalysisFilterResults.R',
								 filteredFile, 'TRUE', location.output))
		polIDsToDrop.lst[[i]] <- readRDS(file.path(location.output,'filterResults',
																					paste0(names(filterSpec)[i],'-filter.RDS')))
		polIDsToDrop <- sort(unique(unlist(polIDsToDrop.lst)))
	} else {
		cat('no such file\n')
	}
	cat(sprintf('PolIDs dropped so far: %i (%i new from this file)\n',
							length(polIDsToDrop),length(polIDsToDrop.lst[[i]])))
}
polIDsToDrop <- sort(unique(unlist(polIDsToDrop.lst)))
saveRDS(polIDsToDrop,file.path(location.output,'droppedPolIDs.RS'))
polIDsToDrop <- readRDS(file.path(location.output,'droppedPolIDs.RS'))
firstThingsToPlot <- c(69,112)
thingsToPlot <- c(firstThingsToPlot,seq(1:length(varsFiles)[-firstThingsToPlot]))
clPlotting <- makeForkCluster(numPlotThreads)
parRes <- parLapplyLB(clPlotting,thingsToPlot,parPlotPolResults,
											varsFiles=varsFiles,
											polIDsToDrop=polIDsToDrop,
											figuresFolder=NULL,
											plotType=2)
stopCluster(clPlotting)


