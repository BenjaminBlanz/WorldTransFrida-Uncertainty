# plot policy analysis output
source('initialise.R')
source('configPolicyAnalysis.R')
source('funPolicyAnalysisPlots.R')
#override location.output
location.output <- 'policy-workOutput/AllPolicies1e6-moreExports'

# save reference file
system(paste('cp','configPolicyAnalysis.R',file.path(location.output,'configPolicyAnalysisUsedForPlotting.R')))

# location.output <- 'policy-workOutput/'

figuresFolder <- file.path(location.output,'figures')
dir.create(figuresFolder,showWarnings = F,recursive = T)

outputFolder <- file.path(location.output,'detectedParmSpace','PerVarFiles-RDS')

# generate missing computable variables ####

generatedVarsMeta <- data.frame("IMPORTANT..NO.WHITESPACE.AT.THE.END.OF.FIELDS."=rep('',3),
																"FRIDA.Module"=rep('',3),
																"FRIDA.Name"=c('real GDP growth rate','surface temperature anomaly growth rate','real GDP per capita'),
																"Index"=rep('[*]',3),
																"FRIDA.FQN"=c('real GDP growth rate','surface temperature anomaly growth rate','real GDP per capita'),
																"Unit" =c('rate','rate','2021c$/p'))
generatedVarsMeta$cleanName <- cleanNames(generatedVarsMeta$FRIDA.FQN)

if(!file.exists(file.path(outputFolder,paste0(generatedVarsMeta$cleanName[1],'.RDS')))){
	tic()
	cat('gdpgr does not exist generating...')
	if(file.exists(file.path(outputFolder,'gdp_real_gdp_in_2021c.RDS'))){
		cat('reading gdp...')
		gdp <- readPerVarFile(file = file.path(outputFolder,'gdp_real_gdp_in_2021c.RDS'))
	} else {
		stop('missing gdp in output files\n')
	}
	cat('calculating...')
	gdpgr <- gdp
	gdpgr[,3:(ncol(gdpgr)-1)] <- (gdpgr[,4:(ncol(gdpgr))] - gdpgr[,3:(ncol(gdpgr)-1)]) /  gdpgr[,3:(ncol(gdpgr)-1)]
	gdpgr[,ncol(gdpgr)] <- NA
	cat('saving...')
	saveRDS(gdpgr,file.path(outputFolder,paste0(generatedVarsMeta$cleanName[1],'.RDS')))
	rm(gdp)
	rm(gdpgr)
	quietgc()
	cat('done\n')
	toc()
}
if(!file.exists(file.path(outputFolder,paste0(generatedVarsMeta$cleanName[2],'.RDS')))){
	tic()
	cat('stagr does not exist generating...')
	if(file.exists(file.path(outputFolder,'energy_balance_model_surface_temperature_anomaly.RDS'))){
		cat('reading sta...')
		sta <- readPerVarFile(file = file.path(outputFolder,'energy_balance_model_surface_temperature_anomaly.RDS'))
	} else {
		stop('missing sta in output files\n')
	}
	cat('calculating...')
	stagr <- sta
	stagr[,3:(ncol(stagr)-1)] <- (stagr[,4:(ncol(stagr))] - stagr[,3:(ncol(stagr)-1)]) /  stagr[,3:(ncol(stagr)-1)]
	stagr[,ncol(stagr)] <- NA
	cat('saving...')
	saveRDS(stagr,file.path(outputFolder,paste0(generatedVarsMeta$cleanName[2],'.RDS')))
	rm(sta)
	rm(stagr)
	quietgc()
	cat('done\n')
	toc()
}
if(!file.exists(file.path(outputFolder,paste0(generatedVarsMeta$cleanName[3],'.RDS')))){
	tic()
	cat('gdppc does not exist generating...')
	if(file.exists(file.path(outputFolder,'demographics_population.RDS')) &&
								 file.exists(file.path(outputFolder,'gdp_real_gdp_in_2021c.RDS'))){
		cat('reading gdp...')
		gdp <- readPerVarFile(file = file.path(outputFolder,'gdp_real_gdp_in_2021c.RDS'))
		cat('reading pop...')
		pop <- readPerVarFile(file = file.path(outputFolder,'demographics_population.RDS'))
	} else {
		stop('missing gdp_real_gdp_in_2021c.RDS or demographics_population.RDS in output files\n')
	}
	cat('calculating...')
	gdppc <- gdp/(pop/1000) # divide by 1000 so that the unit is $/p instead of 1000$/p
	cat('saving...')
	saveRDS(gdppc,file.path(outputFolder,paste0(generatedVarsMeta$cleanName[3],'.RDS')))
	rm(gdp)
	rm(pop)
	rm(gdppc)
	quietgc()
	cat('done\n')
	toc()
}


# read files list ####
varsFiles <- list.files(file.path(outputFolder),pattern='*.RDS')
varsFiles <- varsFiles[varsFiles!='logLike.RDS']
varsMeta <- read.csv(file.path(location.frida.info,name.frida_extra_variables_to_export_list))
varsMeta$cleanName <- cleanNames(varsMeta$FRIDA.FQN)
varsMeta <- rbind(varsMeta,generatedVarsMeta)

# polIDsToDrop.lst <- parLapplyLB(cl,1:length(varsFiles),parFilterResults,
# 														varsFiles = varsFiles,filterSpec=filterSpec,
# 														chunk.size = 1)

# sequential filtering could be faster if we applied the last filter as prefilter
# for the nest, however as long as the data are in data.frames and not data.tables
# the subsetting of the data.frame is so slow it negates any benefit from having to filter
# fewer points. (will have to change the output generation in the main run script to use
# data.tables (no converting is not the solution, takes too long))
polIDsToDrop.lst <- list()
cat('dropping incomplete and inf...')
tic()
cat('reading file to determine incompletes...')
varDat <- readPerVarFile(file.path(outputFolder,varsFiles[1]))
# if the last column is entirely NA this is probably
# a generated variable, drop that col to not mess with the filter
if(sum(is.na(varDat[,ncol(varDat)]))==nrow(varDat)){
	varDat <- varDat[,-ncol(varDat)]
}
numPolIDs <- length(unique(varDat$polID))
cat('filtering...')
polIDsToDrop.lst[[1]] <- unique(varDat$polID[!complete.cases(varDat) | !is.finite(varDat[,ncol(varDat)])])
timing <- toc(quiet=T)
cat('done\n')
rm(varDat)
quietgc()

cat('Filter 0 filtering for incompletes and infinities filter applies in all years\n')
polIDsToDrop <- polIDsToDrop.lst[[1]] 
cat(sprintf('  PolIDs dropped so far: %i (%0.2f%%) (%i in this file %i new) %s\n',
						length(polIDsToDrop), 100*length(polIDsToDrop)/numPolIDs,
						length(polIDsToDrop.lst[[1]]),
						length(polIDsToDrop)-0,
						timing$callback_msg))
for(i in 1:length(filterSpec)){
	filteredFile <- paste0(names(filterSpec)[i],'.RDS')
	cat(sprintf('Filter %i of %i ',i,length(filterSpec)))
	cat(sprintf(' filtering %s %s %s',names(filterSpec)[i],
							filterSpec[[i]]$type,filterSpec[[i]]$level))
	if(!is.null(filterSpec[[i]]$sowID)){
		cat(sprintf(' SOW %i',filterSpec[[i]]$sowID))
	}
	if(!is.null(filterSpec[[i]]$allowedTransgressions)){
		cat(sprintf(' number of allowed transgressions per year %i',filterSpec[[i]]$allowedTransgressions))
	}
	if(!is.null(filterSpec[[i]]$years)){
		rangeFilterYears <- range(filterSpec[[i]]$years)
		numFilterYears <- length(filterSpec[[i]]$years)
		numFilterYearsStr <- ifelse(numFilterYears==(diff(rangeFilterYears)+1),'all',as.character(numFilterYears))
		cat(sprintf(' filter applies in %s years between %i and %i',numFilterYearsStr,rangeFilterYears[1],rangeFilterYears[2]))
	} else {
		cat(' filter applies in all years')
	}
	cat('\n')
	tic()
	if(filteredFile %in% varsFiles){
		# for debugging filterscript
		# varFile <- filteredFile
		# useClusterNumThreads <- min(171/2,detectCores())
		# verbosity <- 9
		# useDesiredFilterSpec <- F
		system(paste('Rscript --max-connections=1024 --no-site-file runPolicyAnalysisFilterResults.R -f',
								 filteredFile, '-c',min(171/2,detectCores()), '-o',location.output))
		polIDsToDrop.lst[[i+1]] <- readRDS(file.path(location.output,'filterResults',
																					paste0(names(filterSpec)[i],'-filter.RDS')))
		polIDsToDrop.old <- polIDsToDrop
		polIDsToDrop <- sort(unique(c(polIDsToDrop,unlist(polIDsToDrop.lst[[i+1]]))))
	} else {
		cat('no such file\n')
	}
	timing <- toc(quiet=T)
	cat(sprintf('  PolIDs dropped so far: %i (%0.2f%%) (%i in this file %i new) %s\n',
							length(polIDsToDrop), 100*length(polIDsToDrop)/numPolIDs,
							length(polIDsToDrop.lst[[i+1]]),
							length(polIDsToDrop)-length(polIDsToDrop.old),
							timing$callback_msg))
	if(length(polIDsToDrop) >= numPolIDs){
		stop('All policies have been filtered out, nothing left to do\n')
	}
}
polIDsToDrop <- sort(unique(unlist(polIDsToDrop.lst)))

# plotting baseline ####
firstThingsToPlot <- c(69,112,131,130,141,106,107)
thingsToPlot <- c(firstThingsToPlot,seq(1:length(varsFiles))[-firstThingsToPlot])
clPlotting <- makeForkCluster(numPlotThreads)
for(plotType in plotTypes){
	cat(sprintf('plotting %i vars with %i threads plot type %s...',
							length(thingsToPlot),numPlotThreads,plotType))
	tic()
	logmax <- log(numInitialJointPol*ifelse(plotType==3,11,1))
	colLevels <- exp(seq(0,logmax,length.out=plot.numColLevels))
	baselinePlotPropsParRes <- parLapplyLB(clPlotting,thingsToPlot,parPlotPolResults,
																	 varsFiles=varsFiles,
																	 polIDsToDrop=polIDsToDrop,
																	 funFigFolder=NULL,
																	 plotType=plotType,
																	 colLevels=colLevels)
	timing <- toc(quiet=T)
	cat(sprintf('done %s\n',timing$callback_msg))
}
stopCluster(clPlotting)
baselinePlotProps <- list()
for(p.i in 1:length(baselinePlotPropsParRes)){
	baselinePlotProps[[thingsToPlot[p.i]]] <- baselinePlotPropsParRes[[p.i]]
}

# desired filtering  and plotting ####
polIDsToDrop <- sort(unique(unlist(polIDsToDrop.lst)))
polIDsToDropDesired.lst <- list()
for(i in 1:length(desiredFilterSpec)){
	filteredFile <- paste0(names(desiredFilterSpec)[i],'.RDS')
	cat(sprintf('Desired filter %i of %i ',i,length(desiredFilterSpec)))
	cat(sprintf(' filtering %s %s %s',names(desiredFilterSpec)[i],
							desiredFilterSpec[[i]]$type,desiredFilterSpec[[i]]$level))
	if(!is.null(desiredFilterSpec[[i]]$sowID)){
		cat(sprintf(' SOW %i',desiredFilterSpec[[i]]$sowID))
	}
	if(!is.null(desiredFilterSpec[[i]]$allowedTransgressions)){
		cat(sprintf(' number of allowed transgressions per year %i',desiredFilterSpec[[i]]$allowedTransgressions))
	}
	if(!is.null(desiredFilterSpec[[i]]$years)){
		rangeFilterYears <- range(desiredFilterSpec[[i]]$years)
		numFilterYears <- length(desiredFilterSpec[[i]]$years)
		numFilterYearsStr <- ifelse(numFilterYears==(diff(rangeFilterYears)+1),'all',as.character(numFilterYears))
		cat(sprintf(' filter applies in %s years between %i and %i',numFilterYearsStr,rangeFilterYears[1],rangeFilterYears[2]))
	} else {
		cat(' filter applies in all years')
	}
	cat('\n')
	tic()
	if(filteredFile %in% varsFiles){
		system(paste('Rscript --max-connections=1024 --no-site-file runPolicyAnalysisFilterResults.R', 
								 '-f', filteredFile, '-c',min(171/2,detectCores()), '-o',location.output,
								 '-d',TRUE))
		polIDsToDropDesired.lst[[i+1]] <- readRDS(file.path(location.output,'filterResults',
																								 paste0(names(desiredFilterSpec)[i],'-desiredFilter.RDS')))
		polIDsToDrop.old <- polIDsToDrop
		polIDsToDrop <- sort(unique(c(polIDsToDrop,unlist(polIDsToDropDesired.lst[[i+1]]))))
	} else {
		cat('no such file\n')
	}
	timing <- toc(quiet=T)
	cat(sprintf('  PolIDs dropped so far: %i (%0.2f%%) (%i in this file %i new) %s\n',
							length(polIDsToDrop), 100*length(polIDsToDrop)/numPolIDs,
							length(polIDsToDropDesired.lst[[i+1]]),
							length(polIDsToDrop)-length(polIDsToDrop.old),
							timing$callback_msg))
	if(length(polIDsToDrop) >= numPolIDs){
		stop('All policies have been filtered out, nothing left to do\n')
	}
	clPlotting <- makeForkCluster(numPlotThreads)
	for(plotType in plotTypes){
		funFigFolder <- file.path(figuresFolder,paste0('plotType',plotType,'-desiredFilters-',i))
		cat(sprintf('plotting plot type %s\n',plotType))
		logmax <- log(numInitialJointPol*ifelse(plotType==3,11,1))
		colLevels <- exp(seq(0,logmax,length.out=plot.numColLevels))
		parRes <- parLapplyLB(clPlotting,thingsToPlot,parPlotPolResults,
													varsFiles=varsFiles,
													polIDsToDrop=polIDsToDrop,
													funFigFolder=funFigFolder,
													plotType=plotType,
													colLevels=colLevels,
													baselinePlotProps=baselinePlotProps)
	}
	stopCluster(clPlotting)
}

# selected run ####
selFilePath <- file.path(outputFolder,paste0(selectedRunSpec$var,'.RDS'))
if(!file.exists(selFilePath)){
	stop(sprintf('File selected for individual run %s does not exist.\n',selFilePath))
} 
varDat <- readPerVarFile(selFilePath)
setForSelection <- varDat[!varDat$polID %in% polIDsToDrop & 
														varDat$sowID==selectedRunSpec$sow,
													c('polID',as.character(selectedRunSpec$year))])
if(selectedRunSpec$optimize=='max'){
	selPolID <- setForSelection$polID[which.max(setForSelection[[as.character(selectedRunSpec$year)]])]
}
if(selectedRunSpec$optimize=='min'){
	selPolID <- setForSelection$polID[which.min(setForSelection[[as.character(selectedRunSpec$year)]])]
}
clPlotting <- makeForkCluster(numPlotThreads)
for(plotType in plotTypes){
	funFigFolder <- file.path(figuresFolder,paste0('plotType',plotType,'-desiredFilters-',i,'-withSelPolID'))
	cat(sprintf('plotting plot type %s\n',plotType))
	logmax <- log(numInitialJointPol*ifelse(plotType==3,11,1))
	colLevels <- exp(seq(0,logmax,length.out=plot.numColLevels))
	parRes <- parLapplyLB(clPlotting,thingsToPlot,parPlotPolResults,
												varsFiles=varsFiles,
												polIDsToDrop=polIDsToDrop,
												funFigFolder=funFigFolder,
												plotType=plotType,
												colLevels=colLevels,
												baselinePlotProps=baselinePlotProps,
												selPolID=selPolID)
}
stopCluster(clPlotting)

# get information on the selPolID
pdp.lst <- readRDS(file.path(location.output,'pdp.lst.RDS'))
pdpMeta <- readRDS(file.path(location.output,'pdpMeta.RDS'))
jointPolicies <- readRDS(file.path(location.output,'jointPolicies.RDS'))
samplePoints <- readRDS(file.path(location.output,'samplePoints.RDS'))

selPolDescStrs <- c()
i <- 0
for(domID in colnames(samplePoints)){
	if(!is.na(samplePoints[selPolID,domID])){
		i <- i+1
		pdpName <- pdpMeta$domain[pdpMeta$domID==domID][1]
		sdmID <- jointPolicies$sdmID[jointPolicies$domID==domID & jointPolicies$dplID==samplePoints[selPolID,domID]]
		if(!is.na(pdpMeta$subdomain[pdpMeta$domID==domID & pdpMeta$sdmID==sdmID])){
			pdpName <- paste0(pdpName,'+',pdpMeta$subdomain[pdpMeta$domID==domID & pdpMeta$sdmID==sdmID])
		} 
		sdpID <- jointPolicies$sdpID[jointPolicies$domID==domID & jointPolicies$dplID==samplePoints[selPolID,domID]]
		selPolDescStrs <-c(selPolDescStrs, pdp.lst[[pdpName]][pdp.lst[[pdpName]]$polID==sdpID,2])
	}
}
write(selPolDescStrs,file.path(figuresFolder,'selRunPolDesc.csv'))
