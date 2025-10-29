# plot policy analysis output
source('initialise.R')
source('configPolicyAnalysis.R')

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




# plot function ####
# In the first run
# plot unfiltered and apply filter and report back the filtered out IDs
# In the second run 
# plot filtered runs
filterResults <- function(varFile,filterSpec,verbosity=0){
	varName <- tools::file_path_sans_ext(varFile)
	if(verbosity>0){cat(sprintf('reading %s...',varName))}
	if(varName %in% names(filterSpec)){
		varDat <- readPerVarFile(file.path(outputFolder,varFile))
		cat('determining filtered...')
		if(sum(is.na(varDat[,ncol(varDat)]))==nrow(varDat)){
			varDat <- varDat[,-ncol(varDat)]
		}
		polIDsToDrop <- unique(varDat$polID[!complete.cases(varDat) | !is.finite(varDat[,ncol(varDat)])])
		varDat <- varDat[!varDat$polID %in% polIDsToDrop,]
		for(year in colnames(varDat)[-c(1,2)]){
			if(filterSpec[[varName]][1] == 'ltabs'){
				polIDsToDrop <- unique(c(polIDsToDrop,varDat$polID[abs(varDat[[year]])<filterSpec[[varName]][2]]))
			} else if(filterSpec[[varName]][1] == 'gtabs'){
				polIDsToDrop <- unique(c(polIDsToDrop,varDat$polID[abs(varDat[[year]])>filterSpec[[varName]][2]]))
			} else if (filterSpec[[varName]][1] == 'ltval'){
				polIDsToDrop <- unique(c(polIDsToDrop,varDat$polID[varDat[[year]]<filterSpec[[varName]][2]]))
			} else if (filterSpec[[varName]][1] == 'gtval'){
				polIDsToDrop <- unique(c(polIDsToDrop,varDat$polID[varDat[[year]]>filterSpec[[varName]][2]]))
			} else {
				stop('unkown filter spec\n')
			}
		}
	} else {
		polIDsToDrop <- c()
	}
	cat('done\n')
	return(polIDsToDrop)
}
parFilterResults <- function(i,varsFiles,filterSpec){
	return(filterResults(varsFiles[i],filterSpec))
}

plotPolResults <- function(varFile,polIDsToDrop=NULL,figuresFolder=NULL,verbosity=0){
	if(is.null(figuresFolder)){
		figuresFolder <- file.path(location.output,'figures')
	}
	dir.create(figuresFolder,showWarnings = F,recursive = T)
	varName <- tools::file_path_sans_ext(varFile)
	if(verbosity>0){cat(sprintf('reading %s...',varName))}
	varDat <- readPerVarFile(file.path(outputFolder,varFile))
	varUnit <- varsMeta$Unit[varsMeta$cleanName==varName]
	varFullName <- varsMeta$FRIDA.FQN[varsMeta$cleanName==varName]
	cat('applying filters...')
	varDat <- varDat[!varDat$polID %in% polIDsToDrop,]
	if(verbosity>0){cat('plotting...')}
	years <- as.numeric(colnames(varDat)[3:ncol(varDat)])
	maxBoundUp <- c() # area reached by the any SOW
	maxBoundDown <- c() # lower bound of area reached by any SOW
	medianBoundUp <- c() # upper bound of median results
	medianBoundDown <- c() # lower bound of median results
	for(year.i in 1:length(years)){
		year <- years[year.i]
		maxBoundUp[year.i] <- max(varDat[[as.character(year)]],na.rm=T)
		maxBoundDown[year.i] <- min(varDat[[as.character(year)]],na.rm=T)
		medianBoundUp[year.i] <- max(varDat[[as.character(year)]][varDat$sowID==5],na.rm=T)
		medianBoundDown[year.i] <- min(varDat[[as.character(year)]][varDat$sowID==5],na.rm=T)
		if(verbosity>1){cat('.')}
	}
	ylims <- quantile(varDat[,seq(3,ncol(varDat),3)],probs=plot.relyrange,na.rm=T)
	maxBoundUpLim <- maxBoundUp
	maxBoundUpLim[maxBoundUp>ylims[2]] <- ylims[2]
	maxBoundDownLim <- maxBoundDown
	maxBoundDownLim[maxBoundDown<ylims[1]] <- ylims[1]
	medianBoundUpLim <- medianBoundUp
	medianBoundUpLim[medianBoundUp>ylims[2]] <- ylims[2]
	medianBoundDownLim <- medianBoundDown
	medianBoundDownLim[medianBoundDown<ylims[1]] <- ylims[1]
	png(file.path(figuresFolder,paste0(varName,'.png')),
			width = plotWidth,height = plotHeight,units = plotUnit,res = plotRes)
	plot(0,0,type='n',
			 xlim=range(years),
			 ylim=ylims,
			 xlab='year',
			 xaxs='i',
			 main=varFullName,ylab=varUnit)
	grid(col='gray',lwd=2)
	polygon(c(years,rev(years)),c(maxBoundUp,rev(maxBoundDown)),
					border = plot.lcol[1],col = plot.col[1],
					lty=plot.lty[1])
	polygon(c(years,rev(years)),c(medianBoundUp,rev(medianBoundDown)),
					border = plot.lcol[2],col = plot.col[2],
					lty=plot.lty[2])
	dev.off()
	if(verbosity>0){cat('done\n')}
}
parPlotPolResults<-function(i,varsFiles,polIDsToDrop,figuresFolder=NULL){
	plotPolResults(varsFiles[i],polIDsToDrop,figuresFolder)	
}

# read files list ####
varsFiles <- list.files(file.path(outputFolder),pattern='*.RDS')
varsFiles <- varsFiles[varsFiles!='logLike.RDS']
varsMeta <- read.csv(file.path(location.frida.info,name.frida_extra_variables_to_export_list))
varsMeta$cleanName <- cleanNames(varsMeta$FRIDA.FQN)

if(sequentialPlotting){
	numFiles <- length(varsFiles)
	polIDsToDrop.lst <- list()
	cat('PlottingPass 1 unfiltered plots, collecting filter outputs\n')
	for(i in 1:length(filters)){
		filteredFile <- paste0(names(filters)[i],'.RDS')
		cat(sprintf('reading for filtering %s ',names(filters)[i]))
		cat(sprintf('filtered file %i of %i\n ',i,length(filters)))
		if(filteredFile %in% varsFiles){
			polIDsToDrop.lst[[i]] <- filterResults(filteredFile,filterSpec = filters, verbosity = 1)
			cat('\n')
		} else {
			cat('no such file\n')
		}
	}
	polIDsToDrop <- sort(unique(unlist(polIDsToDrop.lst)))
	cat('Plotting pass 2 filtered plots.\n')
	for(i in 1:length(varsFiles)){
		varFile <- varsFiles[i]
		cat(sprintf('File %i of %i ',i,numFiles))
		try{
			plotPolResults(varFile,polIDsToDrop = polIDsToDrop,verbosity = 1)
		}
	}
} else {
	library('parallel')
	cl <- makeForkCluster(numPlotThreads)
	parFilterRes <- parLapplyLB(cl,1:length(varsFiles),parFilterResults,
															varsFiles = varsFiles,filterSpec=filters,
															chunk.size = 1)
	polIDsToDrop <- sort(unique(unlist(parFilterRes)))
	saveRDS(polIDsToDrop,file.path(location.output,'droppedPolIDs.RS'))
	
}


























vars <- list()
vars[['sta']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','energy_balance_model_surface_temperature_anomaly.RDS'))
vars[['pop']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','demographics_population.RDS'))
vars[['gdp']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','gdp_real_gdp_in_2021c.RDS'))
vars[['inf']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','inflation_inflation_rate.RDS'))
vars[['rdu']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','gdp_future_year_in_recession.RDS'))

vars$sta1d <- vars$sta
vars$sta1d[,3:(ncol(vars$sta1d)-1)] <- vars$sta1d[,4:(ncol(vars$sta1d))] - vars$sta1d[,3:(ncol(vars$sta1d)-1)]
vars$sta1d[,ncol(vars$sta1d)] <- NA

vars$gdpgr <- vars$gdp
vars$gdpgr[,3:(ncol(vars$gdpgr)-1)] <- (vars$gdpgr[,4:(ncol(vars$gdpgr))] - vars$gdpgr[,3:(ncol(vars$gdpgr)-1)]) /  vars$gdpgr[,3:(ncol(vars$gdpgr)-1)]
vars$gdpgr[,ncol(vars$gdpgr)] <- NA

#define filter
# the filter defines the things we keep
filter <- complete.cases(vars$sta)
for(var in names(filters)){
	filter <- filter & rowSums(abs(vars[[var]][,-c(1,2)]) > filters[[var]],na.rm = T) == 0 
}
# apply filters
numPolID <- length(unique(vars$sta[,'polID']))
polIDsToRemove <- unique(vars$sta[!filter,'polID'])
filterInclSow <- vars$sta[,'polID'] %in% polIDsToRemove
for(var in names(vars)){
	vars[[var]] <- vars[[var]][filterInclSow,]
}
gc(verbose=F)
cat(sprintf('Filters applied: %i removed , %i remain\n',length(polIDsToRemove),numPolID-length(polIDsToRemove)))

# selecting run
# highest GDP in 2150, that stays below 2deg, in median
# idxOfSubSel <- which.max(vars$gdp[['2150']][vars$sta[['2150']]<2 & vars$gdp$sowID==5])
# selRun <- rownames(vars$gdp[vars$sta[['2150']]<2 & vars$gdp$sowID==5,])[idxOfSubSel]

# highest GDP in 2150, that stays below 2deg and has no more than 10 years in recession, in median
idxOfSubSel <- which.max(vars$gdp[['2150']][vars$sta[['2150']]<2 & vars$rdu[['2150']]<10 & vars$gdp$sowID==5])
selRun <- rownames(vars$gdp[vars$sta[['2150']]<2 & vars$rdu[['2150']]<10 & vars$gdp$sowID==5,])[idxOfSubSel]

# highest real gdp in 2075
# idxOfSubSel <- which.max(vars$gdp[['2075']][vars$gdp$sowID==5])
# selRun <- rownames(vars$gdp[vars$gdp$sowID==5,])[idxOfSubSel]

selPolID <- vars$gdp[selRun,1]

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

# plot ####
plotLims <- list()
plotLims$sta <- c(0,7)
plotLims$gdp <- c(0,2e6)
plotLims$pop <- c(0,1.5e4)
plotLims$inf <- c(-0.4,0.4)
plotLims$sta1d <- c(-0.4,0.4)
plotLims$gdpgr <- c(-0.4,0.4)
plotLims$rdu <- c(0,100)
varUnits <- list()
varUnits$sta <- 'degC'
varUnits$gdp <- 'billion constant 2021 $'
varUnits$pop <- 'million people'
varUnits$sta1d <- 'degC'
varUnits$inf <- 'rate'
varUnits$gdpgr <- 'rate'
varUnits$rdu <- 'years'
varNames <- list()
varNames$sta <- 'Surface Temperature Anomaly'
varNames$gdp <- 'Real GDP'
varNames$pop <- 'Population'
varNames$inf <- 'Inflation'
varNames$sta1d <- 'Change in Surface Temperature Anomaly'
varNames$gdpgr <- 'Real GDP Growth Rate'
varNames$rdu <- 'Duration spent in Recession'
years <- seq(2150,1980,-2)
for(var in names(vars)){
	cat(sprintf('plotting %s ...',var))
	png(file.path(figuresFolder,paste0(var,'.png')),width = plt.w,height = plt.h,units = plt.u,res = plt.r)
	plot(0,0,type='n',xlim=range(years),ylim=plotLims[[var]],
			 xlab='year',
			 main=varNames[[var]],ylab=varUnits[[var]])
	grid(col='gray',lwd=2)
	for(year in years){
		boxplot(vars[[var]][[as.character(year)]][vars[[var]]$sowID==5],at=year,width = 0.8,add = T,pch='.',lwd=1,
						axes=F)
	}
	lines(as.numeric(colnames(vars[[var]][,-c(1,2)])),
				vars[[var]][selRun,-c(1,2)])
	dev.off()
	cat('\n')
}


