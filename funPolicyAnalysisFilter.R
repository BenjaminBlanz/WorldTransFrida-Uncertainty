parFilterPolicyAnalysisResults <- function(i,varsFiles,useCluster=T,useDesiredFilterSpec=F,
																					 droppedPolIDs=NULL,verbosity=1,
																					 overrideLocation.output=NULL,returnPolIDsToDrop=T){
	filterPolicyAnalysisResults(varsFiles[i],useCluster,useDesiredFilterSpec,
															droppedPolIDs,verbosity,
															overrideLocation.output,returnPolIDsToDrop)
}
filterPolicyAnalysisResults <- function(varFile,useCluster=T,useDesiredFilterSpec=F,
																				droppedPolIDs=NULL,verbosity=1,
																				overrideLocation.output=NULL,returnPolIDsToDrop=T){
	if(!is.null(overrideLocation.output)){
		location.output <- overrideLocation.output
	}
	outputFolder <- file.path(location.output,'detectedParmSpace','PerVarFiles-RDS')
	writeToFolder <- file.path(location.output,'filterResults')
	if(useDesiredFilterSpec){
		filterSpec <- desiredFilterSpec
		if(file.exists(file.path(writeToFolder,paste0(varName,'-desiredFilterSpec.RDS')))){
			filterSpecCached <- readRDS(file.path(writeToFolder,paste0(varName,'-desiredFilterSpec.RDS')))
		}
	} else if(file.exists(file.path(writeToFolder,paste0(varName,'-filterSpec.RDS')))){
		filterSpecCached <- readRDS(file.path(writeToFolder,paste0(varName,'-filterSpec.RDS')))
	} else {
		filterSpecCached <- NULL
	}
	if(filterSpecsAreEqual(filterSpec,filterSpecCached)){
		if(verbosity>0){cat(sprintf('Valid cached filter results exist for  %s...',varName))}
		if(returnPolIDsToDrop){
			if(verbosity>0){cat('reading...')}
			if(useDesiredFilterSpec){
				polIDsToDrop <- readRDS(file.path(writeToFolder,paste0(varName,'-desiredFilter.RDS')))
			} else {
				polIDsToDrop <- readRDS(file.path(writeToFolder,paste0(varName,'-filter.RDS')))
			}
			if(verbosity>0){cat('done\n')}
			return(polIDsToDrop)
		} else {
			if(verbosity>0){cat('done\n')}
			return()
		}
	} else {
		if(verbosity>0){cat(sprintf('Invalid cached filter results exist for  %s ignoring',varName))}
		filterSpecCached <- NULL
	}
			
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
		# years
		years <- colnames(varDat)[-c(1,2)]
		filterFun <- function(year.i){
			year <- years[year.i]
			if(filterSpec[[varName]]$type == 'ltabs'){ # absolute value less than filter val
				polIDsToDrop <- varDat$polID[abs(varDat[[year]])<filterSpec[[varName]]$level]
			} else if(filterSpec[[varName]]$type == 'gtabs'){ # absolute value greater than filter val
				polIDsToDrop <- varDat$polID[abs(varDat[[year]])>filterSpec[[varName]]$level]
			} else if(filterSpec[[varName]]$type == 'sltabs'){ # absolute value less than filter val in specified SOW
				polIDsToDrop <- varDat$polID[varDat$sowID %in% filterSpec[[varName]]$sowID &
																		 	abs(varDat[[year]])<filterSpec[[varName]]$level]
			} else if(filterSpec[[varName]]$type == 'sgtabs'){ # absolute value greater than filter val in specified SOW
				polIDsToDrop <- varDat$polID[varDat$sowID %in% filterSpec[[varName]]$sowID &
																		 	abs(varDat[[year]])>filterSpec[[varName]]$level]
			}	else if (filterSpec[[varName]]$type == 'ltval'){ # value less than filter val
				polIDsToDrop <- varDat$polID[varDat[[year]]<filterSpec[[varName]]$level]
			} else if (filterSpec[[varName]]$type == 'gtval'){ # value greater than filter val
				polIDsToDrop <- varDat$polID[varDat[[year]]>filterSpec[[varName]]$level]
			} else if (filterSpec[[varName]]$type == 'sltval'){ # value less than filter val in specified SOW
				polIDsToDrop <- varDat$polID[varDat$sowID %in% filterSpec[[varName]]$sowID & 
																		 	varDat[[year]]<filterSpec[[varName]]$level]
			} else if (filterSpec[[varName]]$type == 'sgtval'){ # value greater than filter val in specified SOW
				polIDsToDrop <- varDat$polID[varDat$sowID %in% filterSpec[[varName]]$sowID &
																		 	varDat[[year]]>filterSpec[[varName]]$level]
			} else {
				stop('unkown filter spec\n')
			}
			# if a number of SOW is defined in the filter we only drop if the filter is violated in
			# at least numSOW cases
			if(!is.null(filterSpec[[varName]]$allowedTransgressions)){ 
				dropCounts <- table(polIDsToDrop)
				polIDsToDrop <- as.numeric(names(which(dropCounts>=filterSpec[[varName]]$allowedTransgressions)))
			}
			polIDsToDrop <- unique(polIDsToDrop[!is.na(polIDsToDrop)])
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
		dir.create(writeToFolder,showWarnings = F,recursive = T)
		if(useDesiredFilterSpec){
			saveRDS(polIDsToDrop,file.path(writeToFolder,paste0(varName,'-desiredFilter.RDS')))
		} else {
			saveRDS(polIDsToDrop,file.path(writeToFolder,paste0(varName,'-filter.RDS')))
		}
		if(useDesiredFilterSpec){
			saveRDS(filterSpec[[varName]],file.path(writeToFolder,paste0(varName,'-desiredFilterSpec.RDS')))
		} else {
			saveRDS(filterSpec[[varName]],file.path(writeToFolder,paste0(varName,'-filterSpec.RDS')))
		}
	} else {
		cat(sprintf('Var name %s not in filter spec.\n',varName))
	}
	if(verbosity>0){cat('done\n')}
	if(returnPolIDsToDrop){
		return(polIDsToDrop)
	}
}
filterSpecsAreEqual <- function(filterSpec1,filterSpec2){
	if(length(filterSpec1)!=length(filterSpec2)){
		return(F)
	}
	for(entry in names(filterSpec1)){
		if(length(filterSpec1[[entry]])!=length(filterSpec2[[entry]])){
			return(F)
		}
		for(subentry.i in 1:length(filterSpec1[[entry]])){
			if(filterSpec1[[entry]][subentry.i]!=filterSpec2[[entry]][subentry.i]){
				return(F)
			}
		}
	}
	return(T)
}
