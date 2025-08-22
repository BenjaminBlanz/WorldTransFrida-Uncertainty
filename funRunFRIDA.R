#
# Functions to run frida
#

suppressPackageStartupMessages(require(Rmpfr)) # use to calculate the likelihood from loglikelihood

# prepareSampleParms ####
prepareSampleParms <- function(excludeNames=c(),sampleParms=NULL,integerParms=NULL){
	cat('Specify sampling parameters...')
	if(is.null(sampleParms)){
		cat('reading frida_info...')
		# read in the parameters in frida that have ranges defined
		frida_info <- read.csv(file.path(location.frida.info,name.frida_info))
		columnsThatAreFlags <- c(2,3,4,5,6,7,8,9,10,11)
		# select the parameters to be sampled
		sampleParms <- frida_info[rowSums(frida_info[,columnsThatAreFlags])>0 &
																frida_info$No.Sensi==0 &
																frida_info$Policy==0,
															-columnsThatAreFlags]
		invalidLines <- which(!((sampleParms$Max-sampleParms$Min)>0 &
															sampleParms$Min <= sampleParms$Value &
															sampleParms$Value <= sampleParms$Max))
		suppressWarnings(file.remove('frida_info_errorCases.csv'))
		if(length(invalidLines)>0){
			cat('invalid lines detected, see frida_info_errorCases.csv...')
		}
		write.csv(sampleParms[invalidLines,],'frida_info_errorCases.csv')
	} else {
		invalidLines <- c()
	}
	# deal with manually excluded items
	if(!redoAllCalc &&
		 file.exists(file.path(location.frida.info,name.frida_parameter_exclusion_list)) &&
		 file.size(file.path(location.frida.info,name.frida_parameter_exclusion_list))>0){
		manParExclusionList <- read.csv(file.path(location.frida.info,name.frida_parameter_exclusion_list))
		excludeNames <- unique(c(excludeNames,manParExclusionList$excludedName))
	}
	# deal with excluded Names
	excludedIdc <- which(sampleParms$Variable %in% excludeNames)
	if(length(c(invalidLines,excludedIdc))>0){
		cat('excluding invalid and excluded parms...')
		excludeIdc <- unique(c(invalidLines,excludedIdc))
		sampleParms <- sampleParms[-excludedIdc,]
	}
	# add the integer vars e.g. climate case
	sampleParms$isInteger <- rep(FALSE,nrow(sampleParms))
	if(!is.null(integerParms)){
		cat('adding integer parms...')
		for(p.i in 1:nrow(integerParms)){
			if(!integerParms$Variable[p.i]%in%excludeNames &&
				 !integerParms$Variable[p.i]%in%sampleParms$Variable){
				newIdx <- nrow(sampleParms)+1
				sampleParms[newIdx,c('Variable')] <- integerParms$Variable[p.i] # e.g. 'Climate Units.selected climate case'
				sampleParms[newIdx,c('Value','Min','Max')] <- 
					c(integerParms$Value[p.i],integerParms$Min[p.i],integerParms$Max[p.i])
				sampleParms[newIdx,c('isInteger')] <- TRUE
				sampleParmsRowNames <- rownames(sampleParms)
				if(nrow(sampleParms)>1){
					sampleParmsRowNames[newIdx] <- as.character(as.numeric(sampleParmsRowNames[newIdx-1])+1)
				} else {
					sampleParmsRowNames[newIdx] <- 1
				}
				rownames(sampleParms) <- sampleParmsRowNames
			}
		}
	}
	cat('done\n')
	return(sampleParms)
}

# write firda export vars ####
writeFRIDAExportSpec <- function(varsForExport.fridaNames,location.frida){
	varsForExport.cleanNames <- cleanNames(varsForExport.fridaNames)
	dupe.lst <- split(seq_along(varsForExport.cleanNames), varsForExport.cleanNames)
	nonDupeIdc <- c()
	for(i in 1:length(dupe.lst)){
		nonDupeIdc[i] <- dupe.lst[[i]][1]
	}
	nonDupeIdc <- sort(nonDupeIdc)
	varsForExport.fridaNames <- varsForExport.fridaNames[nonDupeIdc]
	sink(file=file.path(location.frida,'Data',name.fridaExportVarsFile))
	cat(paste0(varsForExport.fridaNames,collapse='\n'))
	sink()
}

# write frida input ####
# uses location.frida and name.fridaInputFile from the global env.
writeFRIDAInput <- function(variables,values,policyMode=F){
	if(disk.free(location.frida)< 2e4){
		stop('less than 20mib in frida location\n')
	}
	if(policyMode){
		sink(file.path(location.frida,'Data',name.fridaInputFile))
		for(domID in names(values)){
			if(!is.na(values[domID])){
				dplID <- values[domID]
				sdmID <- jointPolicies$sdmID[jointPolicies$domID==domID & 
																		 	jointPolicies$dplID==dplID]
				pfID <- pdpMeta$pfID[pdpMeta$domID==domID & pdpMeta$sdmID==sdmID]
				sdpID <- jointPolicies$sdpID[jointPolicies$domID==domID & jointPolicies$dplID==dplID]
				pf <- pdp.lst[[pfID]]
				cat(paste0(pf[which(pf$polID==sdpID),c('policyString')],'\n',collapse = ''))
			}
		}
		sink()
	} else {
		parmValues <- data.frame(Variable=unlist(variables),Value=unlist(values))
		write.table(parmValues,file = file.path(location.frida,'Data',name.fridaInputFile),
								row.names = F,col.names = F,sep=',')
	}
}

# check for free space
disk.free <- function(path = getwd()) {
	if(length(system("which df", intern = TRUE))) {
		cmd <- sprintf("df %s", path)
		exec <- system(cmd, intern = TRUE)
		exec <- strsplit(exec[length(exec)], "[ ]+")[[1]]
		exec <- as.numeric(exec[4])
		structure(exec, names = c("available"))
	} else {
		stop("'df' command not found")
	}
}

# runFridaParmsByIndex ####
# Uses from global env:
#   sampleParms,samplePoints,location.frida, and name.fridaInputFile
# If retNegLogLike also uses from global env:
# 	calDat,resSigma
runFridaParmsByIndex <- function(runid,silent=T,policyMode=F){
	retlist <- vector(mode = "list", length = length(runid))
	cat('\n')
	for(i in runid){
		cat(paste('\r',i))
		if(i <= nrow(samplePoints)){
			sink(file.path(location.frida,'lastRun.txt'))
			cat(row.names(samplePoints)[i],'\n')
			sink()
			writeFRIDAInput(colnames(samplePoints),samplePoints[i,],policyMode=policyMode)
			system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',
									 file.path(location.frida,'FRIDA.stmx')),
						 ignore.stdout = silent,ignore.stderr = silent,wait = T)
			runDat <- read.csv(file.path(location.frida,'Data',name.fridaOutputFile))
			origColNames <- unname(unlist(read.table(file.path(location.frida,'Data',name.fridaOutputFile),
															 sep=',')[1,]))[-1]
			colnames(runDat) <- cleanNames(colnames(runDat))
			# catch failed runs causing NAs in year variable and crash in the rownames assignment
			runDat$year <- seq(runDat$year[1],length.out=nrow(runDat))
			rownames(runDat) <- runDat$year
			runDat <- runDat[,-1]
			if(!policyMode){
				resDat <- calDat-runDat[1:nrow(calDat),colnames(calDat)]
				logLike <- funLogLikelihood(resDat,resSigma)
				# If the logLike is not NA but the run did not complete assign 
				# lowest nonzero value. We use this when narrowing the parms space
				if(is.na(runDat[[1]][nrow(runDat)])||logLike==-Inf){
					logLike <- -.Machine$double.xmax+sum(!is.na(runDat[[1]]))*.Machine$double.eps
				}
			} else {
				logLike <- -.Machine$double.xmax+sum(!is.na(runDat[[1]]))*.Machine$double.eps
			}
			suppressWarnings(parmsIndex<-as.numeric(row.names(samplePoints)[i]))
			if(is.na(parmsIndex)){
				retlist[[i]] <- (list(parmsName=row.names(samplePoints)[i],
															parmsIndex=i,
															runDat=runDat,
															origColNames=origColNames,
															logLike=logLike))
			} else {
				retlist[[i]] <- (list(parmsIndex=parmsIndex,
															runDat=runDat,
															origColNames=origColNames,
															logLike=logLike))
			}
		}
	}
	return(retlist)
}
# runFridaParmsBySamplePoints ####
# the same as above, but runs all samples in samplePoints for pre allocated
# samplePoints per worker.
runFridaParmsBySamplePoints <- function(policyMode=F){
	retlist <- runFridaParmsByIndex(1:nrow(samplePoints),policyMode=policyMode)
	if(writePerWorkerFiles){
		workerID <- ifelse(exists('workerID'),workerID,0)
		workUnit.i <- ifelse(exists('workUnit.i'),workUnit.i,0)
		logLike.df <- saveParOutputToPerVarFiles(parOutput = retlist,workUnit.i = workUnit.i,
														 workerID = workerID)
		if(doNotReturnRunDataSavePerWorkerOnly){
			retlist <- list()
		}
		retlist[['logLike.df']] <- logLike.df
	}
	return(retlist)
}

# runFridaDefaultParms ####
# Uses location.frida, and name.fridaInputFile
# from the global environment
runFridaDefaultParms <- function(silent=T){
	return(runFRIDASpecParms(c(),silent))
}

# runFRIDASpecParms ####
runFRIDASpecParms <- function(parVect,silent=T){
	if(is.null(names(parVect))&length(parVect)>0){
		stop('need names in parVect to write FRIDA input\n')
	}
	if(length(parVect)==0){
		writeFRIDAInput(c(),c())
	}else{
		writeFRIDAInput(names(parVect),parVect)
	}
	system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',#'-s', #to output isdb
							 file.path(location.frida,'FRIDA.stmx')),
				 ignore.stdout = silent,ignore.stderr = silent,wait = T)
	runDat <- read.csv(file.path(location.frida,'Data',name.fridaOutputFile))
	colnames(runDat) <- cleanNames(colnames(runDat))
	rownames(runDat) <- runDat$year
	runDat <- runDat[,-1]
	return(runDat)
}

# cleanNames ####
# takes a vector of e.g. column names and brings them into 
# a comparable standard format
# also drops the trailing 1 of the run id which we do not use
# as we only have single run setups of frida
cleanNames <- function(colNames){
	gsub('time','year',
			 gsub('_+$','',
			 		 gsub('_+','_',
			 		 		 gsub(',','_',
			 		 		 		 gsub('\\$','',
			 		 		 		 		 gsub('_1','',
			 		 		 		 		 		 gsub('\\]','_',
			 		 		 		 		 		 		 gsub('\\[\\*','_',
				 		 		 		 		 		 		 gsub('\\[\\d+','_',
				 		 		 		 		 		 		 		 gsub('[. ]','_',
				 		 		 		 		 		 		 		 		 tolower(colNames)))))))))))
}

# idxOfVarName ####
idxOfVarName <- function(varNames,vecOfVarNames){
	varNames <- cleanNames(varNames)
	vecOfVarNames <- cleanNames(vecOfVarNames)
	idx <- c()
	for(i in 1:length(varNames)){
		idx[i] <- which(vecOfVarNames==varNames[i])
	}
	return(idx)
}

# taken from funchir package
# convert plot region sizes to inches
xdev2in <- function (x = 1) 
{
	x * par("pin")[1L]/diff(par("usr")[1L:2L])
}
ydev2in <- function (y = 1) 
{
	y * par("pin")[2L]/diff(par("usr")[3L:4L])
}
xydev2in <-	function (xy = 1) 
{
	u = par("usr")
	xy * par("pin")/c(u[2L] - u[1L], u[4L] - u[3L])
}
dist.f <- function(oi,yi,ys,offsets,keepOutSize){
	return(min(sqrt((ydev2in(ys-yi))^2+(xdev2in(offsets-oi))^2))-keepOutSize)
}

# funValidRange ####
# returns the first and last index in the variable that has a data point
funValidRange <- function(x){
	validRange <- c(1,length(x))
	if(is.na(x[1])){
		validRange[1] <- which(diff(is.na(x))==-1)[1]+1
	}
	if(is.na(x[length(x)])){
		validRange[2] <- max(which(diff(is.na(x))==1))
	}
	return(validRange)
}

# funLogLikelihood ####
# dmvnorm function from the mvtnorm package for reference
# dmvnorm <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE, 
# 										 checkSymmetry = TRUE) 
# {
# 	if (is.vector(x)) 
# 		x <- matrix(x, ncol = length(x))
# 	p <- ncol(x)
# 	if (!missing(mean)) {
# 		if (!is.null(dim(mean))) 
# 			dim(mean) <- NULL
# 		if (length(mean) != p) 
# 			stop("x and mean have non-conforming size")
# 	}
# 	if (!missing(sigma)) {
# 		if (p != ncol(sigma)) 
# 			stop("x and sigma have non-conforming size")
# 		if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
# 																			check.attributes = FALSE)) 
# 			stop("sigma must be a symmetric matrix")
# 	}
# 	dec <- tryCatch(base::chol(sigma), error = function(e) e)
# 	if (inherits(dec, "error")) {
# 		x.is.mu <- colSums(t(x) != mean) == 0
# 		logretval <- rep.int(-Inf, nrow(x))
# 		logretval[x.is.mu] <- Inf
# 	}
# 	else {
# 		tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
# 		rss <- colSums(tmp^2)
# 		logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * 
# 																												pi) - 0.5 * rss
# 	}
# 	names(logretval) <- rownames(x)
# 	if (log) 
# 		logretval
# 	else exp(logretval)
# }
funLogLikelihood <- function(resid,covmat,treatVarsAsIndep=.GlobalEnv$treatVarsAsIndep){
	if(treatVarsAsIndep){
		perVarLogLike <- rep(NA,ncol(resid))
		for(i in 1:ncol(resid)){
			perVarLogLike[i] <- sum(dnorm(resid[!is.na(resid[,i]),i], 0, sqrt(covmat[i,i]),log = T))
			# treat all non observed values as likely as the average observed one
			# perVarLogLike[i] <- perVarLogLike[i]+sum(is.na(resid[,i]))*perVarLogLike[i]/sum(!is.na(resid[,i]))
		}
		return(sum(perVarLogLike))
	} else {
		if(nrow(resid)==0||sum(is.na(resid))>0){
			return(-Inf)
		}	else {
			logLike <- sum(mvtnorm::dmvnorm(resid,rep(0,ncol(resid)),covmat,log=T,checkSymmetry = F))
			if(is.na(logLike)){
				return(-Inf)
			} else {
				return(logLike)
			}
		}
	}
}

# chunk ####
# cuts a vector into n equal parts
chunk <- function(x,n){
	split(x, cut(seq_along(x), n, labels = FALSE)) 
}

# cluster run ####
clusterRunFridaForSamplePoints <- function(samplePoints,chunkSizePerWorker,
																					 location.output,
																					 calDat,resSigma,
																					 redoAllCalc=F,
																					 plotDatWhileRunning=F,
																					 plotDatPerChunWhileRunning=F,
																					 plotPerChunk=T,
																					 yaxPad=0.2,
																					 baseLL=-29567.06,
																					 skipRunJustRead=F,
																					 numWorkers=length(cl),
																					 calDat.impExtrValue=NULL,
																					 plotMinAlpha = 0.01,
																					 plotBaseName = 'runs-'){
	cat('cluster run...\n')
	# If we are not plotting while running we can directly store from the running
	# worker processes. This is MUCH faster.
	if((plotDatWhileRunning|plotDatPerChunWhileRunning)&writePerWorkerFiles){
		cat('Note: Forcing writePerWorkerFiles to false, as plotting is enabled\n')
		writePerWorkerFiles <<- F
	}
	# prevent skipping work if parOutput exists in the global env.
	if(exists('parOutput')){rm(parOutput)}
	numSample <- nrow(samplePoints)
	logLike <- rep(NA,numSample)
	names(logLike) <- 1:numSample
	workUnitBoundaries <- seq(1,numSample,chunkSizePerWorker*numWorkers)
	# in case the chunkSize is not a perfect divisor of the numSample, add numSample as the 
	# final boundary
	if(workUnitBoundaries[length(workUnitBoundaries)]!=numSample){
		workUnitBoundaries <- c(workUnitBoundaries,numSample)
	}
	# add one to the last work unit boundary, as during running we always deduct one from the next boundary
	workUnitBoundaries[length(workUnitBoundaries)] <- numSample+1
	
	# initialise  cluster
	baseWD <- getwd()
	if(!skipRunJustRead){
		clusterExport(cl,list('location.output','baseWD','sampleParms',
													'chunkSizePerWorker','runFridaParmsBySamplePoints',
													'calDat','resSigma',
													'runFridaParmsByIndex','writePerWorkerFiles'))
	}
	# plot setup 
	if(plotDatWhileRunning & !plotPerChunk){
		if(!(plotCape['X11']|plotCape['aqua'])){
			ncols <- ncol(calDat)
			sqrtNcols <- sqrt(ncols)
			plotCols <- round(sqrtNcols)
			plotRows <- ceiling(sqrtNcols)
			png(file.path(baseWD,location.output,paste0(plotBaseName,'all.png')),
					width=5*plotCols,
					height=5*plotRows+5/4,
					unit='cm',res=150)
		}
		subPlotLocations <- funPlotDat(calDat,calDat.impExtrValue,yaxPad = yaxPad,
																	 shadowIncompleteYears=F)
	}
	# running
	cat(sprintf('  Run of %i runs split up into %i work units of size %i (%i per worker).\n',
							numSample,length(workUnitBoundaries)-1,chunkSizePerWorker*numWorkers,chunkSizePerWorker))
	chunkTimes <- c()
	completeRunsSoFar <- 0
	i <- 0
	while(i<(length(workUnitBoundaries)-1)){
		i <- i+1
		workUnit.i <- i
		clusterExport(cl,list('workUnit.i'),envir=environment())
		if(!redoAllCalc && file.exists(file.path(baseWD,location.output,paste0('workUnit-',i,'.RDS')))){
			cat(sprintf('\r(r) Using existing unit %i',i))
			tryCatch({parOutput <- readRDS(file.path(baseWD,location.output,paste0('workUnit-',i,'.RDS')))},
							 error = function(e){},warning=function(w){})
			if(exists('parOutput')){
				if(length(parOutput)>(workUnitBoundaries[i+1]-workUnitBoundaries[i])){
					lastChunkSize <- length(parOutput)
					cat(sprintf(', existing output has different chunkSize (%i rather than %i), resorting remaining work',
											lastChunkSize,
											chunkSizePerWorker*numWorkers))
					workUnitBoundaries[i+1] <- workUnitBoundaries[i]+lastChunkSize
					workUnitBoundaries <- c(workUnitBoundaries[1:i],
																	seq(workUnitBoundaries[i+1],numSample,chunkSizePerWorker*numWorkers))
					if(workUnitBoundaries[length(workUnitBoundaries)]!=numSample){
						workUnitBoundaries <- c(workUnitBoundaries,numSample)
					}
					workUnitBoundaries[length(workUnitBoundaries)] <- numSample+1
				}
			}
		}
		if(!exists('parOutput')){
			if(skipRunJustRead){
				stop('Missing run files.\n')
			}
			cat(sprintf('\r(r) Running unit %i: samples %i to %i. ',
									i, workUnitBoundaries[i],workUnitBoundaries[i+1]-1))
			if((workUnitBoundaries[i]-1)>0){
				cat(sprintf('So far: Complete runs %i (%.1f%%)',
										completeRunsSoFar,100*completeRunsSoFar/(workUnitBoundaries[i]-1)))
			}
			if(length(chunkTimes>1)){
				cat(sprintf(', time per unit %i s (%.2f r/s, %.2f r/s/thread), expect completion in %s',
										round(mean(chunkTimes,na.rm=T)),
										length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
										chunkSizePerWorker/mean(chunkTimes,na.rm=T),
										if(exists('dseconds')){dseconds(round(mean(chunkTimes,na.rm=T))*
																											(length(workUnitBoundaries)-i))} 
										else{paste0(round(mean(chunkTimes,na.rm=T))*(length(workUnitBoundaries)-i),'s')}))
			}
			tic()
			workUnit <- workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)
			workerWorkUnits <- chunk(workUnit,numWorkers)
			# write the samplePoints of the work units to the worker directories
			for(w.i in workers){
				if(w.i <= length(workerWorkUnits) && !is.null(workerWorkUnits[[w.i]])){
					saveRDS(samplePoints[workerWorkUnits[[w.i]],],
									file.path(name.workDir,paste0(name.workerDirBasename,w.i),'samplePoints.RDS'))
				} else {
					saveRDS(samplePoints[c(),],
									file.path(name.workDir,paste0(name.workerDirBasename,w.i),'samplePoints.RDS'))
				}
			}
			gobble <- clusterEvalQ(cl,{
				samplePoints <- readRDS('samplePoints.RDS')
			})
			parOutput <- clusterEvalQ(cl,runFridaParmsBySamplePoints())
			timing <- toc(quiet=T)
			chunkTimes[i] <- timing$toc-timing$tic
			if(writePerWorkerFiles){
				logLike.df <- data.frame(id=integer(),logLike=double())
				for(r.i in 1:length(parOutput)){
					logLike.df <- rbind(logLike.df,parOutput[[r.i]]$logLike.df)
					parOutput[[r.i]]$logLike.df <- NULL
				}
			}
			cat('\r(s)')
			parOutput <-  unlist(parOutput, recursive = F)
			saveRDS(parOutput,file.path(baseWD,location.output,paste0('workUnit-',i,'.RDS')))
			if(!writePerWorkerFiles){
				saveParOutputToPerVarFiles(parOutput=parOutput, workUnit.i=i)
			}
			cat('\r   ')
		}
		cat('\r(l)')
		
		if(writePerWorkerFiles & exists('logLike.df')){
			logLike[logLike.df$id] <- logLike.df$logLike
			completeRunsSoFar <- sum(logLike > -.Machine$double.xmax+(200*.Machine$double.eps))
		} else {
			for(l in 1:length(parOutput)){
				logLike[parOutput[[l]]$parmsIndex] <- parOutput[[l]]$logLike 
				if(parOutput[[l]]$logLike > -.Machine$double.xmax+(200*.Machine$double.eps)){
					completeRunsSoFar <- completeRunsSoFar + 1
				}
			}
		}
		cat('\r   ')
		if(plotDatWhileRunning&&!plotDatPerChunWhileRunning){
			cat('\r(p)')
			for(dat.i in 1:ncol(calDat)){
				par(mfg = which(subPlotLocations==dat.i,arr.ind = T))
				xlims <- c(min(as.numeric(rownames(calDat)))-0.5,
									 max(as.numeric(rownames(calDat)))+0.5)
				yrange <- range(calDat[[dat.i]],na.rm=T)
				ylims <- c(yrange[1]-abs(diff(yrange))*yaxPad,
									 yrange[2]+abs(diff(yrange))*yaxPad)
				plot(rownames(calDat),calDat[[dat.i]],type='n',
						 xaxt='n',yaxt='n',
						 xaxs='i',yaxs='i',
						 xlim=xlims,
						 ylim=ylims)
				for(l in 1:length(parOutput)){
					runLL <- parOutput[[l]]$logLike
					lines(rownames(parOutput[[l]]$runDat),parOutput[[l]]$runDat[[dat.i]],
								col=adjustcolor(i,min(1,max(plotMinAlpha,
																						1/(abs(runLL-baseLL)+1)
								))))
				}
			}
			cat('\r   ')
		}
		if(plotDatPerChunWhileRunning){
			cat('\r(p)')
			ncols <- ncol(calDat)
			sqrtNcols <- sqrt(ncols)
			plotCols <- round(sqrtNcols)
			plotRows <- ceiling(sqrtNcols)
			png(file.path(baseWD,location.output,paste0(plotBaseName,'-Chunk-',i,'.png')),
					width=5*plotCols,
					height=5*plotRows+5/4,
					unit='cm',res=150)
			subPlotLocations.chk <- funPlotDat(calDat,calDat.impExtrValue,yaxPad = yaxPad,
																				 shadowIncompleteYears=F)
			for(dat.i in 1:ncol(calDat)){
				par(mfg = which(subPlotLocations.chk==dat.i,arr.ind = T))
				xlims <- c(min(as.numeric(rownames(calDat)))-0.5,
									 max(as.numeric(rownames(calDat)))+0.5)
				yrange <- range(calDat[[dat.i]],na.rm=T)
				ylims <- c(yrange[1]-abs(diff(yrange))*yaxPad,
									 yrange[2]+abs(diff(yrange))*yaxPad)
				plot(rownames(calDat),calDat[[dat.i]],type='n',
						 xaxt='n',yaxt='n',
						 xaxs='i',yaxs='i',
						 xlim=xlims,
						 ylim=ylims)
				for(l in 1:length(parOutput)){
					runLL <- parOutput[[l]]$logLike
					lines(rownames(parOutput[[l]]$runDat),parOutput[[l]]$runDat[[dat.i]],
								col=adjustcolor(1,min(1,max(plotMinAlpha,
																						1/(abs(runLL-baseLL)+1)
																						))))
				}
			}
			dev.off()
			cat('\r   ')
		}
		rm(parOutput)
	}
	
	if(plotDatWhileRunning&!plotDatPerChunWhileRunning){
		cat(sprintf('  Saving figure...'))
		plotCape <- capabilities()
		if(!(plotCape['X11']|plotCape['aqua'])){
			dev.off()
		} else{
			dev.print(png,width=5*ncol(subPlotLocations),
								height=5*(nrow(subPlotLocations)-1)+5/4,
								unit='cm',res=150,
								file.path(baseWD,location.output,paste0(plotBaseName,'all.png')))
		}
		cat('done\n')
	}
	if(length(chunkTimes)==0){
		cat(sprintf('\r    all runs read, no calculation necessary. Complete runs %i (%.1f%%)                           \n',
								completeRunsSoFar,100*completeRunsSoFar/numSample))
	} else {
		cat(sprintf('\r    complete runs %i (%.2f%%), average chunk time %i sec (%.2f r/s, %.2f r/s/thread), over all run time %s %s\n',
								completeRunsSoFar,100*completeRunsSoFar/numSample,
								round(mean(chunkTimes,na.rm=T)),
								length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								dseconds(round(sum(chunkTimes,na.rm=T))),
								'                                                                             '))
	}
	# merge all the individual per Var files into one complete one
	mergePerVarFiles()
	return(logLike)
}

loadClusterRuns <- function(location.output){
	runFilesList <- list.files(location.output,pattern = 'workUnit-[0-9]+\\.RDS')
	retList <- c()
	for(f.i in 1:length(runFilesList)){
		parOutput <- readRDS(file.path(baseWD,location.output,paste0('workUnit-',f.i,'.RDS')))
		retList <- c(retList,parOutput)
	}
	return(retList)
}

saveParOutputToPerVarFiles <- function(parOutput, workUnit.i='0', workerID='0',
																			 verbosity=0){
	if(!exists('compressCsv')){compressCsv<-T}
	varNames <- unique(parOutput[[1]]$origColNames)
	workUnitLength <- length(parOutput)
	perVarData <- list()
	logLike <- data.frame(id=rep(NA,workUnitLength),logLike=rep(NA,workUnitLength))
	varsIdc.lst <- list()
	varNamesNoSOW.all <- gsub('\\[\\d+','',varNames)
	varNamesNoSOW <- unique(varNamesNoSOW.all)
	for(v.i in 1:length(varNamesNoSOW)){
		varName <- cleanNames(varNamesNoSOW[v.i])
		varsIdc.lst[[varName]] <- which(varNamesNoSOW.all==varNamesNoSOW[v.i])
		numSOW <- length(varsIdc.lst[[cleanNames(varNamesNoSOW[v.i])]])
		if(numSOW>1){
			perVarData[[varName]] <- data.frame(matrix(NA,ncol=2+nrow(parOutput[[1]]$runDat),nrow=workUnitLength*numSOW))
			colnames(perVarData[[varName]]) <- c('polID','sowID',rownames(parOutput[[1]]$runDat))
		} else {
			perVarData[[varName]] <- data.frame(matrix(NA,ncol=1+nrow(parOutput[[1]]$runDat),nrow=workUnitLength))
			colnames(perVarData[[varName]]) <- c('id',rownames(parOutput[[1]]$runDat))
		}
	}
	varNames <- names(perVarData)
	if(verbosity>0){
		cat('Reading parOutput:\n')
	}
	for(run.i in 1:length(parOutput)){
		if(verbosity>0){
			cat(sprintf('\rReading run %i of %i.',run.i, workUnitLength))
		}
		runDat <- parOutput[[run.i]]$runDat
		for(varName in varNames){
			if(length(varsIdc.lst[[varName]])>1){
				perVarDataIndices <- (run.i+(run.i-1)*(numSOW-1)):((run.i+(run.i-1)*(numSOW-1))+numSOW-1)
				perVarData[[varName]][perVarDataIndices,'polID'] <- parOutput[[run.i]]$parmsIndex
				perVarData[[varName]][perVarDataIndices,'sowID'] <- 1:numSOW
				varDataFromRun <- unname(t(parOutput[[run.i]]$runDat[,varsIdc.lst[[varName]]]))
				# catches varDataFromRun being shorter than the output data frame
				perVarData[[varName]][perVarDataIndices,3:(2+ncol(varDataFromRun))] <- varDataFromRun
			} else{
				perVarData[[varName]][run.i,] <- unname(unlist(c(parOutput[[run.i]]$parmsIndex,runDat[varName])))
			}
		}
		logLike[run.i,] <- c(parOutput[[run.i]]$parmsIndex,parOutput[[run.i]]$logLike)
	}
	if(verbosity>0){
		cat('\nWriting to files\n')
	}
	# outpuperVarOutputTypestTypes from config
	for(outputType in perVarOutputTypes){
		for(varName in varNames){
			if(verbosity>0){
				cat(sprintf('\rWriting var %i of %i: %s %s', run.i, workUnitLength, varName, rep(' ',100)))
			}
			dir.create(file.path(baseWD,location.output,'detectedParmSpace',paste0('PerVarFiles-',outputType),varName),
								 showWarnings = F,recursive = T)
			writePerVarFile(perVarData[[varName]],
											file = file.path(baseWD,location.output,'detectedParmSpace',paste0('PerVarFiles-',outputType),varName,
																			 paste0(varName,'-',workUnit.i,'-',workerID)),
											outputType = outputType, compressCsv = compressCsv)
		}
		if(verbosity>0){
			cat(sprintf('\rWriting logLike %s', rep(' ',100)))
		}
		dir.create(file.path(baseWD,location.output,'detectedParmSpace',paste0('PerVarFiles-',outputType),'logLike'),
						 showWarnings = F,recursive = T)
		writePerVarFile(logLike,
										file.path(baseWD,location.output,'detectedParmSpace',paste0('PerVarFiles-',outputType),'logLike',
															paste0('logLike','-',workUnit.i,'-',workerID)),
										outputType = outputType, compressCsv = compressCsv)
	}
	if(verbosity>0){
		cat('\n')
	}
	return(logLike)
}

workerReadPerVarFiles <- function(i,outputType,perVarSubfolder,fileList){
	readPerVarFile(file.path(perVarSubfolder,fileList[i]))
}
workerMergePerVarFiles <- function(v.i,outputType,outputTypeFolder,varNames,verbosity=0,
																	 compressCsv=T){
	varName <- varNames[v.i]
	perVarSubfolder <- file.path(outputTypeFolder,varName)
	fileList <- list.files(perVarSubfolder)
	if(verbosity>0){cat(sprintf('Processing %i files of %s...',length(fileList),varName))}
	if(verbosity>0){cat('reading and merging...')}
	varData <- readPerVarFile(file.path(perVarSubfolder,fileList[1]),outputType)
	for(f.i in 2:length(fileList)){
		nextData <- readPerVarFile(file.path(perVarSubfolder,fileList[f.i]),outputType)
		# hack to deal with incomplete runs messing up column headers
		# proper fix is in data generation, but this will make old results work
		colnames(nextData) <- colnames(varData)
		varData <- rbind(varData,nextData)
	}
	varData <- sort_by(varData,varData[,1])
	colnames(varData) <- gsub('(^X)([0-9]{4})','\\2',colnames(varData),perl = T)
	if(verbosity>0){cat('writing...')}
	writePerVarFile(varData,file.path(outputTypeFolder,varName),
									outputType=outputType,compressCsv=compressCsv)
	if(verbosity>0){cat('removing split files...')}
	unlink(perVarSubfolder,recursive = T,force = T)
	if(verbosity>0){cat('done\n')}
	gc()
}
mergePerVarFiles <- function(verbosity=1,parStrat=2,compressCsv=T){
	if(verbosity>0){
		cat('Merging per Var files\n')
	}
	outputFolder <- file.path(baseWD,location.output,'detectedParmSpace')
	outputTypeFolders <- basename(list.dirs(outputFolder,recursive = F))
	for(outputTypeFolder in outputTypeFolders){
		if(verbosity>0){cat(paste('Entering',outputTypeFolder,'\n'))}
		outputType <- strsplit(outputTypeFolder,'-')[[1]][2]
		outputTypeFolder <- file.path(baseWD,location.output,'detectedParmSpace',outputTypeFolder)
		varNames <- basename(list.dirs(outputTypeFolder,recursive = F))
		if(verbosity>0){cat(sprintf('Found %i variable sub folder(s)\n',length(varNames)))}
		if(length(varNames)==0){
			next
		}
		if(parStrat==1){
			for(v.i in 1:length(varNames)){
				varName <- varNames[v.i]
				perVarSubfolder <- file.path(outputTypeFolder,varName)
				fileList <- list.files(perVarSubfolder)
				if(verbosity>0){cat(sprintf('(%i of %i) Processing %i files of %s...reading...',
																		v.i, length(varNames),
																		length(fileList),varName))}
				filesContents.lst <- parLapply(cl,1:length(fileList),workerReadPerVarFiles,
																			 outputType=outputType,
																			 perVarSubfolder=perVarSubfolder,
																			 fileList=fileList)
				if(verbosity>0){cat('merging...')}
				varData <- filesContents.lst[[1]]
				for(i in 2:length(filesContents.lst)){
					varData <- rbind(varData,filesContents.lst[[i]])
				}
				varData <- sort_by(varData,varData[,1])
				colnames(varData) <- gsub('(^X)([0-9]{4})','\\2',colnames(varData),perl = T)
				if(verbosity>0){cat('writing...')}
				writePerVarFile(varData,file.path(outputTypeFolder,varName),compressCsv)
				if(verbosity>0){cat('removing split files...')}
				unlink(perVarSubfolder,recursive = T,force = T)
				if(verbosity>0){cat('done\n')}
			}
		} else if(parStrat==2){
			if(verbosity>0){cat('Parallel proccessing all vars...')}
			parLapply(cl,1:length(varNames),workerMergePerVarFiles,
								outputType=outputType,
								outputTypeFolder=outputTypeFolder,
								varNames=varNames,
								compressCsv=compressCsv)
			if(verbosity>0){cat('done\n')}
		} else {
			stop('unkown parStrat\n')
		}
	}
}

readPerVarFile <- function(file,outputType=NULL){
	if(is.null(outputType)){
		outputType <- tools::file_ext(file)
		if(outputType==''){
			stop('No outputType and no file ext\n')
		}
		if(outputType=='gz'){
			outputType <- 'csv'
		}
	}
	# remove file extensions
	fileNoExt <- gsub('\\.gz','',file)
	fileNoExt <- gsub('\\.csv','',fileNoExt)
	fileNoExt <- gsub('\\.RDS','',fileNoExt)
	if(outputType=='csv'){
		if(!file.exists(paste0(fileNoExt,'.csv')) && 
			 file.exists(paste0(fileNoExt,'.csv.gz'))){
			return(read.csv(gzfile(paste0(fileNoExt,'.csv.gz'))))
		}
		return(read.csv(paste0(fileNoExt,'.csv')))
		
	} else if(outputType=='RDS'){
		return(readRDS(paste0(fileNoExt,'.RDS')))
	}
}

writePerVarFile <- function(varData,file,outputType=NULL,compressCsv=T){
	if(is.null(outputType)){
		outputType <- tools::file_ext(file)
		if(outputType==''){
			stop('No outputType and no file ext\n')
		}
		if(outputType=='gz'){
			outputType <- 'csv'
			compressCsv <- T
		}
	}
	# remove file extensions
	fileNoExt <- gsub('\\.gz','',file)
	fileNoExt <- gsub('\\.csv','',fileNoExt)
	fileNoExt <- gsub('\\.RDS','',fileNoExt)
	if(outputType=='csv'){
		write.table(varData,paste0(fileNoExt,'.csv'),
								row.names = F,sep=',')
		if(compressCsv){
			system(paste('gzip',paste0(fileNoExt,'.csv')))
		}
	} else if(outputType=='RDS'){
		saveRDS(varData,paste0(fileNoExt,'.RDS'))
	} else {
		stop('No outputType\n')
	}
}



















