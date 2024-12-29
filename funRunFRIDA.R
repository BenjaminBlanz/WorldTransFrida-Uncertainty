#
# Functions to run frida
#

suppressPackageStartupMessages(require(Rmpfr)) # use to calculate the likelihood from loglikelihood

# prepareSampleParms ####
prepareSampleParms <- function(excludeNames=c()){
	cat('Specify sampling parameters...')
	# read in the parameters in frida that have ranges defined
	frida_info <- read.csv("frida_info.csv")
	columnsThatAreFlags <- c(2,3,4,5,6,7,8,9,10)
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
	# deal with manually excluded items
	if(!redoAllCalc && file.exists('parExclusionList.csv') && file.size('parExclusionList.csv')>0){
		manParExclusionList <- read.csv('parExclusionList.csv')
		excludeNames <- unique(c(excludeNames,manParExclusionList$excludedName))
	}
	# deal with excluded Names
	excludedIdc <- which(sampleParms$Variable %in% excludeNames)
	if(length(c(invalidLines,excludedIdc))>0){
		excludeIdc <- unique(c(invalidLines,excludedIdc))
		sampleParms <- sampleParms[-excludedIdc,]
	}
	# add the climate case
	if(!'Climate Units.selected climate case'%in%excludeNames){
		newIdx <- nrow(sampleParms)+1
		sampleParms[newIdx,c('Variable')] <- 'Climate Units.selected climate case'
		sampleParms[newIdx,c('Value','Min','Max')] <- 
			c(37,0.5,100.5-.Machine$double.eps)
		sampleParmsRowNames <- rownames(sampleParms)
		sampleParmsRowNames[newIdx] <- as.character(as.numeric(sampleParmsRowNames[newIdx-1])+1)
		rownames(sampleParms) <- sampleParmsRowNames
	}
	cat('done\n')
	return(sampleParms)
}

# write firda export vars ####
writeFRIDAExportSpec <- function(varsForExport.fridaNames,location.frida){
	sink(file=file.path(location.frida,'Data',name.fridaExportVarsFile))
	cat(paste0(varsForExport.fridaNames,collapse='\n'))
	sink()
}

# write frida input ####
# uses location.frida and name.fridaInputFile from the global env.
writeFRIDAInput <- function(variables,values){
	if(disk.free(location.frida)< 2e4){
		stop('less than 20mib in frida location\n')
	}
	parmValues <- data.frame(Variable=variables,Value=values)
	write.table(parmValues,file = file.path(location.frida,'Data',name.fridaInputFile),
							row.names = F,col.names = F,sep=',')
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
runFridaParmsByIndex <- function(runid){
	retlist <- vector(mode = "list", length = length(runid))
	for(i in runid){
		if(i <= nrow(samplePoints)){
			writeFRIDAInput(sampleParms$Variable,samplePoints[i,])
			system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',
									 file.path(location.frida,'FRIDA.stmx')),
						 ignore.stdout = T,ignore.stderr = T,wait = T)
			runDat <- read.csv(file.path(location.frida,'Data',name.fridaOutputFile))
			colnames(runDat) <- cleanNames(colnames(runDat))
			rownames(runDat) <- runDat$year
			runDat <- runDat[,-1]
			resDat <- calDat-runDat[1:nrow(calDat),]
			logLike <- funLogLikelihood(resDat[complete.cases(resDat),],resSigma)
			like <- exp(logLike)
			# If the logLike is not NA but the run did not complete assign 
			# lowest nonzero value. We use this when narrowing the parms space
			if(is.na(runDat[[1]][nrow(runDat)])||like==0){
				like <- (sum(!is.na(runDat[[1]]))*.Machine$double.eps)
			}
			retlist[[i]] <- (list(parmsIndex=as.numeric(row.names(samplePoints)[i]),
														runDat=runDat,
														logLike=logLike,
														like=like))
		}
	}
	return(retlist)
}
# runFridaParmsBySamplePoints ####
# the same as above, but runs all samples in samplePoints for pre allocated
# samplePoints per worker.
runFridaParmsBySamplePoints <- function(saveOutPutDontReturn=FALSE,workUnit=NULL){
	retlist <- runFridaParmsByIndex(1:nrow(samplePoints))
	if(saveOutPutDontReturn){
		saveRDS(retlist,file.path(baseWD,location.output,
															paste0('workUnit-',workUnit,'-',workerID,'.RDS')))
	} else {
		return(retlist)
	}
}

# runFridaDefaultParms ####
# Uses location.frida, and name.fridaInputFile
# from the global environment
runFridaDefaultParms <- function(){
	frida_info <- read.csv("frida_info.csv")
	parVect <- frida_info$Value
	names(parVect) <- frida_info$Variable
	return(runFRIDASpecParms(parVect))
}

# runFRIDASpecParms ####
runFRIDASpecParms <- function(parVect){
	if(is.null(names(parVect))){
		stop('need names in parVect to write FRIDA input\n')
	}
	writeFRIDAInput(names(parVect),parVect)
	system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',
							 file.path(location.frida,'FRIDA.stmx')),
				 ignore.stdout = T,ignore.stderr = T,wait = T)
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
			 gsub('_$','',
			 		 gsub('_1','',
			 		 		 gsub('\\[1]','_',
			 		 		 		 gsub('_+','_',
			 		 		 		 		 gsub('[. ]','_',
			 		 		 		 		 		 tolower(colNames)))))))
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
funLogLikelihood <- function(resid,covmat){
	if(treatVarsAsIndep){
		resid[is.na(resid)] <- 0
	}
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
																					 plotDatWhileRunning=F){
	cat('cluster run...\n')
	dir.create(location.output,showWarnings = F,recursive = T)
	numSample <- nrow(samplePoints)
	logLike <- rep(NA,numSample)
	like <- rep(NA,numSample)
	names(logLike) <- 1:numSample
	names(like) <- 1:numSample
	numWorkers <- length(cl)
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
	clusterExport(cl,list('location.output','baseWD','sampleParms',
												'chunkSizePerWorker','runFridaParmsBySamplePoints',
												'calDat','resSigma',
												'runFridaParmsByIndex'))
	# plot setup 
	if(plotDatWhileRunning){
		subPlotLocations <- funPlotDat(calDat,calDat.impExtrValue,yaxPad = yaxPad,
																	 shadowIncompleteYears=F)
	}
	# running
	cat(sprintf('  Run of %i runs split up into %i work units.\n',
							numSample,length(workUnitBoundaries)-1))
	chunkTimes <- c()
	completeRunsSoFar <- 0
	for(i in 1:(length(workUnitBoundaries)-1)){
		if(!redoAllCalc && file.exists(file.path(location.output,paste0('workUnit-',i,'.RDS')))){
			cat(sprintf('\r    Reading existing unit %i',i))
			tryCatch({parOutput <- readRDS(file.path(location.output,paste0('workUnit-',i,'.RDS')))},
							 error = function(e){},warning=function(w){})
		} 
		if(!exists('parOutput')){
			cat(sprintf('\r(r) Running unit %i: samples %i to %i. ',
									i, workUnitBoundaries[i],workUnitBoundaries[i+1]-1))
			if(length(chunkTimes>1)){
				cat(sprintf('So far: Complete runs rate %.2f%%'),
						completeRunsSoFar/(workUnitBoundaries[i+1]-1))
				cat(sprintf(', time per unit %i s (%.2f r/s, %.2f r/s/thread), expect completion in %i sec',
										round(mean(chunkTimes,na.rm=T)),
										length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
										chunkSizePerWorker/mean(chunkTimes,na.rm=T),
										round(mean(chunkTimes,na.rm=T))*(length(workUnitBoundaries)-i)))
			}
			tic()
			workUnit <- workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)
			workerWorkUnits <- chunk(workUnit,numWorkers)
			# write the samplePoints of the work units to the worker directories
			for(w.i in workers){
				if(w.i <= length(workerWorkUnits) && !is.null(workerWorkUnits[[w.i]])){
					saveRDS(samplePoints[workerWorkUnits[[w.i]],],
									file.path('workerDirs',paste0(workDirBasename,w.i),'samplePoints.RDS'))
				} else {
					saveRDS(samplePoints[c(),],
									file.path('workerDirs',paste0(workDirBasename,w.i),'samplePoints.RDS'))
				}
			}
			gobble <- clusterEvalQ(cl,{
				samplePoints <- readRDS('samplePoints.RDS')
			})
			parOutput <- unlist(
				clusterEvalQ(cl,runFridaParmsBySamplePoints()),
				recursive = F)
			timing <- toc(quiet=T)
			chunkTimes[i] <- timing$toc-timing$tic
			cat('\r(s)')
			saveRDS(parOutput,file.path(location.output,paste0('workUnit-',i,'.RDS')))
			cat('\r   ')
		}
		cat('\r(l)')
		chun
		for(l in 1:length(parOutput)){
			logLike[parOutput[[l]]$parmsIndex] <- parOutput[[l]]$logLike
			like[parOutput[[l]]$parmsIndex] <- parOutput[[l]]$like
			if(is.na(parOutput[[l]]$runDat[[1]][nrow(parOutput[[l]]$runDat[[1]])])){
				completeRunsSoFar <- completeRunsSoFar + 1
			}
		}
		cat('\r   ')
		if(plotDatWhileRunning){
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
					lines(rownames(parOutput[[l]]$runDat),parOutput[[l]]$runDat[[dat.i]],
								col=i)
					# col=alpha(i,min(0.01,max(1,0.01*parOutput[[l]]$like/.Machine$double.eps))))
				}
				cat('\r   ')
			}
		}
		rm(parOutput)
	}
	if(length(chunkTimes)==0){
		cat('\r    all runs read, no calculation necessary.                                \n')
	} else {
		cat(sprintf('\r    runs completed average chunk time %i sec (%.2f r/s, %.2f r/s/thread), over all run time %i sec %s\n',
								round(mean(chunkTimes,na.rm=T)),
								length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								round(sum(chunkTimes,na.rm=T)),
								'                                                                             '))
	}
	return(list(logLike=logLike,like=like))
}

