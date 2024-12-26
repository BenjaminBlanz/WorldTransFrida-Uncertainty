#
# Functions to run frida
#

suppressPackageStartupMessages(require(Rmpfr)) # use to calculate the likelihood from loglikelihood

# write firda export vars ####
writeFRIDAExportSpec <- function(varsForExport.fridaNames){
	sink(file=file.path(location.frida,'Data',name.fridaExportVarsFile))
	cat(paste0(varsForExport.fridaNames,collapse='\n'))
	sink()
}

# write frida input ####
# uses location.frida and name.fridaInputFile from the global env.
writeFRIDAInput <- function(variables,values){
	parmValues <- data.frame(Variable=variables,Value=values)
	write.table(parmValues,file = file.path(location.frida,'Data',name.fridaInputFile),
							row.names = F,col.names = F,sep=',')
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

# funStretchSamplePoints ####
funStretchSamplePoints <- function(samplePoints,sampleParms,restretchSamplePoints=F){
	samplePoints <- t(samplePoints)
	if(!restretchSamplePoints){
		# Substract the min and multiply by max-min to strecth the unit interval to the
		# actual sampling range.
		samplePointsStretched <- samplePoints*(sampleParms$Max-sampleParms$Min) + sampleParms$Min
		# plot(samplePointsStretched[1,],samplePointsStretched[2,])
		# abline(v=sampleParms$Value[1],h=sampleParms$Value[2],col='red')
		samplePoints <- samplePointsStretched
		rm(samplePointsStretched)
	} else {
		# stretch the sample points to be left and right of the mean centre value of the
		# description file
		lowIdc <- samplePoints<0.5
		highIdc <- samplePoints>=0.5
		samplePointsLow <- samplePointsHigh <- samplePoints
		samplePointsLow[highIdc] <- NA
		samplePointsLow <- samplePointsLow*2*(sampleParms$Value-sampleParms$Min) + sampleParms$Min
		samplePointsHigh[lowIdc] <- NA
		samplePointsHigh <- (samplePoints-0.5)*2*(sampleParms$Max-sampleParms$Value) + sampleParms$Value
		samplePointsReStretched <- samplePoints
		samplePointsReStretched[lowIdc] <- samplePointsLow[lowIdc]
		samplePointsReStretched[highIdc] <- samplePointsHigh[highIdc]
		# plot(samplePointsReStretched[1,],samplePointsReStretched[2,])
		# abline(v=sampleParms$Value[1],h=sampleParms$Value[2],col='red')
		samplePoints <- samplePointsReStretched
		rm(samplePointsHigh,samplePointsLow,samplePointsReStretched,lowIdc,highIdc)
	}
	# back to vars in cols and samples in rows
	samplePoints<- t(samplePoints)
	rownames(samplePoints) <- 1:nrow(samplePoints)
	colnames(samplePoints) <- sampleParms[,1]
	return(samplePoints)
}


