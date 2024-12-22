#
#	provides a function to determine the likelihood of a given set of frida parameters
#
#
# 2024 Benjamin Blanz
# 

require(caret,quietly = T) # to find linear combinations and remove them in the calib dat
require(matrixcalc,quietly = T) # to test positive definitnes of cov matrix
require(Rmpfr,quietly = T) # for the arbitrary precis math needed in the likelihood

source('config.R')
source('funRunFRIDA.R')

# read calib data ####
cat('reading calibration data ...')
calDat.orig <- read.csv(file.path(location.frida,'Data','Calibration Data.csv'), header=FALSE,
									 colClasses = "character") 
# character class for input as in this step we do not want to change the data in any
# way

# store names of vars for export in frida format (including dots etc)
# but remove trailing spaces as stella does not like those apparently
# also remove first entry, that is time
varsForExport.fridaNames.orig <- gsub(' $','',calDat.orig$V1)[-1]
varsForExport.cleanNames.orig <- cleanNames(varsForExport.fridaNames.orig)

# kick out duplicates
if(any(duplicated(calDat.orig))){
	cat('\n  removing duplicate rows in calibration data\n ',
			paste0(calDat.orig$V1[duplicated(calDat.orig)],collapse = '\n  '))
	calDat.orig <- calDat.orig[!duplicated(calDat.orig),]
}
if(any(duplicated(calDat.orig$V1))){
	stop('\n  There are duplicated variable names with non identical data!\n ',
			 paste0(calDat.orig$V1[duplicated(calDat.orig$V1)],collapse = '\n '))
}
calDat.orig <- t(calDat.orig)
# first row of the transpose table is the column headers
colnames(calDat.orig) <- calDat.orig[1,]
# drop the first column that conatins column headers
calDat.orig <- calDat.orig[-1,]
# drop all rows after the empty row signifying the end of the data
calDat.orig <- calDat.orig[-seq(which(rowSums(calDat.orig!='')==0)[1],nrow(calDat.orig)),]
# write out the transposed and cleaned calDat. Not sure it needs to be this complicaetd
# but all other methods changed the data in situ.
sink(file.path(location.frida,'Data','Calibration Data Cleaned and Transposed.csv'))
cat(paste(paste0('"',colnames(calDat.orig),'"'),collapse=','))
cat('\n')
for(r in 1:nrow(calDat.orig)){
	cat(paste(calDat.orig[r,],collapse=','))
	cat('\n')
}
sink()
# read in the cleaned and transformed file, this time read numbers as numbers etc.
calDat.orig <- read.csv(file.path(location.frida,'Data','Calibration Data Cleaned and Transposed.csv'))
colnames(calDat.orig) <- cleanNames(colnames(calDat.orig))
rownames(calDat.orig) <- calDat.orig$year
calDat.orig <- calDat.orig[,-1]
cat('done\n')

# determine Vars to Exclude ####
cat('Applying rules to exclude variables from likelihood analysis...\n')
# exclude all that have fewer than three data points as we
# can not determine variance in those cases.
numNA <- c()
for(i in 1:ncol(calDat.orig)){
	numNA[i] <- sum(is.na(calDat.orig[,i]))
}
exclude.idc <- which((nrow(calDat.orig)-numNA) < minObsForLike)

cat(sprintf('  Excluding for less than minObsForLike (%i) observations:\n    ',minObsForLike))
cat(paste(colnames(calDat.orig)[exclude.idc],collapse='\n    '))
cat('\n')

# find linear combinations ####
if(removeLinearCombinations){
	calDat <- calDat.orig[,-exclude.idc]
	linearCombos <- caret::findLinearCombos(calDat[complete.cases(calDat),])
	if(length(linearCombos$remove)>0){
		new.exclude.idc <- idxOfVarName(colnames(calDat)[linearCombos$remove],
																		varsForExport.cleanNames.orig)
		cat('  Excluding because they are linear combinations of other vars:\n    ')
		cat(paste(varsForExport.cleanNames.orig[new.exclude.idc],collapse='\n    '))
		cat('\n')
		exclude.idc <- unique(c(exclude.idc,new.exclude.idc))
	}
}

# exclude vars ####
doneExcluding <- F
while(!doneExcluding){
	calDat <- calDat.orig[,-exclude.idc]
	varsForExport.fridaNames <- varsForExport.fridaNames.orig[-exclude.idc]
	doneExcluding <- T
	
	# write Vars to frida input file ####
	sink(file=file.path(location.frida,'Data',name.fridaExportVarsFile))
	cat(paste0(varsForExport.fridaNames,collapse='\n'))
	sink()
	
	# default run ####
	defDat <- runFridaDefaultParms()
	if(sum(colnames(defDat)!=colnames(calDat))!=0){
		stop('Mismatch in the columns of calibration data and model result data')
	}
	
	# default resid ####
	resDat <- defDat[1:nrow(calDat),]-calDat
	
	# check for zero var vars in resid ####
	for(i in 1:ncol(resDat)){
		if(var(resDat[,i],na.rm=T)==0){
			cat(sprintf('  Excluding %s for zero variance in resid\n',
					colnames(resDat)[i]))
			doneExcluding <- F
			exclude.idc <- c(exclude.idc,
											 which(varsForExport.cleanNames.orig==colnames(resDat)[i]))
		}
	}
	
	# check for perfect corr in resid ####
	if(doneExcluding){
		resDat.cor <- cor(resDat,use='complete.obs')
		perfCors <- which((resDat.cor-diag(1,ncol(resDat)))==1,arr.ind=T)
		if(nrow(perfCors)>0){
			for(i in 1:nrow(perfCors)){
				if(!(colnames(resDat)[perfCors[i,2]] %in% varsForExport.cleanNames.orig[exclude.idc])){
					cat(sprintf('  Excluding %s for perfect cor. with %s\n',
											colnames(resDat)[perfCors[i,1]],colnames(resDat)[perfCors[i,2]]))
					exclude.idc <- c(exclude.idc,
													 which(varsForExport.cleanNames.orig==colnames(resDat)[perfCors[i,1]]))
					doneExcluding <- F
				}
			}
		}
	}
	
	# check for further linear combos ####
	if(removeLinearCombinations){
		if(doneExcluding){
			linearCombos <- caret::findLinearCombos(resDat[complete.cases(resDat),])
			if(length(linearCombos$remove)>0){
				new.exclude.idc <- idxOfVarName(colnames(calDat)[linearCombos$remove],
																				varsForExport.cleanNames.orig)
				cat('  Excluding because they are linear combinations of other vars:\n    ')
				cat(paste(varsForExport.cleanNames.orig[new.exclude.idc],collapse='\n    '))
				cat('\n\n')
				exclude.idc <- unique(c(exclude.idc,new.exclude.idc))
				doneExcluding <- F
			}
		}
	}
}
calDat <- calDat.orig[,-exclude.idc]
cat('...done\n')
cat(sprintf('After exclusion we are left with %i out of %i variables\n',ncol(resDat),length(varsForExport.cleanNames.orig)))

# impute missing vars ####
if(imputeMissingVars||extrapolateMissingVarMethod!='n'){
	calDat.beforeImputeExtr <- calDat.orig[,-exclude.idc]
	ncols <- ncol(calDat)
	require(imputeTS)
	for(i in 1:ncol(calDat)){
		## setup
		x <- calDat.beforeImputeExtr[[i]]
		validRange <- funValidRange(x)
		## impute 
		if(imputeMissingVars){
			x[validRange[1]:validRange[2]] <- na_interpolation(x[validRange[1]:validRange[2]])
		}
		## extrapolate
		if(extrapolateMissingVarMethod=='f'){
			if(validRange[1]>1){
				x[1:(validRange[1]-1)] <- x[validRange[1]]
			}
			if(validRange[2]<length(x)){
				x[(validRange[2]+1):length(x)] <- x[validRange[2]]
			}
		}
		## write back
		calDat[[i]] <- x
	}
	calDat.impExtrValue <- calDat.beforeImputeExtr
	calDat.impExtrValue <- is.na(calDat.beforeImputeExtr)&!is.na(calDat)
}

# data Plots ####
completeIdc <- which(complete.cases(calDat))
incompleteIdc <- which(!complete.cases(calDat))
incompleteYears <- as.numeric(rownames(calDat)[incompleteIdc])
firstCompleteIdx <- completeIdc[1]
lastCompleteIdx <- completeIdc[length(completeIdc)]
years <- as.numeric(rownames(calDat))
if(plotWhileRunning){
	sqrtNcols <- sqrt(ncols)
	plotCols <- round(sqrtNcols)
	plotRows <- ceiling(sqrtNcols)
	rm(sqrtNcols)
	# par(mfrow=c(plotRows,plotCols))
	heights <- c(rep(2,plotRows),0.5)
	heights <- heights/sum(heights)
	layout(
		rbind(
			matrix((plotCols+1):(plotCols+plotRows*plotCols),
						 byrow=T,nrow=plotRows),
			t(1:plotCols)),
		widths = 1,
		heights = heights)
	par(mar=c(0,1,1,0.1))
	for(i in 1:plotCols){
		oldMar <- par('mar')
		par(mar=c(0,oldMar[2],0,oldMar[4]))
		plot(0,0,type='n',axes=F,
				 xlim=as.numeric(range(rownames(calDat))),
				 ylim=c(0,1))
		# abline(h=0)
		axis(1,pos=1)
		if(i == 1){
			mtext('Gray areas do not have complete cases',
						3,line = -4,adj=0)
			mtext('Red lines highlight the vars which most limit the complete cases window',
						3,line = -6,adj=0)
			mtext(sprintf('Complete cases: %i',nrow(calDat)-length(incompleteIdc)),
						3,line = -8,adj=0)
		}
		par(mar=oldMar)
	}
	for(i in 1:ncol(calDat)){
		plot(rownames(calDat),calDat[[i]],type='n',
				 xaxt='n',yaxt='n')
		for(year in incompleteYears){
			year <- as.numeric(year)
			rect(year-0.5,par('usr')[3],year+0.5,par('usr')[4],
					 density=-1,col='gray',border=NA)
		}
		box()
		# add label specifying var indext to top left
		text(years[1],max(calDat[[i]],na.rm=T),i,adj=c(0,1))
		# highlight if this var has its last obs on the boundary of incomplete obs
		validRange <- funValidRange(calDat[[i]])
		if(firstCompleteIdx!=1 && validRange[1]==firstCompleteIdx){
			abline(v=years[validRange[1]-1]+0.5,col='red',lwd=3)
		}
		if(lastCompleteIdx!=nrow(calDat) && validRange[2]==lastCompleteIdx){
			abline(v=years[validRange[2]+1]-0.5,col='red',lwd=3)
		}
		# add all points
		points(rownames(calDat),calDat[[i]])
		# add imputed points in red
		if(imputeMissingVars||extrapolateMissingVarMethod!='n'){
			points(rownames(calDat)[calDat.impExtrValue[,i]],
						 calDat[[i]][calDat.impExtrValue[,i]],
						 col='red',pch=20)
		}
		# add axis ticks
		axis(1,labels=F,tcl=0.3)
		yAxVals <- axis(2,labels=F,tcl=0.3)
		axTextCex <- 1
		# min yax label
		adjVal <- c(0.5,0)
		yPosVal <- min(yAxVals)
		if(ydev2in(yPosVal-par('usr')[3]) < 
			 strwidth(as.character(min(yAxVals)),
			 				 # font=par('font'),
			 				 family = 'Times New Roman',
			 				 units='inch',
			 				 cex=axTextCex)){
			adjVal <- c(0,0)
			yPosVal <- par('usr')[3]
		}
		text(par('usr')[1]-0.01*diff(par('usr')[1:2]),
				 yPosVal,
				 min(yAxVals),
				 xpd=T,cex=axTextCex,srt=90,adj=adjVal)
		# max yax label
		adjVal <- c(0.5,0)
		yPosVal <- max(yAxVals)
		if(ydev2in(yPosVal-par('usr')[4]) < 
			 strwidth(as.character(max(yAxVals)),
			 				 font=par('font'),
			 				 units='inch',
			 				 cex=axTextCex)){
			#0.1*(par('usr')[4]-par('usr')[3])){
			adjVal <- c(1,0)
			yPosVal <- par('usr')[4]
		}
		text(par('usr')[1]-0.01*diff(par('usr')[1:2]),
				 yPosVal,
				 max(yAxVals),
				 xpd=T,cex=axTextCex,srt=90,adj=adjVal)
		# highlight the labeld points
		points(rep(par('usr')[1],2),range(yAxVals),col=1,xpd=F,pch=18,cex=2.1)
		#
		text(mean(par('usr')[1:2]),par('usr')[4]+diff(par('usr')[3:4])*0.01,
				 colnames(calDat)[i],
				 cex=0.8,xpd=T,adj=c(0.5,0))
	}
}


# covar ####
cat('Determining the distribution of the residuals in the default case...\n')
resDat.cv <- NA
# this is not the correct way to do it. May not yield a positive definit matrix.
# but if it is valid uses more of the observations, so try it.
# See the following for why this is bad.
# https://www.r-bloggers.com/2015/06/pairwise-complete-correlation-considered-dangerous/
cat('  trying pariwise complete obs...')
try({resDat.cv <- cov(resDat,use='pairwise.complete.obs')})
if(sum(is.na(resDat.cv))>0||is.negative.definite(resDat.cv)){
	cat('fail\n')
	resDat.cv <- NA
} else {
	cat('success\n')
}
if(sum(is.na(resDat.cv))>0){
	cat(sprintf('  trying for complete obs, using %i obs...',sum(complete.cases(resDat))))
	tryCatch({resDat.cv <- cov(resDat,use='complete.obs')},
					 error=function(e){resDat.cv <- NA})
	if(sum(is.na(resDat.cv))>0){
		cat('fail\n')
		resDat.cv <- NA
	} else {
		cat('success\n')
	}
}
# if that there are no comple obs, give this a try but make it be quiet.
if(sum(is.na(resDat.cv))>0){
	cat('fail\n  trying FitGMMM...')
	sinkFile <- file('FitGMMM_babble','w')
	sink(sinkFile,type=c('message'))
	sink('FitGMM_babble2')
	try({
		resDat.mvn <- FitGMM(as.matrix(resDat),
												 init_means = list(a=colMeans(resDat,na.rm=T)),
												 init_covs = list(a=diag(1,nrow = ncol(resDat), ncol=ncol(resDat))),
												 #fix_means = T,
												 maxit = 1e5,
												 eps = 1e-6)
		resDat.mu <- resDat.mvn@Mean
		resDat.sigma <- resDat.mvn@Covariance
	})
	sink(NULL,type='message')
	sink()
	close(sinkFile)
	if(sum(is.na(resDat.cv))>0){
		cat('fail\n')
		resDat.cv <- NA
	} else {
		cat('success\n')
	}
}
if(sum(is.na(resDat.cv))>0){
	stop('The covariance matrix could not be calculated\n')
}
if(is.negative.definite(resDat.cv)){
	stop('The estimated covariance matrix is negative definite.\n')
}
if(is.singular.matrix(resDat.cv)){
	cat('Estimated covarianve matrix is singular.
			We have linear combinations in calDat but can continue.\n')
} else {
	cat('Estimated covariance matrix is fine\n')
}

# likelihood ####
precisBit <- 2048
resDat <- mpfr(as.matrix(resDat),precisBit)
resDat.cv <- mpfr(as.matrix(resDat.cv),precisBit)
funLikelihood(resDat,resDat.cv)

