#
#	provides a function to determine the likelihood of a given set of frida parameters
#
#
# 2024 Benjamin Blanz
# 

suppressPackageStartupMessages(require(caret,quietly = T)) # to find linear combinations and remove them in the calib dat
suppressPackageStartupMessages(require(matrixcalc,quietly = T)) # to test positive definitnes of cov matrix
suppressPackageStartupMessages(require(imputeTS)) # used for interpolating missing values
suppressPackageStartupMessages(require(Rmpfr)) # use to calculate the likelihood from loglikelihood

source('config.R')
source('funRunFRIDA.R')
source('funPlot.R')

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
calDat <- calDat.orig <- calDat.orig[,-1]
cat('done\n')

# impute missing vars ####
if(imputeMissingVars||extrapolateMissingVarMethod!='n'){
	ncols <- ncol(calDat)
	for(i in 1:ncol(calDat)){
		## setup
		x <- calDat.orig[[i]]
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
	calDat.impExtrValue <- calDat.orig
	calDat.impExtrValue <- is.na(calDat.orig)&!is.na(calDat)
}
calDat.afterImpute <- calDat

# determine Vars to Exclude ####
cat('Applying rules to exclude variables from likelihood analysis...\n')
# exclude all that have fewer than three data points as we
# can not determine variance in those cases.
numNA <- c()
for(i in 1:ncol(calDat.afterImpute)){
	numNA[i] <- sum(is.na(calDat.afterImpute[,i]))
}
exclude.idc <- which((nrow(calDat.afterImpute)-numNA) < minObsForLike)

cat(sprintf('  Excluding for less than minObsForLike (%i) observations:\n    ',minObsForLike))
cat(paste(colnames(calDat.afterImpute)[exclude.idc],collapse='\n    '))
cat('\n')
calDat <- calDat.afterImpute[,-exclude.idc]

# exclude vars ####
doneExcluding <- F
while(!doneExcluding){
	doneExcluding <- T
	# calc resid
	calDat <- calDat.afterImpute[,-exclude.idc]
	varsForExport.fridaNames <- varsForExport.fridaNames.orig[-exclude.idc]
	writeFRIDAExportSpec(varsForExport.fridaNames)
	defDat <- runFridaDefaultParms()
	if(sum(colnames(defDat)!=colnames(calDat))!=0){
		stop('Mismatch in the columns of calibration data and model result data')
	}
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
}

# check for perfect corr in resid ####
# calc resid
calDat <- calDat.afterImpute[,-exclude.idc]
varsForExport.fridaNames <- varsForExport.fridaNames.orig[-exclude.idc]
writeFRIDAExportSpec(varsForExport.fridaNames)
defDat <- runFridaDefaultParms()
if(sum(colnames(defDat)!=colnames(calDat))!=0){
	stop('Mismatch in the columns of calibration data and model result data')
}
resDat <- defDat[1:nrow(calDat),]-calDat
resDat.cor <- cor(resDat,use='complete.obs')
perfCors <- which((resDat.cor-diag(1,ncol(resDat)))==1,arr.ind=T)
if(nrow(perfCors)>0){
	for(i in 1:nrow(perfCors)){
		if(!(colnames(resDat)[perfCors[i,2]] %in% varsForExport.cleanNames.orig[exclude.idc])){
			cat(sprintf('  Excluding %s for perfect cor. with %s\n',
									colnames(resDat)[perfCors[i,1]],colnames(resDat)[perfCors[i,2]]))
			exclude.idc <- c(exclude.idc,
											 which(varsForExport.cleanNames.orig==colnames(resDat)[perfCors[i,1]]))
		}
	}
}

# find linear combinations ####
##  in calDat ####
# new resid
calDat <- calDat.afterImpute[,-exclude.idc]
varsForExport.fridaNames <- varsForExport.fridaNames.orig[-exclude.idc]
writeFRIDAExportSpec(varsForExport.fridaNames)
defDat <- runFridaDefaultParms()
if(sum(colnames(defDat)!=colnames(calDat))!=0){
	stop('Mismatch in the columns of calibration data and model result data')
}
resDat <- defDat[1:nrow(calDat),]-calDat
if(removeLinearCombinations){
	linearCombos <- caret::findLinearCombos(resDat[complete.cases(resDat),])
	if(length(linearCombos$remove)>0){
		new.exclude.idc <- idxOfVarName(colnames(calDat)[linearCombos$remove],
																		varsForExport.cleanNames.orig)
		cat('  Excluding because they are linear combinations of other vars:\n    ')
		cat(paste(varsForExport.cleanNames.orig[new.exclude.idc],collapse='\n    '))
		cat('\n')
		exclude.idc <- unique(c(exclude.idc,new.exclude.idc))
	}
}

## in calDat.cv ####
# new resid
calDat <- calDat.afterImpute[,-exclude.idc]
varsForExport.fridaNames <- varsForExport.fridaNames.orig[-exclude.idc]
writeFRIDAExportSpec(varsForExport.fridaNames)
defDat <- runFridaDefaultParms()
if(sum(colnames(defDat)!=colnames(calDat))!=0){
	stop('Mismatch in the columns of calibration data and model result data')
}
resDat <- defDat[1:nrow(calDat),]-calDat

tryCatch({resDat.cv <- cov(resDat,use='complete.obs')},
				 error=function(e){resDat.cv <- NA})
if(sum(is.na(resDat.cv))==0){
	if(is.singular.matrix(resDat.cv)){
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
}
calDat <- calDat.afterImpute[,-exclude.idc]
cat('...done\n')
cat(sprintf('After exclusion we are left with %i out of %i variables\n',ncol(resDat),length(varsForExport.cleanNames.orig)))

# remove excluded idc ###
calDat.orig <- calDat.orig[,-exclude.idc]
calDat <- calDat.afterImpute[,-exclude.idc]
calDat.impExtrValue <- calDat.impExtrValue[,-exclude.idc]
rm(calDat.afterImpute)

# clean def run ####
varsForExport.fridaNames <- varsForExport.fridaNames.orig[-exclude.idc]
writeFRIDAExportSpec(varsForExport.fridaNames)
defDat <- runFridaDefaultParms()
resDat <- defDat[1:nrow(calDat),]-calDat

# data Plots ####
if(plotWhileRunning){
	funPlotDat(calDat,calDat.impExtrValue,defDat)
	dev.print(pdf,
						file.path(location.output,'calDatPlot.pdf'))
	funPlotDat(resDat,calDat.impExtrValue)
	dev.print(pdf,
						file.path(location.output,'resDatPlot.pdf'))
}

# covar ####
if(sum(colnames(defDat)!=colnames(calDat))!=0){
	stop('Mismatch in the columns of calibration data and model result data')
}


cat('Determining the distribution of the residuals in the default case...\n')
resDat.cv <- NA
# this is not the correct way to do it. May not yield a positive definit matrix.
# but if it is valid uses more of the observations, so try it.
# See the following for why this is bad.
# https://www.r-bloggers.com/2015/06/pairwise-complete-correlation-considered-dangerous/
cat('  trying pariwise complete obs...')
try({resDat.cv <- cov(resDat,use='pairwise.complete.obs')})
if(sum(is.na(resDat.cv))>0||!is.positive.definite(resDat.cv)){
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
# if there are no comple obs, give this a try but make it be quiet.
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
	stop('The estimated covariance matrix is negative definite.\nLikelihood can not be calculated\n')
}
if(is.singular.matrix(resDat.cv)){
	cat('Estimated covarianve matrix is singular.\nMight be a problem.\n')
} else {
	cat('Estimated covariance matrix is fine\n')
}

# likelihood ####
defLike <- exp(mpfr(funLogLikelihood(resDat[complete.cases(resDat),],resDat.cv),64))
cat(paste0('Likelihood in the default case: ',formatMpfr(defLike)))
if(is.infinite(defLike)||is.na(defLike)){
	stop('Bad default log likelihood\n')
}

# save run prep ####
writeFRIDAExportSpec(varsForExport.fridaNames.orig[-exclude.idc])
saveRDS(resDat.cv,file.path(location.output,'sigma.RDS'))
saveRDS(list(calDat=calDat,calDat.impExtrValue=calDat.impExtrValue),file.path(location.output,'calDat.RDS'))

