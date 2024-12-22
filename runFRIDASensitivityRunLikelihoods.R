#
#	provides a function to determine the likelihood of a given set of frida parameters
#
#
# 2024 Benjamin Blanz
# 

require(caret) # to find linear combinations and remove them in the calib dat
require(matrixcalc) # to test positive definitnes of cov matrix

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
# exclude all that have fewer than three data points as we
# can not determine variance in those cases.
numNA <- c()
for(i in 1:ncol(calDat.orig)){
	numNA[i] <- sum(is.na(calDat.orig[,i]))
}
exclude.idc <- which((nrow(calDat.orig)-numNA) < minObsForLike)

cat(sprintf('Excluding for less than minObsForLike (%i) observations:\n',minObsForLike))
cat(paste(colnames(calDat.orig)[exclude.idc],collapse='\n'))
cat('\n\n')


# find linear combinations ####
calDat <- calDat.orig[,-exclude.idc]
linearCombos <- caret::findLinearCombos(calDat[complete.cases(calDat),])
if(length(linearCombos$remove)>0){
	new.exclude.idc <- idxOfVarName(colnames(calDat)[linearCombos$remove],
																	varsForExport.cleanNames.orig)
	cat('Excluding because they are linear combinations of other vars:\n')
	cat(paste(varsForExport.cleanNames.orig[new.exclude.idc],collapse='\n'))
	cat('\n\n')
	exclude.idc <- unique(c(exclude.idc,new.exclude.idc))
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
	
	# try to fill in missing data in calDat ####
	# require(imputeTS)
	# calDat.infilled <- na_interpolation(calDat)
	
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
			cat(sprintf('Excluding %s for zero variance in resid\n',
					colnames(resDat)[i]))
			doneExcluding <- F
			exclude.idc <- c(exclude.idc,
											 which(varsForExport.cleanNames.orig==colnames(resDat)[i]))
		}
	}
	
	# check for further linear combos ####
	if(doneExcluding){
		linearCombos <- caret::findLinearCombos(resDat[complete.cases(resDat),])
		if(length(linearCombos$remove)>0){
			new.exclude.idc <- idxOfVarName(colnames(calDat)[linearCombos$remove],
																			varsForExport.cleanNames.orig)
			cat('Excluding because they are linear combinations of other vars:\n')
			cat(paste(varsForExport.cleanNames.orig[new.exclude.idc],collapse='\n'))
			cat('\n\n')
			exclude.idc <- unique(c(exclude.idc,new.exclude.idc))
			doneExcluding <- F
		}
	}
}


## covar ####
resDat.cv <- cov(resDat,use='complete.obs')
resDat.cor <- cor(resDat,use='complete.obs')

if(is.singular.matrix(resDat.cv)){
	cat('Estimated covarianve matrix is singular. We have uncessary variables in calDat\n')
}

if(!is.positive.definite(resDat.cv)){
	stop('The estimated covariance matrix is not positive definite.\n')
}




if(F){

resDat.mvn <- FitGMM(as.matrix(resDat),
										 init_means = list(a=colMeans(resDat,na.rm=T)),
										 init_covs = list(a=diag(1,nrow = ncol(resDat), ncol=ncol(resDat))),
										 # fix_means = T,
										 maxit = 1e5,
										 eps = 1e-10)
resDat.mu <- resDat.mvn@Mean
resDat.sigma <- resDat.mvn@Covariance



# calculating the covariance matrix.
# this is not the correct way to do it, may not yield a positive definit matrix.
resDat.cv <- cov(resDat,use='pairwise.complete.obs')





# Joint Model Likelihood Function ####
# 
# The model likelihood is likelihood that the residuals of the model are a 
# multivariate normal distribution. This is the standard assumption in most 
# fit procedures, including OLS.
# In addition to the model parameters we need the residual variance and covariance
# parameters. In principle these are also fit parameters. However, it is reasonable
# to treat them as nuisance parameters and use the values from the calibration result
# for these parameters for all likelihood evaluations.
# 

# vornoi in higher dimensions 
#https://cran.r-project.org/web/packages/geometry/index.html

# 
# Likelihood is that the 
# residuals of the model are a multivariate normal distribution with model parameters
# and residual variances and covariance specified from jParVect.
# Equations taken from: https://mathworld.wolfram.com/BivariateNormalDistribution.html
# 
# jParVect contains coefficeints of the mrate model, 
# followed by the sd of the mortality residuals,
# followed by the coefficients of the growth model,
# followed by the sd of the growth residuals
# followed by the covariance of the two residuals.
#          Idx: 1  2  3  4      5  6  7  8  9      10      
# jParVect == c(m1,m2,m3,msigma,g1,g2,g3,g4,gsigma,jcov)
jnegLLikelihood.f <- function(jParVect){
	mort.resid <- mort.resid.f(jParVect[1:3])
	growth.resid <- growth.resid.f(jParVect[5:8])
	x <- matrix(c(mort.resid,growth.resid),ncol=2,byrow = F)
	cov.mat <- matrix(c(jParVect[4]^2,jParVect[10],jParVect[10],jParVect[9]^2),nrow=2)
	return(-sum(dmvnorm(x,c(0,0),cov.mat,log=T)))
}

#note that the parVect here does not have the sigmas or jcov.
jnegLLikelihoodFixedCovMat.f <- function(jParVectNoSigma){
	mort.resid <- mort.resid.f(jParVectNoSigma[1:3])
	growth.resid <- growth.resid.f(jParVectNoSigma[4:7])
	x <- matrix(c(mort.resid,growth.resid),ncol=2,byrow = F)
	return(-sum(dmvnorm(x,c(0,0),cov.mat,log=T)))
}

# Growth Model Likelihood Function ####
# ignoring 
gnegLLikelihood.f <- function(gParVect){
	growth.resid <- growth.resid.f(gParVect)
	return(-sum(dnorm(growth.resid,0,sd(growth.resid))))
}

}
