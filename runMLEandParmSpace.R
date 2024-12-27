source('initialise.R')

# config ####
cat('Config...')
source('config.R')
# read covariance matrix used for baseNegLL
if(treatVarsAsIndep&&
	 file.exists(file.path(location.output,'sigma-indepParms.RDS'))){
	resSigma <- readRDS(file.path(location.output,'sigma-indepParms.RDS'))
} else if(file.exists(file.path(location.output,'sigma.RDS'))){
	resSigma <- readRDS(file.path(location.output,'sigma.RDS'))
	if(treatVarsAsIndep){
		# get the diagonal elements
		resSigma.var <- diag(resSigma)
		# make a diagonal matrix with those elements
		resSigma <- diag(resSigma.var)
	}
} else {
	stop('Missing covariance matrix file. Run runInitialiseData.R first.\n')
}
resSigma.names <- array(paste('s',
															rep(1:nrow(resSigma),ncol(resSigma)),
															rep(1:nrow(resSigma),each=nrow(resSigma)),
															sep='_'),
												dim=dim(resSigma))
# read calibration data
if(file.exists(file.path(location.output,'calDat.RDS'))){
	calDat.lst <- readRDS(file.path(location.output,'calDat.RDS'))
	calDat <- calDat.lst$calDat
	calDat.impExtrValue <- calDat.lst$calDat.impExtrValue
	calDat.orig <- calDat.lst$calDat.orig
	calDat.withAllVars <- calDat.lst$calDat.withAllVars
} else {
	stop('Missing calDat file. Run runInitialiseData.R first.\n')
}

# specify sampling parameters ####
# reads frida_info.csv and outputs the SampleParms
# also removes parms we will not sample
# and complains about invalid lines in frida_info.csv
sampleParms.orig <- sampleParms <- prepareSampleParms()
saveRDS(sampleParms,file.path(location.output,'sampleParms.RDS'))


# mle and like ####
# Optimisation of parameters (min neg log likelihood) is performed including
# the covariance properties. The evaluation of likelihood of each of the parameters
# for the uncertainty representation is performed with covariance matrix fixed to the
# MLE.

# starting value is the value column from frida_info.csv a known good run.
parVect <- sampleParms$Value
names(parVect) <- sampleParms$Variable
# jParVect contains all fit parameters, including covariance matrix
# parVect contains the sampled fit parameters
# Note that
# jParVect == c(parVect,covarVect)
if(treatVarsAsIndep){
	resSigmaVect <- diag(resSigma)
	names(resSigmaVect) <- diag(resSigma.names)
} else {
	resSigmaVect <- as.vector(resSigma[!lower.tri(resSigma)])
	names(resSigmaVect) <- as.vector(resSigma.names[!lower.tri(resSigma)])
}
jParVect <- c(parVect,resSigmaVect)

# start cluster ####
source('clusterHelp.R')

# MLE and Sensi Loop ####
parscale <- rep(NA,length(jParVect))
names(parscale) <- names(jParVect)
if(file.exists(file.path(location.output,'parscale.RDS'))){
	parscale.old <- readRDS(file.path(location.output,'parscale.RDS'))
	matches <- which(names(parscale) %in% names(parscale.old))
	parscale[matches] <- parscale.old[matches]
}
ordersOfMagGuesses.parvect <- funOrderOfMagnitude(sampleParms$Max-sampleParms$Min)
ordersOfMagGuesses.resSigmaVect <- funOrderOfMagnitude(resSigmaVect)-6
ordersOfMagGuesses <- c(ordersOfMagGuesses.parvect,ordersOfMagGuesses.resSigmaVect)

# used by the funFindParScale function
ordersOfMagLimits <- c(min(ordersOfMagGuesses)-2,max(ordersOfMagGuesses)+4)
ordersOfMag <- seq(ordersOfMagLimits[1],ordersOfMagLimits[2])
responseTolerance <- 0.01

newMaxFound <- T
# while(newMaxFound){
	# MLE ####
	cat('running fit procedure...')
	# Optimisation of parameters (min neg log likelihood) is performed including
	# the covariance properties. The evaluation of likelihood of each of the parameters
	# for the uncertainty representation is performed with covariance matrix fixed to the
	# MLE.
	
	# determine parscale ####
	cat('Determining parscales...\n')
	baseNegLL <- jnegLLikelihood.f(jParVect)
	
	iterations <- 0
	parallelParscale <- T
	useOrdersOfMagGuesses <- T
	while(iterations < 2 && sum(is.na(parscale)|is.infinite(parscale))>0){
		parsToDet <- which(is.na(parscale)|is.infinite(parscale))
		cat(sprintf('Determining the parscale of %i parameters. %i parameters with already known parscale.%s\n',
								length(parsToDet),length(parscale)-length(parsToDet),
								if(useOrdersOfMagGuesses){' Using guesses.'}else{' Not using guesses.'}))
		if(parallelParscale){
			clusterExport(cl,list('baseNegLL',
														'ordersOfMagLimits','ordersOfMag','responseTolerance',
														'orderOfMagNegLLErrorFun','funFindParScale',
														'jnegLLikelihood.f','ordersOfMagGuesses',
														'ordersOfMagLimits',
														'calDat','resSigma',
														'jParVect'))
			gobble <- clusterEvalQ(cl,source(file.path(baseWD,'funParmSpace.R')))
			parParscaleOutput <- parLapplyLB(cl,parsToDet,funFindParScale,
																			 useOrdersOfMagGuesses=useOrdersOfMagGuesses)
			parscale[parsToDet] <- unlist(parParscaleOutput)
			names(parscale) <- names(jParVect)
		} else {
			for(par.i in 1:length(jParVect)){
				if(is.na(parscale[par.i])|is.infinite(parscale[par.i])){
					parscale[par.i] <- funFindParScale(par.i,useOrdersOfMagGuesses=useOrdersOfMagGuesses)
				}
			}
		}
		# cat('saving ParScale...')
		# saveRDS(parscale,file.path(location.output,'parscale.RDS'))
		# cat('done\n')
		# try those that did not succeed with the guess again with the full range
		useOrdersOfMagGuesses <- F
		iterations <- iterations+1
	}
	parscale.parvect <- parscale[1:nrow(sampleParms)]
	parscale.resSigmaVect <- parscale[(nrow(sampleParms)+1):length(jParVect)]
	cat('done\n')
	
	# check for bad behaviour in parscale ####
	#TODO: deal with bad behaviour in parscale: Kick them out
	# only the entries in parVect can be excluded. The entries in resSigmaVect need to 
	# be delt with. E.g. by using the guess values. The maximum likelihood vars (diag
	# elements of the covmat can always be determined as the variance of those obs.
	problemCasesIdc <- which(is.infinite(parscale)|is.na(parscale))
	problemCasesIdc.parVect <- which(is.infinite(parscale.parvect)|is.na(parscale.parvect))
	problemCasesIdc.resSigmaVect <- which(is.infinite(parscale.resSigmaVect)|is.na(parscale.resSigmaVect))
	cat(sprintf('%i parscales could not be determined.',length(problemCasesIdc)))
	if(length(problemCasesIdc.resSigmaVect)>0){
		cat(sprintf('  %i in resSigmaVect, guesses will be used',
								length(problemCasesIdc.resSigmaVect)))
		parscale.resSigmaVect[problemCasesIdc.resSigmaVect] <- 
			10^ordersOfMagGuesses.resSigmaVect[problemCasesIdc.resSigmaVect]
	} else {
		cat('  No problem cases in resSigmaVect\n')
	}
	if(length(problemCasesIdc.parVect)>0){
		cat(sprintf('  %i in parVect, these parms will be dropped\n',length(problemCasesIdc.parVect)))
		parscale.parvect <- parscale.parvect[-problemCasesIdc.parVect]
		excludeParmNames <- sampleParms$Variable[problemCasesIdc.parVect]
		cat(paste(excludeParmNames,collapse='\n'))
		cat('\n')
		sampleParms <- prepareSampleParms(excludeNames = excludeParmNames)
		if(file.exists('parExclusionList.csv')&&file.size('parExclusionList.csv')>0){
			oldExclusionList <- read.csv('parExclusionList.csv')
			exclusionList <- unique(c(oldExclusionList$excludedName,excludeParmNames))
		} else {
			exclusionList <- data.frame(excludedName=excludeParmNames)
		}
		write.csv(exclusionList,'parExclusionList.csv')
		parVect <- sampleParms$Value
		names(parVect) <- sampleParms$Variable 
		jParVect <- c(parVect,resSigmaVect)
	}
	parscale.all <- parscale
	parscale <- c(parscale.parvect,parscale.resSigmaVect)
	
	# save parscale ####
	cat('saving ParScale...')
	saveRDS(parscale.all,file.path(location.output,'parscale.RDS'))
	sampleParms$parscale <- parscale.parvect
	write.csv(sampleParms,file.path(location.output,'sampleParmsParscale.csv'))
	
	
	# optimize parameters ####
	sv <- jParVect
	optimOutput <- array(NA,dim=c(1,length(jParVect)+9))
	colnames(optimOutput) <- c(names(jParVect),
														 c('value','fevals','gevals','niter','convcode',
														 	'kkt1','kkt2','xtime','check value'))
	optimOutput <- as.data.frame(optimOutput)
	optimOutput[1,] <- c(jParVect,baseNegLL,rep('',8))
	rownames(optimOutput) <- 'sv'
	oldVal <- 0
	newVal <- 1
	iteration <- 0
	all.methods <- T
	while(abs(oldVal-newVal)>1e-12){
		iteration <- iteration+1
		cat(sprintf('Running likelihood maximization (min neg log like) iteration %i...',
								iteration))
		oldVal <- newVal
		# sv <- sv * 1.1
		# specifying limits breaks the parscale info for bobyqa!
		optRes <- optimx(sv,jnegLLikelihood.f,method=c('bobyqa'),
										 control=list(all.methods=all.methods,
										 						 parscale = parscale,
										 						 # fnscale = newVal,
										 						 dowarn=F,
										 						 # trace=9,
										 						 kkt=F,
										 						 maxit = 2e4, # recommended: 10*length(jParVect)^2,
										 						 reltol = 1e-15))
		svNegLLike <-c ()
		for(opt.i in 1:nrow(optRes)){
			sv.i <- unlist(as.vector(optRes[opt.i,1:length(jParVect)]))
			svNegLLike[opt.i] <- jnegLLikelihood.f(sv.i)
		}
		maxMethod <- which.min(svNegLLike)
		sv <- unlist(as.vector(optRes[maxMethod,1:length(jParVect)]))
		cat(sprintf('%10f %10f\n',
								optRes$value[1],svNegLLike[maxMethod]))
		newOptimOutputRowNums <- (nrow(optimOutput)+1):((nrow(optimOutput))+nrow(optRes))
		optimOutput[newOptimOutputRowNums,] <- 
			base::cbind(optRes,svNegLLike)
		rownames(optimOutput)[newOptimOutputRowNums] <-
			paste(rep(iteration,nrow(optRes)),rownames(optRes))
		write.csv(optimOutput,file.path(location.output,'optRes.csv'))
		newVal <- optRes$value[maxMethod]
		all.methods <- F
	}
	if(sum(is.na(sv[1:length(jParVect)]))==0|optRes$value<baseNegLL){
		jParVect.names <- names(jParVect)
		jParVect <- sv[1:length(jParVect)]
		names(jParVect) <- jParVect.names
		parVect <- jParVect[1:length(parVect)]
		resSigmaVect <- jParVect[(length(parVect)+1):length(jParVect)]
		if(treatVarsAsIndep){
			resSigma <- diag(resSigmaVect)
		} else {
			resSigma <- array(NA,dim=rep(ncol(calDat),2))
			resSigma[!lower.tri(resSigma)]<- resSigmaVect
			resSigma[lower.tri(resSigma)] <- t(resSigma)[lower.tri(resSigma)]
		}
		saveRDS(jParVect,file.path(location.output,'jParVectAfterOptim.RDS'))
		saveRDS(resSigma,file.path(location.output,
															 paste0('sigma',
															 			 ifelse(treatVarsAsIndep,'-indepParms',''),
															 			 '.RDS')))
		cat('completed optimization\n')
	} else {
		stop('failed optimization\n')
	}
	
	# coef range ####
	# manually specify borders
	max.coefs.prior <- max.coefs <- sampleParms$Max
	min.coefs.prior <- min.coefs <- sampleParms$Min
	# minimize and maximize each parameter with others free, until density is 
	# equal to pdensEps
	cat('determining coef sample range...\n')
	# boundary value in log likelihood
	lpdensEps <- -negLLike(parVect) - log(likeCutoffRatio)
	#limit to 4 c. due to mem. reqs. and also len(coefs)
	idcToMod <- 1:length(parVect)
	## for testing
	# for(i in 1:length(parVect)){
	# 	cat('\n\nmax parm ',i,'\n')
	# 	findDensValBorder(i,parVect=parVect,lpdensEps=lpdensEps,
	# 										ceterisParibusPars=ceterisParibusPars,
	# 										tol=1e-10,max=T,trace=1,idcToMod=idcToMod,
	# 										maxiter=1e3,
	# 										parscale=parVect.parscale)
	# 	cat('\nmin parm ',i,'\n')
	# 	findDensValBorder(i,parVect=parVect,lpdensEps=lpdensEps,
	# 										ceterisParibusPars=ceterisParibusPars,
	# 										tol=1e-10,max=F,trace=1,idcToMod=idcToMod,
	# 										maxiter=1e3,
	# 										parscale=parVect.parscale)
	# }
	
	cat('  determining min par values...')
	min.coefs <- unlist(parLapply(cl,1:length(parVect),findDensValBorder,
																parVect=parVect,lpdensEps=lpdensEps,
																ceterisParibusPars=treatVarsAsIndep,
																tol=1e-10,max=F,idcToMod=idcToMod,
																parscale=parVect.parscale))
	names(min.coefs) <- names(parVect)
	# fallback values in case borders could not be determined:
	notDeterminedMinBorders <- which((is.infinite(min.coefs)+(parVect==min.coefs))>=1)
	min.coefs[notDeterminedMinBorders] <- min.coefs.prior[notDeterminedMinBorders]
	cat(sprintf('done %i failures\n',length(notDeterminedMinBorders)))
	
	cat('  determining max par values...')
	max.coefs <- unlist(parLapply(cl,1:length(parVect),findDensValBorder,
																parVect=parVect,lpdensEps=lpdensEps,
																ceterisParibusPars=treatVarsAsIndep,
																tol=1e-10,max=T,idcToMod=idcToMod,
																parscale=parVect.parscale))
	names(max.coefs) <- names(parVect)
	# fallback values in case borders could not be determined:
	notDeterminedMaxBorders <- which((is.infinite(max.coefs)+(max.coefs==parVect))>=1)
	max.coefs[notDeterminedMaxBorders] <- max.coefs.prior[notDeterminedMaxBorders]
	cat(sprintf('done %i failures\n',length(notDeterminedMinBorders)))
	
	cat('  saving par ranges...')
	sampleParms$oldMin <- sampleParms.orig$Min
	sampleParms$oldMax <- sampleParms.orig$Max
	sampleParms$Max <- max.coefs
	sampleParms$Min <- min.coefs
	write.csv(sampleParms,file.path(location.output,'sampleParmsParscale.csv'))
	cat('...done   \n')
	
	# look at the min max and opt coefs
	# if(FALSE){
	# mat.coefs <- matrix(c(min.coefs,parVect,max.coefs),nrow=3,byrow=T)
	# colnames(mat.coefs) <- names(parVect)
	# rownames(mat.coefs) <- c('min.coefs','parVect','max.coefs')
	# View(mat.coefs)
	# }
	# reselect coefs
	# max.coefs <- par.vals[6,]
	# min.coefs <- par.vals[4,]
	
	#### llike ####
	
	
	if(max(like.arr[maxInd,'llike'],likeMaxVals.max,na.rm = T) > -jnegLLikelihoodFixedCovMat.f(parVect)){
		if(fitType=='sTime'){
			jParVect <- c(like.arr[maxInd,1:3],
										jParVect[4],
										like.arr[maxInd,4:7],
										jParVect[9:10])
		} else if (fitType=='sTime2'||fitType=='sTime3'){
			jParVect <- c(like.arr[maxInd,1:3],
										jParVect[4],
										like.arr[maxInd,4:5],
										jParVect[7:8])
		} else if (fitType=='nsplRa0b0'||fitType=='nsplBH'){
			jParVect <- c(like.arr[maxInd,1:2],
										jParVect[3])
		} else if (fitType=='nsplRa0b0MR1'){
			jParVect <- c(like.arr[maxInd,1],
										like.arr[maxInd,2:3],
										jParVect[4])
		} else if (fitType=='nsplTVRa2b2MR1'){
			jParVect <- c(like.arr[maxInd,1],
										like.arr[maxInd,2:7],
										jParVect[8])
		}
		newMaxFound <- T
		cat('Found greater likelihood pars in sampling, rerunning fit procedure\n')
	} else if(runMaxLikeForAllSamples && n.sample.full!=n.sample.runMaxLikeForAllSamples){
		runMaxLikeForAllSamples <- F
		cat('No greater likelihood found in per sample opt, now rerunning with full samples\n')
	} else if(runMaxLikeForAllSamples){
		runMaxLikeForAllSamples <- F
		cat('No greater likelihood found in per sample opt, now rerunning for fit likelihood\n')
	} else {
		newMaxFound <- F	
	}
}
