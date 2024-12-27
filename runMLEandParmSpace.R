source('initialise.R')

# config ####
cat('Config...')
source('config.R')
if(file.exists(file.path(location.output,'sigma.RDS'))){
	resSigma <- readRDS(file.path(location.output,'sigma.RDS'))
	resSigma.names <- array(paste('s',
																rep(1:nrow(resSigma),ncol(resSigma)),
																rep(1:nrow(resSigma),each=nrow(resSigma)),
																sep='_'),
													dim=dim(resSigma))
	if(treatVarsAsIndep){
		# get the diagonal elements
		resSigma.var <- diag(resSigma)
		# make a diagonal matrix with those elements
		resSigma <- diag(resSigma.var)
	}
} else {
	stop('Missing covariance matrix file. Run runInitialiseData.R first.\n')
}
if(file.exists(file.path(location.output,'calDat.RDS'))){
	calDat.lst <- readRDS(file.path(location.output,'calDat.RDS'))
	calDat <- calDat.lst$calDat
	calDat.impExtrValue <- calDat.lst$calDat.impExtrValue
	calDat.orig <- calDat.orig
	calDat.withAllVars <- calDat.lst$calDat.withAllVars
} else {
	stop('Missing calDat file. Run runInitialiseData.R first.\n')
}

# specify sampling parameters ####
# reads frida_info.csv and outputs the SampleParms
# also removes parms we will not sample
# and complains about invalid lines in frida_info.csv
sampleParms <- prepareSampleParms()
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
	
	#determine parscale
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
		sampleParms <- prepareSampleParms(excludeNames = excludeParmNames)
		parVect <- sampleParms$Value
		names(parVect) <- sampleParms$Variable 
		jParVect <- c(parVect,resSigmaVect)
	}
	
	# save parscale ####
	cat('saving ParScale...')
	saveRDS(parscale,file.path(location.output,'parscale.RDS'))
	sampleParms$parscale <- parscale.parvect
	write.csv(sampleParms,file.path(location.output,'sampleParmsParscale.csv'))
	
	
	# optimize parameters ####
	sv <- jParVect
	optimOutput <- array(NA,dim=c(1,length(jParVect)+8))
	colnames(optimOutput) <- c(names(jParVect),
														 c('value','fevals','gevals','niter','convcode','kkt1','kkt2','xtime'))
	optimOutput <- as.data.frame(optimOutput)
	optimOutput[1,] <- c(jParVect,baseNegLL,rep(NA,7))
	rownames(optimOutput) <- 'sv'
	oldVal <- 0
	newVal <- 1
	iteration <- 0
	while(abs(oldVal-newVal)>1e-12){
		iteration <- iteration+1
		cat(sprintf('Running likelihood maximization (min neg log like) iteration %i...',
								iteration))
		oldVal <- newVal
		# sv <- sv * 1.1
		optRes <- optimx(sv,jnegLLikelihood.f,method=c('bobyqa'),
										 lower=c(sampleParms$Min,resSigmaVect-abs(parscale.resSigmaVect)*1e5),
										 upper=c(sampleParms$Max,resSigmaVect+abs(parscale.resSigmaVect)*1e5),
										 control=list(all.methods=F,
										 						 parscale = parscale,
										 						 # fnscale = newVal,
										 						 dowarn=F,
										 						 # trace=9,
										 						 kkt=F,
										 						 maxit = 2e4,
										 						 reltol = 1e-15))
		optimOutput[(nrow(optimOutput)+1):nrow(optRes),] <- optRes
		write.csv(optimOutput,file.path(location.output,'optRes.csv'))
		newVal <- optRes$value
		rownames(optRes) <- NULL
		sv <- unlist(optRes[1:length(jParVect)])
		cat(sptintf('%f\n',optRes$value[1]))
	}
	jParVect <- sv
	cat('completed optimization\n')
	
	# something else ####
	cat('optRes:\n')
	print(sv)
	jParVect.names <- names(jParVect)
	jParVect <- sv
	if(fitType=='sTime'){
		jParVect <- unlist(optRes[1:10])
		names(jParVect) <- jParVect.names
		mParVect <- jParVect[1:3]
		gParVect <- jParVect[5:8]
		cov.mat <- matrix(c(jParVect[4]^2,jParVect[10],jParVect[10],jParVect[9]^2),nrow=2)
		parVect.parscale <- parscale[c(1:3,5:8)]
	} else if (fitType=='sTime2'||fitType=='sTime3'){
		jParVect <- unlist(optRes[1:8])
		names(jParVect) <- jParVect.names
		mParVect <- jParVect[1:3]
		gParVect <- jParVect[5:6]
		cov.mat <- matrix(c(jParVect[4]^2,jParVect[8],jParVect[8],jParVect[7]^2),nrow=2)
		parVect.parscale <- parscale[c(1:3,5:6)]
	} else if (fitType=='nsplRa0b0'||fitType=='nsplBH'){
		jParVect <- unlist(optRes[1:3])
		names(jParVect) <- jParVect.names
		mParVect <- c()
		gParVect <- jParVect[1:2]
		cov.mat <- jParVect[3]
		parVect.parscale <- parscale[1:3]
	}else if (fitType=='nsplRa0b0MR1'){
		jParVect <- unlist(optRes[1:4])
		names(jParVect) <- jParVect.names
		mParVect <- jParVect[1]
		gParVect <- jParVect[2:3]
		cov.mat <- jParVect[4]
		parVect.parscale <- parscale[1:3]
	} else if (fitType=='nsplTVRa2b2MR1'){
		jParVect <- unlist(optRes[1:8])
		names(jParVect) <- jParVect.names
		mParVect <- jParVect[1]
		gParVect <- jParVect[2:7]
		cov.mat <- jParVect[8]
		parVect.parscale <- parscale[c(1,2:7)]
	} else {
		stop('Unkown fitType\n')
	}
	saveRDS(jParVect,file.path(storeCalcLocation2,paste0('jParVect-',fitType,'-',yrsForFit.name,'.RDS')))
	parVect <- c(mParVect,gParVect)
	
	
	sgslDat.l$jParVect <- jParVect
	sgslDat.l$mParVect <- mParVect
	sgslDat.l$gParVect <- gParVect
	sgslDat.l$coefs.gm.mle <- parVect
	sgslDat.l$llike.gm.mle <- -jnegLLikelihood.f(jParVect)
	sgslDat.l$sigma.gm.mle <- cov.mat
	sgslDat.l$m.pred.mle <- mort.f(mParVect)
	sgslDat.l$m.resid.mle <- sgslDat.l$m - sgslDat.l$m.pred.mle
	sgslDat.l$g.pred.mle <- growth.f(gParVect)
	sgslDat.l$g.resid.mle <- sgslDat.l$g - sgslDat.l$g.pred.mle
	cat('done\n')
	
	#### coef range ####
	# manually specify borders
	if(fitType=='sTime'){
		distances  <- c(1.25, #m1
										0.00075, #m2
										0.3e-6, #m3
										0.2, #g1
										4e-6, #g2
										4e-9, #g3
										2e-12) #g4
	} else if(fitType=='sTime2'||fitType=='sTime3'){
		distances  <- c(1.25, #m1
										0.00075, #m2
										0.3e-6, #m3
										0.2, #g1
										4e-6) #g2
	} else if(fitType=='nsplBH'){
		distances  <- c(9.24, 
										2e5) 
	}else if(fitType=='nsplRa0b0'){
		distances  <- c(9.24, #ga1
										2.03e-5) #gb1
	} else if(fitType=='nsplRa0b0MR1'){
		distances  <- c(0.36, #m1
										9.24, #ga1
										2.03e-5) #gb1
	} else if(fitType=='nsplTVRa2b2MR1'){
		distances  <- c(0.36, #m1
										9.24, #ga1
										4.68e-3, #ga2
										2.43e-6, #ga3
										2.03e-5, #gb1
										1.02e-8, #gb2
										5.07e-12) #gb3
	}
	if(n.sample==1){
		max.coefs <- parVect+distances
		min.coefs <- parVect-distances
	} else {
		# minimize and maximize each parameter with others free, until density is 
		# equal to pdensEps
		cat('determining coef sample range...')
		
		# boundary value in log likelihood
		lpdensEps <- -negLLike(parVect) - log(likeCutoffRatio)
		max.coefs <- parVect
		min.coefs <- parVect
		#limit to 4 c. due to mem. reqs. and also len(coefs)
		if(speratelyDetermineCoefRangeForMortAndGrowth){
			for(i in 1:length(parVect)){
				idcToMod <- list()
				if(i <= length(mParVect)){
					idcToMod[[i]] <- 1:length(mParVect)
				} else {
					idcToMod[[i]] <- (length(mParVect)+1):length(parVect)
				}
			}
		} else {
			idcToMod <- 1:length(parVect)
		}
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
		clRoot <- makeForkCluster(min(length(parVect),detectCores())) 
		cat('\rdetermining coef sample range... min\n')
		min.coefs <- unlist(parLapply(clRoot,1:length(parVect),findDensValBorder,
																	parVect=parVect,lpdensEps=lpdensEps,
																	ceterisParibusPars=ceterisParibusPars,
																	tol=1e-10,max=F,idcToMod=idcToMod,
																	parscale=parVect.parscale))
		names(min.coefs) <- names(parVect)
		# fallback values in case borders could not be determined:
		min.coefs[(is.infinite(min.coefs)+(parVect==min.coefs))>=1] <- 
			parVect[(is.infinite(min.coefs)+(parVect==min.coefs))>=1]-
			distances[(is.infinite(min.coefs)+(parVect==min.coefs))>=1]
		print(min.coefs)
		cat('\rdetermining coef sample range... max\n')
		max.coefs <- unlist(parLapply(clRoot,1:length(parVect),findDensValBorder,
																	parVect=parVect,lpdensEps=lpdensEps,
																	ceterisParibusPars=ceterisParibusPars,
																	tol=1e-10,max=T,idcToMod=idcToMod,
																	parscale=parVect.parscale))
		names(max.coefs) <- names(parVect)
		# fallback values in case borders could not be determined:
		max.coefs[(is.infinite(max.coefs)+(max.coefs==parVect))>=1] <- 
			parVect[(is.infinite(max.coefs)+(max.coefs==parVect))>=1]+
			distances[(is.infinite(max.coefs)+(parVect==max.coefs))>=1]
		print(max.coefs)
		stopCluster(clRoot)
		cat('\rdetermining coef sample range...done   \n')
		
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
	}
	
	#### llike ####
	if(runMaxLikeForAllSamples){
		n.sample <- n.sample.runMaxLikeForAllSamples
	} else {
		n.sample <- n.sample.full
	}
	par.vals <- array(NA,dim=c(n.sample,n.par))
	like.arr <- array(NA, dim=c(n.sample^n.par,n.par+2))
	colnames(like.arr) <- c(names(parVect),'llike','prob')
	if(!exists('samplingType') || samplingType=='raster'){
		for(i in 1:n.par){
			# par.vals[,i] <- seq(min.coefs[i],max.coefs[i],length.out=n.sample)
			par.vals[,i] <- sort(c(parVect[i],seq(min.coefs[i],max.coefs[i],length.out=n.sample-1)))
			# like.arr[2:(n.sample+1),i] <- runif(n.sample,min.coefs[i],max.coefs[i])
			# like.arr[,i] <- rnorm(n.sample,parVect[i],mean((parVect[i]-min.coefs[i])^2,
			# 																								(parVect[i]-max.coefs[i])^2))
		}
	} else if (samplingType=='orthogonal-latin-hypercube'){
		
	}
	sgslDat.l$mleCoefsAdded <- TRUE
	for(i in 1:n.par){
		like.arr[,i] <- rep(rep(unname(par.vals[,i]),each=n.sample^(i-1)),
												n.sample^(n.par-i))
	}
	# cl <- makePSOCKcluster(detectCores()/2)
	# clusterEvalQ(cl,source('runEarlyWarningPaper2-FitExploration-JointModelLikelihood.R'))
	# clusterExport(cl,list('n.par','like.arr','data','negLLike','parlike.f','cov.mat'))
	cl <- makeForkCluster(floor(detectCores()),renice=10)
	
	#### run opt par for all samplees ####
	likeMaxVals.max <- NA
	if(runMaxLikeForAllSamples){
		cat('Optimising Coef for each sample (n.sample=',n.sample,')...')
		parMaxLikelyFun <- function(i){
			sv <- c(like.arr[i,1:3],
							jParVect[4],
							like.arr[i,4:7],
							jParVect[9:10])
			optRes <- NA
			try({
				optRes <- optimx(sv,jnegLLikelihood.f,method=c('Nelder-Mead'),
												 control=list(all.methods=F,
												 						 parscale = parscale,
												 						 # fnscale = newVal,
												 						 dowarn=F,
												 						 kkt=F,
												 						 maxit = 2e4,
												 						 reltol = 1e-20))$value
			})
			return(unlist(optRes))
		}
		likeMaxVals <- unlist(parLapply(cl,1:(dim(like.arr)[1]),parMaxLikelyFun))
		maxInd <- which.min(likeMaxVals)
		likeMaxVals.max <- -likeMaxVals[maxInd]
		cat('done\n')
	} else {
		cat('calculating coef llike...')	# for testing# for(i in 1:(n.sample+1)){
		# 	like.arr[i,n.par+1] <- -negLLike(like.arr[i,1:n.par])
		
		# Note that
		#          Idx: 1  2  3  4  5  6  7  8  9  10      
		# jParVect == c(m1,m2,m3,ms,g1,g2,g3,g4,gs,jcov)
		# parVect  ==  c(m1,m2,m3,g1,g2,g3,g4)
		parlike.f <- function(i){
			-negLLike(like.arr[i,1:n.par])
		}
		like.arr[,'llike'] <- unlist(parLapply(cl,1:(dim(like.arr)[1]),parlike.f))
		cat('done\n')
		# check if a highler likelihood than at jParVect was found
		maxInd <- which.max(like.arr[,'llike'])
	}
	
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
