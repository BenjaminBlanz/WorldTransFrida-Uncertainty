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
} else {
	stop('Missing covariance matrix file. Run runInitialiseData.R first.\n')
}
if(file.exists(file.path(location.output,'calDat.RDS'))){
	calDat.lst <- readRDS(file.path(location.output,'calDat.RDS'))
	calDat <- calDat.lst$calDat
	calDat.impExtrValue <- calDat.lst$calDat.impExtrValue
} else {
	stop('Missing calDat file. Run runInitialiseData.R first.\n')
}

# specify sampling parameters ####
# reads frida_info.csv and outputs the SampleParms
# also removes parms we will not sample
# and complains about invalid lines in frida_info.csv
if(file.exists(file.path(location.output,'sampleParms.RDS'))){
	cat('Reading sampling parameters...')
	sampleParms <- readRDS(file.path(location.output,'sampleParms.RDS'))
	cat('done\n')
} else {
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
	# cat('invalid parm specs min val max\n')
	# for(i in invalidLines){
	# 	cat(sprintf('%5i %-100s %10.f %10.f %10.f\n',
	# 							i,
	# 							sampleParms$Variable[i],
	# 							sampleParms$Min[i],
	# 							sampleParms$Value[i],
	# 							sampleParms$Max[i]))
	# }
	suppressWarnings(file.remove('frida_info_errorCases.csv'))
	if(length(invalidLines)>0){
		cat('invalid lines detected, see frida_info_errorCases.csv...')
	}
	write.csv(sampleParms[invalidLines,],'frida_info_errorCases.csv')
	sampleParms <- sampleParms[-invalidLines,]
	cat('done\n')
}
sampleParms.orig <- sampleParms
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
resSigmaVect <- as.vector(resSigma[!lower.tri(resSigma)])
names(resSigmaVect) <- as.vector(resSigma.names[!lower.tri(resSigma)])
jParVect <- c(parVect,resSigmaVect)

# start cluster ####
source('clusterHelp.R')

# MLE and Sensi Loop ####
parscale <- rep(NA,length(jParVect))
names(parscale) <- names(jParVect)
if(file.exists(file.path(location.output,'parscale.RDS'))){
	parscale.old <- readRDS('parscale.RDS')
	matches <- which(names(parscale) %in% names(parscale.old))
	parscale[matches] <- parscale.old[matches]
}
newMaxFound <- T
while(newMaxFound){
	# MLE ####
	cat('running fit procedure...')
	sv <- jParVect
	oldVal <- 0
	newVal <- 1
	# Optimisation of parameters (min neg log likelihood) is performed including
	# the covariance properties. The evaluation of likelihood of each of the parameters
	# for the uncertainty representation is performed with covariance matrix fixed to the
	# MLE.
	
	#determine parscale
	cat('Determining parscales...\n')
	orderOfMagNegLLErrorFun <- function(delta,par.i){
		jParVect.i <- jParVect
		jParVect.i[par.i] <- jParVect[par.i] + delta
		return(abs(baseNegLL-jnegLLikelihood.f(jParVect.i))-1)
	}
	baseNegLL <- jnegLLikelihood.f(jParVect)
	ordersOfMagGuesses <- c((floor(log10(abs(sampleParms$Max-sampleParms$Min)))-2),
													rep(-20,length(resSigmaVect)))
	ordersOfMagGuesses.orig <- ordersOfMagGuesses
	ordersOfMagLimits <- c(min(ordersOfMagGuesses)-2,max(ordersOfMagGuesses)+4)
	ordersOfMag <- seq(ordersOfMagLimits[1],ordersOfMagLimits[2])
	responseTolerance <- 0.01
	funFindParScale <- function(par.i,niter=100,guessOrderOfMag=NULL){
		if(is.null(guessOrderOfMag)){
			guessOrderOfMag <- min(ordersOfMag)
		}
		if(is.vector(guessOrderOfMag)){
			guessOrderOfMag <- guessOrderOfMag[par.i]
		}
		cat(sprintf('%4i %-50s ... magscale:     ',
								par.i,substr(names(jParVect)[par.i],1,50)))
		ordersOfMagDeltRes <- c()
		ordersOfMagNegLLResp <- c()
		for(ord.i in max(which(ordersOfMag<=guessOrderOfMag)):length(ordersOfMag)){
			cat(sprintf('\b\b\b\b\b\b\b\b\b\b%+10.1e',10^ordersOfMag[ord.i]))
			ordersOfMagDeltRes[ord.i] <- secant(orderOfMagNegLLErrorFun,x0=10^ordersOfMag[ord.i],
																					x1=10^ordersOfMag[ord.i]*1.1,
																					niter=niter,
																					doWarn=F,tol=1e-2,par.i=par.i)
			ordersOfMagNegLLResp[ord.i] <- orderOfMagNegLLErrorFun(ordersOfMagDeltRes[ord.i],par.i)
			if(!is.nan(ordersOfMagNegLLResp[ord.i])&&
				 abs(ordersOfMagNegLLResp[ord.i])<responseTolerance){
				break
			}
		}
		if(length(ordersOfMagDeltRes)==0){
			cat(sprintf('\r%4i %-50s ... %+e                     \n',
									par.i,substr(names(jParVect)[par.i],1,50),NA))
			return(NA)
		} else {
			retScale <- ordersOfMagDeltRes[which.min(ordersOfMagNegLLResp)]
			cat(sprintf('\r%4i %-50s ... %+e                     \n',
									par.i,substr(names(jParVect)[par.i],1,50),retScale))
			return(ordersOfMagDeltRes[which.min(ordersOfMagNegLLResp)])
		}
	}
	iterations <- 0
	parallelParscale <- T
	while(iterations < 2 && sum(is.na(parscale)|is.infinite(parscale))>0){
		if(parallelParscale){
			clusterExport(cl,list('baseNegLL',
														'ordersOfMagLimits','ordersOfMag','responseTolerance',
														'orderOfMagNegLLErrorFun','funFindParScale',
														'jnegLLikelihood.f','ordersOfMagGuesses',
														'calDat','resSigma',
														'jParVect'))
			gobble <- clusterEvalQ(cl,source(file.path(baseWD,'funParmSpace.R')))
			parsToDet <- which(is.na(parscale)|is.infinite(parscale))
			parParscaleOutput <- parLapplyLB(cl,parsToDet,funFindParScale,
																			 chunk.size = 1,
																			 guessOrderOfMag=ordersOfMagGuesses)
			parscale[parsToDet] <- unlist(parParscaleOutput)
			names(parscale) <- names(jParVect)
		} else {
			for(par.i in 1:length(jParVect)){
				if(is.na(parscale[par.i])){
					parscale[par.i] <- funFindParScale(par.i,guessOrderOfMag=ordersOfMagGuesses)
				}
			}
		}
		# try those that did not succeed with the guess again with the full range
		ordersOfMagGuesses <- ordersOfMagLimits[1]
		iterations <- iterations+1
	}
	cat('done\n')
	
	# check for bad behaviour in parscale ####
	
	
	problemCases <- which(is.infinite(parscale)|is.na(parscale))
	parscale[problemCases] <- ordersOfMagGuesses.orig[problemCases]
	
	# save parscale ####
	cat('saving ParScale...')
	saveRDS(parscale,file.path(location.output,'parscale.RDS'))
	
	stop()
	
	while(abs(oldVal-newVal)>1e-12){
		oldVal <- newVal
		# sv <- sv * 1.1
		optRes <- optimx(sv,jnegLLikelihood.f,method=c('Nelder-Mead'),
										 control=list(all.methods=F,
										 						 parscale = parscale,
										 						 # fnscale = newVal,
										 						 dowarn=F,
										 						 kkt=F,
										 						 maxit = 2e4,
										 						 reltol = 1e-15))
		newVal <- optRes$value
		rownames(optRes) <- NULL
		# print(optRes[1,c(1:12,15)])
		if(fitType=='sTime'){
			sv <- c(unlist(optRes[1,1:10]))
		} else if(fitType=='sTime2'||fitType=='sTime3'){
			sv <- c(unlist(optRes[1,1:8]))
		} else if (fitType=='nsplRa0b0'||fitType=='nsplBH'){
			sv <- c(unlist(optRes[1,1:3]))
		} else if (fitType=='nsplRa0b0MR1'){
			sv <- c(unlist(optRes[1,1:4]))
		} else if (fitType=='nsplTVRa2b2MR1'){
			sv <- c(unlist(optRes[1,1:8]))
		}
	}
	cat('optRes:\n')
	print(sv)
	jParVect.names <- names(jParVect)
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
