source('initialise.R')

# config ####
cat('Config...')
source('config.R')
source('runInitialiseData.R')
continue <- readline(paste0('Output location created. Move any files to be used here\n',
														location.output,'\nHit ENTER when done.\n'))
source('setupTMPFS.R')
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
resSigma.names <- array(paste('s',
															rep(colnames(calDat),ncol(resSigma)),
															rep(colnames(calDat),each=nrow(resSigma)),
															sep='_X_'),
												dim=dim(resSigma))

# specify sampling parameters ####
# reads frida_info.csv and outputs the SampleParms
# also removes parms we will not sample
# and complains about invalid lines in frida_info.csv
integerParms <- read.csv(file.path(location.frida.info,name.frida_integer_parms))
excludedParmsForBeingIntegers <- integerParms$Variable
sampleParms.orig <- sampleParms <- prepareSampleParms(excludeNames=excludedParmsForBeingIntegers)

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
if(!redoAllCalc){
	names(parscale) <- names(jParVect)
	if(file.exists(file.path(location.output,'parscale.RDS'))){
		parscale.old <- readRDS(file.path(location.output,'parscale.RDS'))
		matches <- which(names(parscale) %in% names(parscale.old))
		parscale[matches] <- parscale.old[matches]
	}
}
ordersOfMagGuesses.parvect <- funOrderOfMagnitude(sampleParms$Max-sampleParms$Min)
ordersOfMagGuesses.resSigmaVect <- funOrderOfMagnitude(resSigmaVect)-6
ordersOfMagGuesses <- c(ordersOfMagGuesses.parvect,ordersOfMagGuesses.resSigmaVect)

# used by the funFindParScale function
ordersOfMagLimits <- c(min(ordersOfMagGuesses)-2,max(ordersOfMagGuesses)+4)
ordersOfMag <- seq(ordersOfMagLimits[1],ordersOfMagLimits[2])
responseTolerance <- 0.01

#
frida_info <- read.csv(file.path(location.frida.info,name.frida_info))

newMaxFound <- T
iterationNewMax <-0
while(newMaxFound){
	iterationNewMax <- iterationNewMax+1
	cat(sprintf('running everything iteration %i...\n',iterationNewMax))
	# Optimisation of parameters (min neg log likelihood) is performed including
	# the covariance properties. The evaluation of likelihood of each of the parameters
	# for the uncertainty representation is performed with covariance matrix fixed to the
	# MLE.
	
	baseNegLL <- jnegLLikelihood.f(jParVect)
	
	if(forceParBounds){
		#cat('Forced using frida_info bounds\n')
	} else if(!redoAllCalc&&file.exists(file.path(location.output,'sampleParmsParscaleRanged.RDS'))){
		cat('loading existing sampleParmsParscaleRanged\n')
		sampleParms <- readRDS(file.path(location.output,'sampleParmsParscaleRanged.RDS'))
		parVect <- sampleParms$Value
		names(parVect) <- sampleParms$Variable 
		jParVect <- c(parVect,resSigmaVect)
	} else {
		# determine parscale ####
		cat('Determining parscales...\n')
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
		
		## check for bad behaviour in parscale ####
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
			sampleParms <- prepareSampleParms(excludeNames = c(excludeParmNames,excludedParmsForBeingIntegers))
			if(file.exists(file.path(location.frida.info,name.frida_parameter_exclusion_list))&&file.size(file.path(location.frida.info,name.frida_parameter_exclusion_list))>0){
				oldExclusionList <- read.csv(file.path(location.frida.info,name.frida_parameter_exclusion_list))
				exclusionList <- data.frame(excludedName=unique(c(oldExclusionList$excludedName,excludeParmNames)))
			} else {
				exclusionList <- data.frame(excludedName=excludeParmNames)
			}
			write.csv(exclusionList,file.path(location.frida.info,name.frida_parameter_exclusion_list))
			parVect <- sampleParms$Value
			names(parVect) <- sampleParms$Variable 
			jParVect <- c(parVect,resSigmaVect)
		}
		parscale.all <- parscale
		parscale <- c(parscale.parvect,parscale.resSigmaVect)
		
		## save parscale ####
		cat('saving ParScale...')
		saveRDS(parscale.all,file.path(location.output,'parscale.RDS'))
		sampleParms$parscale <- parscale.parvect
		write.csv(sampleParms,file.path(location.output,'sampleParmsParscale.csv'))
		# MLE ####
		if(!skipParMLE){
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
			all.methods <- T # use all methods on the first iteration then use whichever was the best
			methods <- c('bobyqa')
			while(abs(oldVal-newVal)>1e-12&&iteration<1e3){
				iteration <- iteration+1
				cat(sprintf('Running likelihood maximization (min neg log like) iteration %i...',
										iteration))
				oldVal <- newVal
				# sv <- sv * 1.1
				# specifying limits breaks the parscale info for bobyqa!
				lower <- c(sampleParms$Min,resSigmaVect-abs(parscale.resSigmaVect)*100)
				names(lower) <- names(jParVect)
				which(lower==sv)
				upper <- c(sampleParms$Min,resSigmaVect+abs(parscale.resSigmaVect)*100)
				optRes <- optimx(sv,jnegLLikelihood.f,method=methods,
												 lower = lower,
												 upper = upper,
												 control=list(all.methods=all.methods,
												 						 parscale = 1/parscale,
												 						 # fnscale = newVal,
												 						 dowarn=F,
												 						 # trace=9,
												 						 kkt=F,
												 						 maxit = 10*length(jParVect)^2,
												 						 reltol = 1e-15))
				svNegLLike <-c ()
				for(opt.i in 1:nrow(optRes)){
					sv.i <- unlist(as.vector(optRes[opt.i,1:length(jParVect)]))
					svNegLLike[opt.i] <- jnegLLikelihood.f(sv.i)
				}
				maxMethod <- which.min(svNegLLike)
				methods <- rownames(optRes[which(!is.na(optRes[,1]))])
				sv <- unlist(as.vector(optRes[maxMethod,1:length(jParVect)]))
				cat(sprintf('%10f %10f\n',
										optRes$value[1],svNegLLike[maxMethod]))
				newOptimOutputRowNums <- (nrow(optimOutput)+1):((nrow(optimOutput))+nrow(optRes))
				optimOutput[newOptimOutputRowNums,] <- 
					base::cbind(optRes,svNegLLike)
				rownames(optimOutput)[newOptimOutputRowNums] <-
					paste(rep(iteration,nrow(optRes)),rownames(optRes))
				write.csv(optimOutput,file.path(location.output,'optRes.csv'),)
				saveRDS(optRes,file.path(location.output,'optRes.RDS'))
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
		}
	}
	# coef range ####
	## par bounds ####
	idcOfSampleParmsInFridaInfo <- c()
	for(p.i in 1:nrow(sampleParms)){
		idcOfSampleParmsInFridaInfo[p.i] <- which(frida_info$Variable==sampleParms$Variable[p.i])
	}
	parBounds <- frida_info[idcOfSampleParmsInFridaInfo,c('Min','Max')]
	lpdensEps <- -negLLike(parVect) - log(likeCutoffRatio)
	notDeterminedBorders <- array(TRUE,dim=c(length(parVect),2))
	colnames(notDeterminedBorders) <- c('Min','Max')
	border.coefs <- notDeterminedBorders
	if(forceParBounds){
		cat('Forcing coefs sample range to be equal tovalues frida_info\n')
		border.coefs <- sampleParms[,c('Min','Max')]
	} else if (!file.exists(file.path(location.output,'sampleParmsParscaleRanged.RDS'))) {
		# minimize and maximize each parameter with others free, until density is 
		# equal to pdensEps
		cat('determining coef sample range...\n')
		# boundary value in log likelihood
		idcToMod <- 1:length(parVect)
		if(ignoreParBounds){
			parBounds <- array(c(rep(-.Machine$double.xmax,length(parVect)),rep(.Machine$double.xmax,length(parVect))),
												 dim=c(length(parVect),2))
		}
		rownames(parBounds) <- sampleParms$Variable
		colnames(parBounds) <- c('Min','Max')
		## for testing
		# test.i <- 1
		# #min bound
		# for(test.i in 1:nrow(sampleParms)){
		# 	findDensValBorder(test.i,
		# 										parVect=parVect,lpdensEps=lpdensEps,
		# 										ceterisParibusPars=treatVarsAsIndep,
		# 										tol=rangeTol,max=F,idcToMod=idcToMod,
		# 										parscale=parscale.parvect,
		# 										bounds=parBounds,
		# 										trace = 9,
		# 										niter=1e2)
		# }
		# #max bound
		# findDensValBorder(test.i,
		# 									parVect=parVect,lpdensEps=lpdensEps,
		# 									ceterisParibusPars=treatVarsAsIndep,
		# 									tol=rangeTol,max=T,idcToMod=idcToMod,
		# 									parscale=parscale.parvect,
		# 									bounds=parBounds,
		# 									trace = 9,
		# 									niter=1e2)
		
		## range find ####
		for(direction in c('Min','Max')){
			cat(sprintf('  determining %s par values...',tolower(direction)))
			clusterExport(cl,list('calDat','treatVarsAsIndep'))
			border.coefs[which(notDeterminedBorders[,direction]),direction] <- 
				unlist(parLapplyLB(cl,which(notDeterminedBorders[,direction]),findDensValBorder,
													 parVect=parVect,lpdensEps=lpdensEps,
													 ceterisParibusPars=treatVarsAsIndep,
													 tol=rangeTol,max=(direction=='Max'),idcToMod=idcToMod,
													 parscale=parscale.parvect,
													 bounds=parBounds,
													 niter=1e3,# set niter so that the errors at least in the indep case are small
													 workerStagger = T)) 
			names(border.coefs[,direction]) <- names(parVect)
			# fallback values in case borders could not be determined:
			notDeterminedBorders[,direction] <- (is.infinite(border.coefs[,direction])+(parVect==border.coefs[,direction]))>=1
			border.coefs[,direction][notDeterminedBorders[,direction]] <- sampleParms[[direction]][notDeterminedBorders[,direction]]
			cat(sprintf('done. %i failures\n',sum(notDeterminedBorders[,direction])))
			write.csv(notDeterminedBorders,file.path(location.output,'notDeterminedBorders.csv'))
			# check that the min val actually has the desired like
			# this check only works for the independent case, as we do not retain the information
			# what the values of the other parameters where during range finding
			
			cat('\nsaving...')
			if(direction=='Min'){
				sampleParms[[direction]] <- pmax(border.coefs[,direction],parBounds[,'Min'])
			} else {
				sampleParms[[direction]] <- pmin(border.coefs[,direction],parBounds[,'Max'])
			}
			write.csv(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.csv'))
			saveRDS(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.RDS'))
			cat('done\n')	
		}
	}
	## make borders symmetric ####
	if(symmetricRanges%in%c('Max','Min')){
		cat('Symmetrifying parameter ranges\n')
		if(symmetricRanges=='Max'){
			sampleParms$distance <- pmax(sampleParms$Value-sampleParms$Min,
																	 sampleParms$Max-sampleParms$Value)
		} else {
			sampleParms$distance <- pmin(sampleParms$Value-sampleParms$Min,
																	 sampleParms$Max-sampleParms$Value)
		}
		# those that would have a distance of zero, we do not reassign
		sampleParms$MaxAfterDet <- sampleParms$Max
		sampleParms$MinAfterDet <- sampleParms$Min
		if(allowAssymetricToAvoidZeroRanges){
			sampleParms$Max[sampleParms$distance!=0] <- sampleParms$Value[sampleParms$distance!=0]+sampleParms$distance[sampleParms$distance!=0]
			sampleParms$Min[sampleParms$distance!=0] <- sampleParms$Value[sampleParms$distance!=0]-sampleParms$distance[sampleParms$distance!=0]
		} else {
			sampleParms$Max <- sampleParms$Value+sampleParms$distance
			sampleParms$Min <- sampleParms$Value-sampleParms$distance
		}
	}
	## read manual borders ####
	if(ignoreParBounds){
		cat('Not reading manual ranges, as ignoreParBounds==TRUE\n')
	} else {
		manualBorders <- read.csv(file.path(location.frida.info,name.frida_external_ranges))
		cat(sprintf('applying manual ranges for %i parameters...',nrow(manualBorders)))
		if(nrow(manualBorders)>0){
			for(r.i in 1:nrow(manualBorders)){
				sp.i <- which(sampleParms$Variable==manualBorders$Variable[r.i])
				if(!is.na(manualBorders$Min[r.i])){
					sampleParms$Min[sp.i] <- border.coefs[sp.i,'Min'] <- manualBorders$Min[r.i]
					notDeterminedBorders[sp.i,'Min'] <- FALSE
				}
				if(!is.na(manualBorders$Max[r.i])){
					sampleParms$Max[sp.i] <- border.coefs[sp.i,'Max'] <- manualBorders$Max[r.i]
					notDeterminedBorders[sp.i,'Max'] <- FALSE
				}
			}
		}
		write.csv(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.csv'))
		saveRDS(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.RDS'))
		cat('done\n')
	}
	## check for errors at the borders ####
	if(checkBorderErrors || kickParmsErrorRangeDet){
		borderLogLikeError <- array(NA,dim=c(length(parVect),2))
		colnames(borderLogLikeError) <- c('Min','Max')
		parVect <- sampleParms$Value
		names(parVect) <- sampleParms$Variable
		for(direction in c('Min','Max')){
			if(treatVarsAsIndep){
				cat(sprintf('Checking for likelihood at %s failures...',tolower(direction)))
				borderLogLikeError[,direction] <- unlist(parLapplyLB(cl,1:length(parVect),rangeCheckFun,
																														 parVect=parVect,
																														 border.coefs=border.coefs[,direction],
																														 lpdensEps=lpdensEps))
			}
			sampleParms[[paste0(direction,'NotDeterminedBorder')]] <- notDeterminedBorders[,direction]
			sampleParms[[paste0(direction,'BorderLogLikeError')]] <- borderLogLikeError[,direction] 
			sampleParms[[paste0(direction,'BoundByAuthors')]] <- sampleParms[[direction]]==parBounds[,direction]
			if(kickParmsErrorRangeDet){
				sampleParms[[paste0(direction,'KickParmsErrorRangeDet')]] <- 
					abs(borderLogLikeError[,direction]) > kickParmsErrorRangeDet.tolerance
			} else {
				sampleParms[[paste0(direction,'KickParmsErrorRangeDet')]] <- FALSE
			}
			cat('done\n')	
		}
	} else {
		for(direction in c('Min','Max')){
			sampleParms[[paste0(direction,'NotDeterminedBorder')]] <- notDeterminedBorders[,direction]
			sampleParms[[paste0(direction,'BorderLogLikeError')]] <- NA
			sampleParms[[paste0(direction,'KickParmsErrorRangeDet')]] <- FALSE
			sampleParms[[paste0(direction,'BoundByAuthors')]] <- sampleParms[[direction]]==parBounds[,direction]
		}
	}
	write.csv(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.csv'))
	saveRDS(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.RDS'))
	
	
	# write to frida_info like file for comparison to input
	frida_info.toModify <- read.csv(file.path(location.frida.info,name.frida_info))
	frida_info.toModify$includedInSampleParms <- frida_info.toModify$Variable %in% sampleParms$Variable
	idcOfSampleParmsInFridaInfo <- c()
	for(p.i in 1:nrow(sampleParms)){
		idcOfSampleParmsInFridaInfo[p.i] <- which(frida_info.toModify$Variable==sampleParms$Variable[p.i])
	}
	frida_info.toModify$newMin <- NA
	frida_info.toModify$newMin[idcOfSampleParmsInFridaInfo] <- sampleParms$Min
	frida_info.toModify$newMax <- NA
	frida_info.toModify$newMax[idcOfSampleParmsInFridaInfo] <- sampleParms$Max
	frida_info.toModify$beforeSymMin <- NA
	frida_info.toModify$beforeSymMin[idcOfSampleParmsInFridaInfo] <- sampleParms$MinAfterDet
	frida_info.toModify$beforeSymMax <- NA
	frida_info.toModify$beforeSymMax[idcOfSampleParmsInFridaInfo] <- sampleParms$MaxAfterDet
	frida_info.toModify$llikeErrorAtNewMin <- NA
	frida_info.toModify$llikeErrorAtNewMin[idcOfSampleParmsInFridaInfo] <- sampleParms$MinBorderLogLikeError
	frida_info.toModify$llikeErrorAtNewMax <- NA
	frida_info.toModify$llikeErrorAtNewMax[idcOfSampleParmsInFridaInfo] <- sampleParms$MaxBorderLogLikeError
	frida_info.toModify$MinBoundByAuthors <- NA
	frida_info.toModify$MinBoundByAuthors[idcOfSampleParmsInFridaInfo] <- sampleParms$MinBoundByAuthors
	frida_info.toModify$MaxBoundByAuthors <- NA
	frida_info.toModify$MaxBoundByAuthors[idcOfSampleParmsInFridaInfo] <- sampleParms$MaxBoundByAuthors
	frida_info.toModify$MinNotDetermined <- NA
	frida_info.toModify$MinNotDetermined[idcOfSampleParmsInFridaInfo] <- sampleParms$MinNotDeterminedBorder
	frida_info.toModify$MaxNotDetermined <- NA
	frida_info.toModify$MaxNotDetermined[idcOfSampleParmsInFridaInfo] <- sampleParms$MaxNotDeterminedBorder
	frida_info.toModify$MinKickedParmsErrorRangeDet <- NA
	frida_info.toModify$MinKickedParmsErrorRangeDet[idcOfSampleParmsInFridaInfo] <- sampleParms$MinKickParmsErrorRangeDet
	frida_info.toModify$MaxKickedParmsErrorRangeDet <- NA
	frida_info.toModify$MaxKickedParmsErrorRangeDet[idcOfSampleParmsInFridaInfo] <- sampleParms$MaxKickParmsErrorRangeDet
	write.csv(frida_info.toModify,file.path(location.output,'frida_info_ranged.csv'))

	# Kick out parameters with errors in the range determination and kickParmsErrorRangeDet was true
	if(kickParmsErrorRangeDet){
		cat(sprintf('Kicking out %i parameters for errors in range determination\n',
								sum(sampleParms$MinKickParmsErrorRangeDet|sampleParms$MaxKickParmsErrorRangeDet)))
		if(sum(sampleParms$MinKickParmsErrorRangeDet|sampleParms$MaxKickParmsErrorRangeDet)==nrow(sampleParms)){
			stop('would kick out all parms\n')
		}
		if(length(which(sampleParms$MinKickParmsErrorRangeDet))>0){
			sampleParms <- sampleParms[-which(sampleParms$MinKickParmsErrorRangeDet),]
		}
		if(length(which(sampleParms$MaxKickParmsErrorRangeDet))>0){
			sampleParms <- sampleParms[-which(sampleParms$MaxKickParmsErrorRangeDet),]
		}
		write.csv(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.csv'))
		saveRDS(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.RDS'))
	}
	
	# Sample the Parmeter Space ####
	parVect <- sampleParms$Value
	names(parVect) <- sampleParms$Variable
	maxLLike <- -negLLike(parVect)
	if(-baseNegLL!=maxLLike){stop('call ghostbusters\n')}
	
	## sample points ####
	# add the integer parms back
	sampleParms <- prepareSampleParms(sampleParms = sampleParms,integerParms = integerParms)
	samplePoints <- generateSobolSequenceForSampleParms(sampleParms,numSample,
																											restretchSamplePoints,
																											ignoreExistingResults = redoAllCalc,
																											integerParms = integerParms)
	if(ncol(samplePoints) != nrow(sampleParms) || nrow(samplePoints)!=numSample){
		cat('Invalid sample points regenerating\n')
		samplePoints <- generateSobolSequenceForSampleParms(sampleParms,numSample,
																												restretchSamplePoints,
																												ignoreExistingResults = T,
																												integerParms = integerParms)
	}
	
	write.csv(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.csv'))
	saveRDS(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.RDS'))
	saveRDS(samplePoints,file.path(location.output,'samplePoints.RDS'))
	# write.csv(samplePoints,file.path(location.output,'samplePoints.csv'))
	
	
	## write export spec ####
	extraVarNamesForExport <- read.csv(file.path(location.frida.info,name.frida_extra_variables_to_export_list))$FRIDA.FQN
	extraVarNamesForExport <- extraVarNamesForExport[nchar(extraVarNamesForExport)>4]
	writeFRIDAExportSpec(varsForExport.fridaNames = unique(c(varsForExport.fridaNames.orig,extraVarNamesForExport)),
											 location.frida)
	
	## evaluate sample points ####	
	logLikes <- clusterRunFridaForSamplePoints(samplePoints,chunkSizePerWorker,
																						 calDat=calDat,
																						 resSigma=resSigma,
																						 location.output=file.path(location.output,'detectedParmSpace'),
																						 redoAllCalc=redoAllCalc,
																						 plotDatWhileRunning=F,
																						 plotDatPerChunWhileRunning=plotDatPerChunWhileRunning,
																						 baseLL=-baseNegLL)
	logLikes[logLikes==-Inf] <- -.Machine$double.xmax
	if(plotWhileRunning){
		plotCape <- capabilities()
		if(!(plotCape['X11']|plotCape['aqua'])){
			pdf(file.path(location.output,'detectedParmSpace',paste0('logLikesDensity-',iterationNewMax,'.pdf')),
					width=10,	height=10)
		}
		histDat <- hist(logLikes,plot=F)
		plot(0,type='n',
				 main='Distribution log likelihoods of sample points',
				 xlab='log likelihood',
				 xlim=c(min(logLikes,maxLLike),max(logLikes,maxLLike)),
				 ylim=c(-1,max(histDat$counts)*1.04),
				 yaxs='i')
		box(col='gray')
		plot(histDat,add=T)
		abline(v=maxLLike,col='red')
		mtext('Vertical line is log likelihood of the best guess',3,0.1)
		if(!(plotCape['X11']|plotCape['aqua'])){
			dev.off()
		} else {
			dev.print(pdf,width=10,
								height=10,
								unit='cm',res=150,
								file.path(location.output,'detectedParmSpace',paste0('logLikesDensity-',iterationNewMax,'.pdf')))
		}
	}
	
	maxInd <- which.max(logLikes)
	if(logLikes[maxInd] > maxLLike){
		parVect <- samplePoints[maxInd,]
		newMaxFound <- T
		redoAllCalc <- F
		skipParMLE <- F
		cat('Found greater likelihood pars in sampling, rerunning with fit procedure\n')
	} else {
		cat('No greater likelihood found in sampling.\n')
		newMaxFound <- F	
	}
	
	# disable looping
	newMaxFound <- F
}
