cat('\nRunning runMLEandParmSpace.R\n\n')

source('initialise.R')

# config ####
cat('Config...')
source('config.R')
sink(file.path(location.output,'log.txt'),append=T)
cat(paste0(
	'\n###############################################################\n',
	format(Sys.time(), "%c"),
	location.output,
	'\n###############################################################\n'))
sink()
sink(file.path(location.output,'log.txt'),append=T,split=T)
# stop right here if anything the config specifies is missing, rather than
# silently running with e.g. the policy file of whatever ran here before
source('configValidator.R')
source('runInitialiseData.R')
continue <- readline(paste0('Output location created. Move any files to be used here\n',
														location.output,'\nHit ENTER when done.\n'))
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
# The parameters whose range is handed to us in frida_external_ranges.csv. Read
# here rather than where the ranges are applied, because the parscale
# determination and the range determination both skip these and both run first.
if(ignoreParBounds||forceParBounds){
	externalRanges <- data.frame(Variable=character(0),Min=numeric(0),Max=numeric(0))
} else {
	externalRanges <- read.csv(file.path(location.frida.info,name.frida_external_ranges))
}
externalRangeParmNames <- externalRanges$Variable[externalRanges$Variable%in%sampleParms$Variable]

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
# An externally ranged parameter never has its border determined, so the only
# thing left that would use its parscale is the MLE optimisation. Determining one
# costs a sweep over every order of magnitude, twice over for the ones that
# cannot be determined at all, so when we are not optimising we do not.
parscaleSkip <- rep(FALSE,length(jParVect))
names(parscaleSkip) <- names(jParVect)
if(skipParMLE){
	# index over the whole of jParVect. A logical index only as long as parVect
	# would be recycled over the resSigmaVect tail and skip covariance entries
	# that have nothing to do with any external range.
	parscaleSkip[names(jParVect)%in%externalRangeParmNames] <- TRUE
}

# start cluster ####
source('clusterHelp.R')
gobble <- clusterEvalQ(cl,source(file.path(baseWD,'config.R')))

# Development instrumentation, see developmentTools/countModelRuns.R. Inert
# unless FRIDA_COUNT_MODEL_RUNS is set in the environment, and harmless in a
# checkout that has no developmentTools directory: the marks below become no-ops.
if(identical(toupper(Sys.getenv('FRIDA_COUNT_MODEL_RUNS')),'TRUE')&&
	 file.exists('developmentTools/countModelRuns.R')){
	source('developmentTools/countModelRuns.R')
	devToolsCountModelRuns(cl)
} else {
	devToolsMarkSection <- function(...){invisible(NULL)}
	devToolsReportModelRuns <- function(...){invisible(NULL)}
}

# MLE and Sensi Loop ####
# the run specific config copies in the output directories predate this flag
if(!exists('redoFailedParscales')){
	redoFailedParscales <- F
}
# likewise, and these default to what the border search did before the knobs
# existed, so an old config keeps its old behaviour
if(!exists('rangeRootTol')){
	rangeRootTol <- NA
}
if(!exists('rangeRootMaxIter')){
	rangeRootMaxIter <- 1e3
}
# What a determination depends on, so a cache written from something else is not
# silently reused. Rebuilt where it is needed rather than held in a variable,
# because jParVect can be rebuilt under kickParmsParScaleDet.
funCurrentDeterminationKey <- function(baseNegLL=NULL){
	funDeterminationKey(location.frida,location.frida.info,name.frida_info,
											calDat,resSigma,names(jParVect),
											baseNegLL=baseNegLL,
											settings=list(treatVarsAsIndep=treatVarsAsIndep,
													 likeCutoffRatio=likeCutoffRatio,
													 rangeTol=rangeTol,
													 ignoreParBounds=ignoreParBounds,
													 forceParBounds=forceParBounds,
													 rangeRootTol=rangeRootTol,
													 rangeRootMaxIter=rangeRootMaxIter))
}
# A parscale a previous run could not determine is a result to keep, not work to
# redo. See funReadCachedParscale.
if(redoAllCalc){
	parscale <- rep(NA_real_,length(jParVect))
	names(parscale) <- names(jParVect)
	parscaleCachedNotDetermined <- rep(FALSE,length(jParVect))
	names(parscaleCachedNotDetermined) <- names(jParVect)
} else {
	parscale.cached <- funReadCachedParscale(file.path(location.output,'parscale.RDS'),
																					 names(jParVect),
																					 redoFailedParscales=redoFailedParscales,
																					 currentKey=funCurrentDeterminationKey())
	parscale <- parscale.cached$parscale
	parscaleCachedNotDetermined <- parscale.cached$notDetermined
	if(sum(parscaleCachedNotDetermined)>0){
		cat(sprintf('%i parameters had no determinable parscale last time, keeping that result. Set redoFailedParscales to retry them.\n',
								sum(parscaleCachedNotDetermined)))
	}
}
ordersOfMagGuesses.parvect <- funOrderOfMagnitude(sampleParms$Max-sampleParms$Min)
names(ordersOfMagGuesses.parvect) <- sampleParms$Variable
ordersOfMagGuesses.resSigmaVect <- funOrderOfMagnitude(resSigmaVect)-6
ordersOfMagGuesses <- c(ordersOfMagGuesses.parvect,ordersOfMagGuesses.resSigmaVect)

# used by the funFindParScale function ####
# The fallback sweep, the one a parameter gets when the guess around its own order
# of magnitude found nothing. It used to be a single global range, the lowest guess
# anywhere minus two up to the highest guess anywhere plus four, which mixes the
# sampled parameters against the residual variances and comes out 37 orders wide
# for both. The parameters that reach this sweep are the ones that fail, so they
# walk all 37.
#
# Per parameter instead. A step size larger than the range the parameter is
# sampled over is not a meaningful answer, so the sweep does not go above the
# order of that range, which is where the guess pass already ended. It extends
# downwards, which is where a finer scale might still be found. The residual
# variances have no author range to cap against and keep the old headroom.
ordersOfMagLimits <- array(NA_real_,dim=c(length(jParVect),2),
													 dimnames=list(names(jParVect),c('min','max')))
ordersOfMagLimits.parvectIdc <- 1:nrow(sampleParms)
ordersOfMagLimits.resSigmaIdc <- (nrow(sampleParms)+1):length(jParVect)
# a parameter whose Max equals its Min has no order of magnitude, and would drag
# every other parameter's floor to -Inf with it
ordersOfMagLimits.finite <- function(x){
	x <- x[is.finite(x)]
	if(length(x)==0){
		return(NA_real_)
	}
	return(x)
}
ordersOfMagLimits[ordersOfMagLimits.parvectIdc,'min'] <-
	min(ordersOfMagLimits.finite(ordersOfMagGuesses.parvect))-2
ordersOfMagLimits.maxParvect <- ordersOfMagGuesses.parvect+1
# A parameter whose author range has zero width has no order of magnitude of its
# own to cap against, so it keeps the ceiling the global range used to give it
# rather than being declared undeterminable by the bookkeeping.
ordersOfMagLimits.maxParvect[!is.finite(ordersOfMagLimits.maxParvect)] <-
	max(ordersOfMagLimits.finite(ordersOfMagGuesses))+4
ordersOfMagLimits[ordersOfMagLimits.parvectIdc,'max'] <- ordersOfMagLimits.maxParvect
ordersOfMagLimits[ordersOfMagLimits.resSigmaIdc,'min'] <-
	min(ordersOfMagLimits.finite(ordersOfMagGuesses.resSigmaVect))-2
ordersOfMagLimits[ordersOfMagLimits.resSigmaIdc,'max'] <-
	max(ordersOfMagLimits.finite(ordersOfMagGuesses.resSigmaVect))+4
cat(sprintf('Fallback parscale sweep is %.0f orders wide on average, was %.0f when it was one global range.\n',
						mean(ordersOfMagLimits[,'max']-ordersOfMagLimits[,'min']+1,na.rm=TRUE),
						max(ordersOfMagLimits.finite(ordersOfMagGuesses))+4-
							(min(ordersOfMagLimits.finite(ordersOfMagGuesses))-2)+1))
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
	
	# The determination a previous run left behind, if it is one we can use. Read
	# before the branch rather than inside it: a cache we cannot use has to fall
	# through to redetermining, and R cannot fall out of a branch it has taken.
	sampleParms.cached <- NULL
	if(!redoAllCalc&&!forceParBounds&&
		 file.exists(file.path(location.output,'sampleParmsParscaleRanged.RDS'))){
		sampleParms.cached <- funReadCachedRangedSampleParms(
			file.path(location.output,'sampleParmsParscaleRanged.RDS'),
			currentKey=funCurrentDeterminationKey(baseNegLL=baseNegLL))
	}
	if(forceParBounds){
		cat('Forced using frida_info bounds\n')
	} else if(!is.null(sampleParms.cached)){
		cat('loading existing sampleParmsParscaleRanged\n')
		sampleParms <- prepareSampleParms(excludeNames=excludedParmsForBeingIntegers,
																		sampleParms = sampleParms.cached)
		# What follows works from the determination results, not from the ranges a
		# previous run derived from them. Putting Min and Max back to what the
		# determination produced is what lets the external range overrides and the
		# symmetrification run again against the config as it is now, rather than being
		# applied a second time on top of themselves.
		sampleParms$Min <- sampleParms$MinAfterDet
		sampleParms$Max <- sampleParms$MaxAfterDet
		parVect <- sampleParms$Value
		names(parVect) <- sampleParms$Variable 
		jParVect <- c(parVect,resSigmaVect)
		parscale.parvect <- sampleParms$parscale
		parscaleSkip.parvect <- sampleParms$parscaleStatus=='skippedExternalRange'
		parscaleNotDetermined.parVect <- sampleParms$parscaleStatus=='notDetermined'
	} else {
		# determine parscale ####
		devToolsMarkSection('parscale determination',cl)
		cat('Determining parscales...\n')
		if(sum(parscaleSkip)>0){
			cat(sprintf('Skipping parscale determination for %i parameters with external ranges.\n',
									sum(parscaleSkip)))
		}
		iterations <- 0
		parallelParscale <- T
		useOrdersOfMagGuesses <- T
		# Hoisted out of the loop below: none of this changes between the two passes,
		# and under a psock cluster every export serialises calDat to every worker.
		if(parallelParscale){
			clusterExport(cl,list('baseNegLL',
														'ordersOfMagLimits','responseTolerance',
														'orderOfMagNegLLErrorFun','funFindParScale',
														'jnegLLikelihood.f','ordersOfMagGuesses',
														'calDat','resSigma',
														'jParVect'))
			gobble <- clusterEvalQ(cl,source(file.path(baseWD,'funParmSpace.R')))
		}
		while(iterations < 2 && sum((is.na(parscale)|is.infinite(parscale))&!parscaleSkip&!parscaleCachedNotDetermined)>0){
			parsToDet <- which((is.na(parscale)|is.infinite(parscale))&!parscaleSkip&!parscaleCachedNotDetermined)
			cat(sprintf('Determining the parscale of %i parameters. %i parameters with already known parscale.%s\n',
									length(parsToDet),length(parscale)-length(parsToDet)-sum(parscaleSkip),
									if(useOrdersOfMagGuesses){' Using guesses.'}else{' Not using guesses.'}))
			if(parallelParscale){
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
		parscaleSkip.parvect <- parscaleSkip[1:nrow(sampleParms)]
		cat('done\n')
		
		## check for bad behaviour in parscale ####
		# only the entries in parVect can be excluded. The entries in resSigmaVect need to 
		# be delt with. E.g. by using the guess values. The maximum likelihood vars (diag
		# elements of the covmat can always be determined as the variance of those obs.
		# A parameter whose determination we skipped is not a problem case, we already
		# know what its range is going to be.
		problemCases <- (is.infinite(parscale)|is.na(parscale))&!parscaleSkip
		parscaleNotDetermined.parVect <- problemCases[1:nrow(sampleParms)]
		problemCasesIdc.parVect <- which(parscaleNotDetermined.parVect)
		problemCasesIdc.resSigmaVect <- which(problemCases[(nrow(sampleParms)+1):length(jParVect)])
		cat(sprintf('%i parscales could not be determined.\n',sum(problemCases)))
		if(length(problemCasesIdc.resSigmaVect)>0){
			cat(sprintf('  %i in resSigmaVect, guesses will be used\n',
									length(problemCasesIdc.resSigmaVect)))
			parscale.resSigmaVect[problemCasesIdc.resSigmaVect] <- 
				10^ordersOfMagGuesses.resSigmaVect[problemCasesIdc.resSigmaVect]
		} else {
			cat('  No problem cases in resSigmaVect\n')
		}
		# Why each parameter has the parscale it has, kept on sampleParms because a
		# later run that reuses this determination has to rebuild the same distinction
		# and the lines printed here are long gone by then.
		sampleParms$parscaleStatus <- ifelse(parscaleSkip.parvect,'skippedExternalRange',
																				 ifelse(parscaleNotDetermined.parVect,'notDetermined',
																				 			 'determined'))
		scaleErrorParmNames <- sampleParms$Variable[problemCasesIdc.parVect]
		if(length(problemCasesIdc.parVect)>0){
			if(kickParmsParScaleDet){
				cat(sprintf('  %i in parVect, these parms will be dropped\n',length(problemCasesIdc.parVect)))
				cat(paste(scaleErrorParmNames,collapse='\n'))
				cat('\n')
				exclusionList <- data.frame(excludedName=scaleErrorParmNames)
				write.csv(exclusionList,file.path(location.output,name.frida_parameter_exclusion_list))
				# everything indexed against parVect has to lose the same entries
				keptIdc <- which(!parscaleNotDetermined.parVect)
				parscale.parvect <- parscale.parvect[keptIdc]
				parscaleSkip.parvect <- parscaleSkip.parvect[keptIdc]
				parscaleNotDetermined.parVect <- parscaleNotDetermined.parVect[keptIdc]
				sampleParms <- prepareSampleParms(excludeNames = c(scaleErrorParmNames,excludedParmsForBeingIntegers))
				sampleParms$parscaleStatus <- ifelse(parscaleSkip.parvect,'skippedExternalRange','determined')
				parVect <- sampleParms$Value
				names(parVect) <- sampleParms$Variable 
				jParVect <- c(parVect,resSigmaVect)
			} else {
				# They are kept and sampled over the ranges their authors gave them in
				# frida_info, the same fallback a failed border determination gets. This
				# used to drop them regardless of the setting, while saying it did not.
				cat(sprintf('  %i in parVect, these parms keep the ranges their authors gave them in frida_info.\n',length(problemCasesIdc.parVect)))
				cat(paste(scaleErrorParmNames,collapse='\n'))
				cat('\n')
			}
		}
		# The parameters we determine no range for still need a parscale entry, the
		# vector is indexed positionally against parVect everywhere it is used and a
		# short one would silently misalign. The order of magnitude of the author range
		# is the same stand in the resSigmaVect problem cases above get, and it keeps
		# the vector numeric for the optimisation.
		parscale.fillIdc <- which(parscaleNotDetermined.parVect|parscaleSkip.parvect)
		if(length(parscale.fillIdc)>0){
			parscale.parvect[parscale.fillIdc] <- 
				10^funOrderOfMagnitude(sampleParms$Max-sampleParms$Min)[parscale.fillIdc]
		}
		parscale.all <- parscale
		parscale <- c(parscale.parvect,parscale.resSigmaVect)
		
		## save parscale ####
		cat('saving ParScale...')
		# The status travels with the values, so a later run can tell a parameter that
		# was never tried from one that was tried and could not be determined.
		parscaleStatus.jParVect <- ifelse(parscaleSkip,'skippedExternalRange',
																			ifelse(problemCases,'notDetermined','determined'))
		# names off parscale.all, not off jParVect: kickParmsParScaleDet rebuilds
		# jParVect shorter above, while parscale.all still spans what was determined
		names(parscaleStatus.jParVect) <- names(parscale.all)
		saveRDS(list(parscale=parscale.all,
								 status=parscaleStatus.jParVect,
								 key=funCurrentDeterminationKey(baseNegLL=baseNegLL),
								 savedAt=Sys.time()),
						file.path(location.output,'parscale.RDS'))
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
	devToolsMarkSection('range finding',cl)
	## par bounds ####
	parBounds <- funParBoundsForSampleParms(sampleParms,frida_info)
	notDeterminedBorders <- array(TRUE,dim=c(length(parVect),2))
	colnames(notDeterminedBorders) <- c('Min','Max')
	# the borders themselves. Numeric from the start: every path below writes
	# numbers into it, and a matrix left logical would compare TRUE against
	# parameter values in the fallback test further down.
	border.coefs <- array(NA_real_,dim=dim(notDeterminedBorders),
											dimnames=dimnames(notDeterminedBorders))
	if(!is.null(sampleParms.cached)){
		# Reusing a determination means reusing what it found, including which borders
		# it could not determine. Left at the TRUE they are initialised to, every
		# parameter would look like a border failure here, and the flags written out
		# below would say so.
		notDeterminedBorders[,'Min'] <- sampleParms$MinNotDeterminedBorder
		notDeterminedBorders[,'Max'] <- sampleParms$MaxNotDeterminedBorder
		border.coefs[,'Min'] <- sampleParms$MinAfterDet
		border.coefs[,'Max'] <- sampleParms$MaxAfterDet
		rangeDetSkip <- parscaleNotDetermined.parVect|parscaleSkip.parvect
	}
	if(is.null(sampleParms.cached)||checkBorderErrors||kickParmsErrorRangeDet){
		# costs a frida run, so only for the searches that actually use it
		# kept, not just folded into lpdensEps: the border search needs it to skip the
		# probe at the starting point, which is the same for every parameter and both
		# directions
		lpdensAtParVect <- -negLLike(parVect)
		lpdensEps <- lpdensAtParVect - log(likeCutoffRatio)
	}
	if(forceParBounds){
		cat('Forcing coefs sample range to be equal tovalues frida_info\n')
		border.coefs <- as.matrix(sampleParms[,c('Min','Max')])
	} else if (is.null(sampleParms.cached)) {
		# minimize and maximize each parameter with others free, until density is 
		# equal to pdensEps
		cat('determining coef sample range...\n')
		# boundary value in log likelihood
		idcToMod <- 1:length(parVect)
		if(ignoreParBounds){
			parBounds <- array(c(rep(-.Machine$double.xmax,length(parVect)),
													 rep(.Machine$double.xmax,length(parVect))),
												 dim=c(length(parVect),2))
		}
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
		# A parameter without a parscale cannot have its border determined, the search
		# needs a scale to step with, and one with an external range has no reason to.
		# Both already have their answer, the range their authors gave them, which is
		# what the fallback below assigns. Marking them infinite here puts them through
		# that same fallback instead of a second code path.
		rangeDetSkip <- parscaleNotDetermined.parVect|parscaleSkip.parvect
		# Once, rather than once per direction. resSigma is in the list because the MLE
		# above rewrites it and the workers evaluate the likelihood against their own
		# copy, which was otherwise left at whatever clusterHelp.R last sent them.
		clusterExport(cl,list('calDat','resSigma','treatVarsAsIndep'))
		# Min and Max used to run as two pools, one after the other. That put a barrier
		# in the middle and left workers idle at the tail of each while a straggler
		# finished a search that ran to the iteration limit. The two directions are
		# independent, so they go into one pool of every border to be found.
		borderTasks <- list()
		for(direction in c('Min','Max')){
			for(td in which(notDeterminedBorders[,direction]&!rangeDetSkip)){
				borderTasks[[length(borderTasks)+1]] <-
					list(parIdx=td,max=(direction=='Max'),direction=direction)
			}
		}
		# Longest first. parLapplyLB balances the load, but it hands the tasks out in
		# the order it was given them, so a border that needs a thousand iterations can
		# be picked up last and hold the pool open on its own. How many parscale steps
		# lie between the value and the bound is the cost estimate available here
		# without paying for one: a border a long way out in units of the step size
		# takes more steps to reach. It is a heuristic, but the order it replaces is
		# parameter index, which is unrelated to cost.
		borderTaskCost <- function(tsk){
			scale <- parscale.parvect[tsk$parIdx]
			span <- abs(parBounds[tsk$parIdx,ifelse(tsk$max,2,1)]-parVect[tsk$parIdx])
			if(!is.finite(scale)||scale==0||!is.finite(span)){
				# unknown cost goes out early rather than last
				return(Inf)
			}
			return(span/scale)
		}
		if(length(borderTasks)>0){
			cat(sprintf('  determining %i borders over %i workers...',
									length(borderTasks),length(cl)))
			borderTasks <- borderTasks[order(sapply(borderTasks,borderTaskCost),
																			 decreasing=TRUE)]
			borderResults <- parLapplyLB(cl,borderTasks,funBorderTask,
																	 parVect=parVect,lpdensEps=lpdensEps,
																	 ceterisParibusPars=treatVarsAsIndep,
																	 tol=rangeTol,idcToMod=idcToMod,
																	 parscale=parscale.parvect,
																	 bounds=parBounds,
																	 niter=1e3,# set niter so that the errors at least in the indep case are small
																	 rootTolFactor=rangeRootTol,
																	 rootMaxIter=rangeRootMaxIter,
																	 lpdensAtParVect=lpdensAtParVect,
																	 workerStagger = T)
			for(task.i in seq_along(borderTasks)){
				border.coefs[borderTasks[[task.i]]$parIdx,borderTasks[[task.i]]$direction] <-
					as.numeric(borderResults[[task.i]])
			}
			cat('done\n')
		}
		for(direction in c('Min','Max')){
			border.coefs[rangeDetSkip,direction] <- Inf
			names(border.coefs[,direction]) <- names(parVect)
			# fallback values in case borders could not be determined:
			notDeterminedBorders[,direction] <- 
				(is.infinite(border.coefs[,direction])+(parVect==border.coefs[,direction]))>=1
			border.coefs[,direction][notDeterminedBorders[,direction]] <- 
				sampleParms[[direction]][notDeterminedBorders[,direction]]
			cat(sprintf('  %s: %i determined, %i failed, %i skipped (no parscale), %i skipped (external range)\n',
									tolower(direction),
									sum(!notDeterminedBorders[,direction]),
									sum(notDeterminedBorders[,direction]&!rangeDetSkip),
									sum(parscaleNotDetermined.parVect),
									sum(parscaleSkip.parvect)))
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
	## record the determined borders ####
	sampleParms$MaxAfterDet <- sampleParms$Max
	sampleParms$MinAfterDet <- sampleParms$Min
	## read manual borders ####
	# Applied before the symmetrification rather than after it, so that
	# symmetrifyExternalRanges is free to decide whether an external range gets
	# symmetrified. Applying them afterwards, as this used to, made that decision
	# for us: an external range could never be symmetrified.
	if(ignoreParBounds || forceParBounds){
		cat('Not reading manual ranges, as ignoreParBounds||forceParBounds==TRUE\n')
	} else {
		manualBorders <- externalRanges[externalRanges$Variable %in% sampleParms$Variable,]
		cat(sprintf('applying manual ranges for %i parameters...',nrow(manualBorders)))
		if(nrow(manualBorders)>0){
			for(r.i in 1:nrow(manualBorders)){
				sp.i <- which(sampleParms$Variable==manualBorders$Variable[r.i])
				if(length(sp.i)==1){
					if(!is.na(manualBorders$Min[r.i])){
						sampleParms$Min[sp.i] <- border.coefs[sp.i,'Min'] <- manualBorders$Min[r.i]
						notDeterminedBorders[sp.i,'Min'] <- FALSE
					}
					if(!is.na(manualBorders$Max[r.i])){
						sampleParms$Max[sp.i] <- border.coefs[sp.i,'Max'] <- manualBorders$Max[r.i]
						notDeterminedBorders[sp.i,'Max'] <- FALSE
					}
				} else if (length(sp.i)>1){
					stop('Multiple parms with same name\n')
				}
				# if the parm is not present at all do nothing
			}
		}
		write.csv(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.csv'))
		saveRDS(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.RDS'))
		cat('done\n')
	}
	## make borders symmetric ####
	sampleParms <- funSymmetrifyRanges(sampleParms,parBounds,notDeterminedBorders,
																			 externalRangeParmNames=externalRangeParmNames,
																			 symmetricRanges=symmetricRanges,
																			 allowAssymetricToAvoidZeroRanges=allowAssymetricToAvoidZeroRanges,
																			 symmetricRangesBoundByAuthors=symmetricRangesBoundByAuthors,
																			 symmetrifyExternalRanges=symmetrifyExternalRanges,
																			 symmetrifyFallbackAuthorRanges=symmetrifyFallbackAuthorRanges)
	## check for errors at the borders ####
	devToolsMarkSection('border checks',cl)
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
			sampleParms[[paste0(direction,'BorderLogLikeError')]] <- rep(NA,nrow(sampleParms))
			sampleParms[[paste0(direction,'KickParmsErrorRangeDet')]] <- rep(FALSE,nrow(sampleParms))
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
	frida_info.toModify$parscaleStatus <- NA
	frida_info.toModify$parscaleStatus[idcOfSampleParmsInFridaInfo] <- sampleParms$parscaleStatus
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
	
	# The determination is complete, so the cache it leaves behind can be keyed.
	# Written here and nowhere else on purpose: the intermediate saves above happen
	# while the determination is still running, and a run interrupted among them
	# should leave a cache that fails this check rather than one that looks whole.
	saveRDS(funCurrentDeterminationKey(baseNegLL=baseNegLL),
					file.path(location.output,'sampleParmsParscaleRanged.key.RDS'))

	devToolsMarkSection('determination done',cl)
	devToolsReportModelRuns(cl)

	# Sample the Parmeter Space ####
	parVect <- sampleParms$Value
	names(parVect) <- sampleParms$Variable
	maxLLike <- -negLLike(parVect)
	if(-baseNegLL!=maxLLike) {
		cat(sprintf('Would have called ghostbusters... -baseNegLL=%10f, maxLLike=%10f\n', -baseNegLL, maxLLike))
	}
	
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
	# apply baseline parms
	if(!is.na(name.baselineParmFile)&&name.baselineParmFile!=''){
		baselineParms <- read.csv(file.path(location.frida.configs,name.baselineParmFile))
		baselineParmsNames <- unname(unlist(read.table(file.path(location.frida.configs,name.baselineParmFile),nrows = 1,sep=',')))
		baselineParmsNames <- gsub('([^\\]])$','\\1\\[1\\]',baselineParmsNames,perl=T)
		baselineParmsNames <- gsub('\\[\\*','\\[1',baselineParmsNames,perl=T)
		if(nrow(baselineParms)>1){
			stop('baseline parms may only be a single set (line)\n')
		}
		for(par.i in 1:length(colnames(baselineParms))){
			if(!baselineParmsNames[par.i]%in%colnames(samplePoints)){
				newcol <- array(baselineParms[par.i],dim=c(nrow(samplePoints),1))
				colnames(newcol) <- baselineParmsNames[par.i]
				samplePoints <- base::cbind(samplePoints,newcol)
			}
		}
	}
	cat('saving sampleParms and samplePoints...')
	write.csv(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.csv'))
	saveRDS(sampleParms,file.path(location.output,'sampleParmsParscaleRanged.RDS'))
	saveRDS(samplePoints,file.path(location.output,'samplePoints.RDS'))
	write.csv(samplePoints,file.path(location.output,'samplePoints.csv'))
	cat('done\n')
	
	## write export spec ####
	extraVarNamesForExport <- read.csv(file.path(location.frida.info,name.frida_extra_variables_to_export_list))$FRIDA.FQN
	extraVarNamesForExport <- extraVarNamesForExport[nchar(extraVarNamesForExport)>4]
	writeFRIDAExportSpec(varsForExport.fridaNames = unique(extraVarNamesForExport),
											 location.frida)
	source('clusterHelp.R') #make sure the workers also get the updated export spec
	
	## evaluate sample points ####	
	logLikes <- clusterRunFridaForSamplePoints(samplePoints,chunkSizePerWorker,
																						 calDat=calDat,
																						 resSigma=resSigma,
																						 location.output=file.path(location.output,'detectedParmSpace'),
																						 redoAllCalc=redoAllCalc,
																						 plotDatWhileRunning=F,
																						 plotDatPerChunWhileRunning=plotDatPerChunWhileRunning,
																						 baseLL=-baseNegLL)
	logLikes[logLikes==-Inf] <- logLike.failedRun
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

# stop cluster ####
# The sample point evaluation at the end of the loop above is the last thing that uses
# the cluster. The scripts that run after this one share this R session, so a cluster
# left running here would keep its workers alive for the rest of the job, and
# setupTMPFS.R would take a still existing cl as the sign that it has nothing to do.
if(exists('cl')){
	cat('stopping cluster...')
	tryCatch(stopCluster(cl),error=function(e){})
	rm(cl)
	cat('done\n')
}
