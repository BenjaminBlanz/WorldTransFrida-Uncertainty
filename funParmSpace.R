# The following functions work together to 
# find the borders of the density where density is equal to pdensEps
# repeat two steps, find root in current parm
# maximize density for root of current parm using other parms
# This function relies on the negLLike function being present in the global env.
# This allows the user to specify the type of likelihood function.

# likelihood ####
jnegLLikelihood.f <- function(jParVect){
	parVect <- jParVect[1:nrow(sampleParms)]
	if(treatVarsAsIndep){
		resSigma <- diag(jParVect[(nrow(sampleParms)+1):length(jParVect)])
	} else {
		resSigma <- array(NA,dim=rep(ncol(calDat),2))
		resSigma[!lower.tri(resSigma)]<- jParVect[(nrow(sampleParms)+1):length(jParVect)]
		resSigma[lower.tri(resSigma)] <- t(resSigma)[lower.tri(resSigma)]
	}
	runDat <- runFRIDASpecParms(parVect)
	calDatInRunDat <- which(colnames(calDat)%in%colnames(runDat))
	if(length(calDatInRunDat)>0){
		resDat <- calDat[calDatInRunDat]-runDat[1:nrow(calDat),colnames(calDat)[calDatInRunDat]]
		lLikelihood <- funLogLikelihood(resDat,resSigma)
	} else {
		lLikelihood <- rep(1,ncol(runDat))
	}
	# A run that did not complete gets marked. Used when narrowing the parms space.
	# A run the model could not finish must never come back with a likelihood that 
	# looks like a good fit. Testing the last
	# row of runDat for an NA does not see those runs: stella writes a short
	# output file when a run stops early, and its last row is a perfectly good
	# year. 
	if(!funRunReachedFinalYear(runDat)){
		lLikelihood <- logLike.failedRun+(sum(!is.na(runDat[[1]]))*logLike.quasiEps)
	}
	return(-lLikelihood)
}
negLLike <- function(parVect){
	runDat <- runFRIDASpecParms(parVect)
	calDatInRunDat <- which(colnames(calDat)%in%colnames(runDat))
	if(length(calDatInRunDat)>0){
		resDat <- calDat[calDatInRunDat]-runDat[1:nrow(calDat),colnames(calDat)[calDatInRunDat]]
		lLikelihood <- funLogLikelihood(resDat,resSigma)
	} else {
		lLikelihood <- rep(1,ncol(runDat))
	}
	# A run that did not complete gets marked. Used when narrowing the parms space.
	# A run the model could not finish must never come back with a likelihood that 
	# looks like a good fit. Testing the last
	# row of runDat for an NA does not see those runs: stella writes a short
	# output file when a run stops early, and its last row is a perfectly good
	# year. 
	if(!funRunReachedFinalYear(runDat)){
		lLikelihood <- logLike.failedRun+(sum(!is.na(runDat[[1]]))*logLike.quasiEps)
	}
	return(-lLikelihood)
}
jnegLLikelihood.gr <- function(jParVect){
	nllike <- jnegLLikelihood.f(jParVect)
	grad <- c()
	for(i in 1:length(jParVect)){
		jParVectD <- jParVect
		jParVectD[i] <- jParVect[i]+parscale[i]/100
		#TODO: finish this
	}
}

likeGoalDiffFun <- function(par,parVect,parIdx,lpdensEps, ...){
	parVect[parIdx] <- par
	llike <- -negLLike(parVect, ...)
	return(llike-lpdensEps)
}

# density border search ####
densMaxGivenParFun <- function(otherPars,parVect,parIdx,idcToMod, ...){
	parVect[setdiff(idcToMod,parIdx)] <- otherPars
	return(negLLike(parVect, ...))
}
# if ceterisParibusPars is TRUE, the densValBorder is found for the selecteed par Idx,
# keeping all other pars at the values in parVect, i.e. this will not account
# for rotated elipsoid parameter distributions, but just for the slice through
# the likelihood at parVect.
# 
# idcToMod: Specify the indices that should be varied together with parIdx in the 
#           search. Defaults to all.
findDensValBorder <- function(parIdx,parVect,lpdensEps,ceterisParibusPars=F,
															maxiter=1e4,tol=1e-4,max=F,
															trace=0, idcToMod=1:length(parVect),
															parscale=rep(1,length(parVect)),
															bounds=NULL,
															niter=1000,
															rootTolFactor=NA,
															rootMaxIter=NULL,
															lpdensAtParVect=NULL,
															workerStagger=FALSE,
															...){
	if(workerStagger){
		# The stagger exists so the workers do not all reach for the filesystem in the
		# same instant when the pool starts.
		if(!exists('workerHasStaggered',envir=globalenv())){
			Sys.sleep(workerID*0.04)
			assign('workerHasStaggered',TRUE,envir=globalenv())
		}
	}
	if(length(parIdx)>1){
		stop('only one parIdx at a time\n')
	}
	if(is.list(idcToMod)){
		idcToMod <- idcToMod[[parIdx]]
	}
	# The objective is a stella run whose output is read back from a csv with some loss 
	# of precision. A tolerance scaled to the parameter's own step size is one the
	# model can actually resolve. rootTolFactor of NA gives an absolute 1e-16.
	rootTol <- if(is.finite(rootTolFactor)&&
								is.finite(parscale[parIdx])&&parscale[parIdx]!=0){
		abs(parscale[parIdx])*rootTolFactor
	} else {
		1e-16
	}
	if(is.null(rootMaxIter)){
		rootMaxIter <- niter
	}
	# The log density at the parameters as they came in. The caller measured it to
	# derive lpdensEps. Only usable while nothing in parVect has moved.
	parVectUnchanged <- TRUE
	fAtParVect <- if(is.null(lpdensAtParVect)){NULL}else{lpdensAtParVect-lpdensEps}
	likeGoalDiffAt <- function(x,parVect,parIdx,lpdensEps,...){
		if(!is.null(fAtParVect)&&parVectUnchanged&&
			 isTRUE(as.numeric(x)==as.numeric(parVect[parIdx]))){
			return(fAtParVect)
		}
		return(likeGoalDiffFun(x,parVect,parIdx,lpdensEps,...))
	}
	idcToMod.base <- idcToMod
	# One pass per widening of the set of parameters that move with parIdx. With a
	# single element there is nothing to widen to, but the border for parIdx itself
	# still has to be found, so the loop runs once.
	idcsToModSeq <- if(length(idcToMod.base)<2){
		seq_along(idcToMod.base)
	} else {
		2:length(idcToMod.base)
	}
	if(length(idcsToModSeq)==0){
		stop('idcToMod is empty, there is nothing to search\n')
	}
	for(idcsToMod.i in idcsToModSeq){
		idcToMod <- idcToMod.base[c(1:idcsToMod.i)]
		if(trace>0&&!ceterisParibusPars){
			cat('Running with idcToMod ',idcToMod,'\n')
		}
		par.val <- parVect[parIdx]
		likeAtMax  <- lpdensEps+5
		likeAtMaxOld  <- lpdensEps
		iter <- 1
		while(likeAtMax-lpdensEps > tol && 
					# likeAtMax-likeAtMaxOld > tol &&
					iter <= maxiter){
			# find root
			if(max){ # Maximizing
				if(!is.null(bounds)){
					if(par.val>=bounds[parIdx,2]){
						return(bounds[parIdx,2])
					}
					root.range <- c(par.val,bounds[parIdx,2])
					bound <- bounds[parIdx,2]
				} else {
					root.range <- c(par.val,par.val+abs(par.val)*10)
					bound <- NULL
				}
			} else { # Minimizing
				if(!is.null(bounds)){
					if(par.val<=bounds[parIdx,1]){
						return(bounds[parIdx,1])
					}
					root.range <- c(bounds[parIdx,1],par.val)
					bound <- bounds[parIdx,1]
				} else {
					root.range <- c(par.val-abs(par.val)*10,par.val)
					bound <- NULL
				}
			}
			# if there is no sign change between the endpoints of root.range, use secant's
			# method otherwise use uniroot
			# Both of these are stella runs, and both branches below would otherwise
			# evaluate the same two points a second time. uniroot takes them as
			# f.lower/f.upper. The secant branch can only reuse f.lo when maximising:
			# it rebuilds root.range first, and only in that direction does the new
			# root.range[1] stay equal to the point f.lo was measured at.
			f.lo <- likeGoalDiffAt(root.range[1],parVect,parIdx,lpdensEps,...)
			f.hi <- likeGoalDiffAt(root.range[2],parVect,parIdx,lpdensEps,...)
			if(((f.lo>0)-(f.hi>0))==0){
				if(max){
					root.range <- c(par.val,par.val+parscale[parIdx])
				} else {
					root.range <- c(par.val-parscale[parIdx],par.val)
				}
				par.val.old <-par.val
				par.val <- secant(likeGoalDiffFun,
													root.range[1],root.range[1]*1.001,
													f0=if(max){f.lo}else{NULL},
													parVect=parVect,
													parIdx=parIdx,
													lpdensEps=lpdensEps,
													doWarn = F,
													bound = bound,
													trace=trace,
													niter=niter,...)
				if(max){
					if(par.val < par.val.old){
						# hail mary
						par.val <- par.val.old
					}
				} else {
					par.val <- min(par.val,par.val.old)
				}
				if(is.infinite(par.val)){
					return(par.val)
				}
			} else {
				par.val <- suppressWarnings(uniroot(likeGoalDiffFun,
																						root.range,
																						f.lower = f.lo, f.upper = f.hi,
																						parVect=parVect,
																						parIdx=parIdx,
																						lpdensEps=lpdensEps,
																						tol = rootTol, maxiter = rootMaxIter,...)$root)
			}
			if(ceterisParibusPars){
				return(par.val)
			} else {
				parVect[parIdx] <- par.val
				parVectUnchanged <- FALSE
				if(trace>0){
					cat('iter ',iter,' ',par.val,' : ')
				}
				# maximize density at root using other parms
				otherIdx <- setdiff(idcToMod,parIdx)
				if(length(otherIdx)==0){
					# nothing to reoptimise against, so this is the ceteris paribus answer
					return(par.val)
				}
				otherPars <- parVect[otherIdx]
				res <- suppressWarnings(optimx(otherPars,densMaxGivenParFun,
																			 method = 'Nelder-Mead',
																			 control=list(dowarn = F,
																			 						 parscale=parscale[otherIdx]),
																			 parVect=parVect,
																			 parIdx=parIdx, 
																			 idcToMod=idcToMod,...))
				parVect[otherIdx] <- unlist(res[1:length(otherIdx)])
				parVectUnchanged <- FALSE
				likeAtMaxOld <- likeAtMax
				likeAtMax <- -res$value
				if(trace>0){
					cat(parVect,' ',likeAtMax-lpdensEps,'\n')
				}
				if(likeAtMax>likeAtMaxOld){
					# The step improved the likelihood, so the point we moved to sits above the
					# contour we are trying to find. Reoptimising every parameter brings it back
					# down.
					if(trace>0){
						cat('Likelihood Imporovement After Step Reoptimizing\n')
					}
					res <- suppressWarnings(optimx(parVect,negLLike,
					                               method = 'Nelder-Mead',
					                               control=list(dowarn = F,
					                                            parscale=parscale),...))
					parVect <- unlist(res[1:length(parVect)])
					parVectUnchanged <- FALSE
					par.val <- parVect[parIdx]
					likeAtMax <- -res$value
					if(trace>0){
						cat(parVect,' ',likeAtMax-lpdensEps,'\n')
					}
				}
			}
			iter <- iter+1
		} 
	}
	return(par.val)
}

funBorderTask <- function(task,...){
	findDensValBorder(task$parIdx,max=task$max,...)
}

# Every call to fun here is a stella run, so the values are carried rather than
# recomputed: x0 and x1 are the previous iteration's x1 and x2. The tolerance 
# test's x2 becomes the next x1. A caller that has already evaluated the starting 
# points can pass them as f0 and f1.
#
# The returned root carries the value of fun at that point as the attribute
# 'fval', so a caller does not have to evaluate it again. It is absent on the
# returns that never evaluated the point they hand back.
secant <- function(fun, x0, x1, tol=1e-07, niter=1e4, doWarn=T, trace=0,
									 bound=NULL,hasToBePositive=FALSE,f0=NULL,f1=NULL,...){
	if(is.null(bound)){
		bound <- sign(x1-x0)*Inf
	}
	withFval <- function(x,fval){
		attr(x,'fval') <- fval
		return(x)
	}
	if(is.null(f0)){
		f0 <- fun(x0,...)
	}
	if(is.null(f1)){
		f1 <- fun(x1,...)
	}
	for ( i in 1:niter ){
		x2 <- x1-f1*(x1-x0)/(f1-f0)
		if(trace>0){
			cat(sprintf('secant x0: %10f f0: %10f x1: %10f f1: %10f  x2: %10f\n',
									x0,f0,x1,f1,x2))
		}
		if(is.infinite(x2)||is.nan(x2)){
			return(bound)
		}
		if(hasToBePositive && x2 < 0){
			return(NA)
		}
		f2 <- fun(x2,...)
		if(abs(f2) < tol || abs(x2)>abs(bound)){
			return(withFval(x2,f2))
		}
		if(x0==x2){
			if(doWarn){
				warning("In secant cycle detected\n")
			}
			return(withFval(x2,f2))
		}
		x0 <- x1
		f0 <- f1
		x1 <- x2
		f1 <- f2
	}
	if(doWarn){
		warning("In secant exceeded allowed number of iterations\n")
	}
	return(withFval(x2,f2))
}

# parameter scale ####
rangeCheckFun <- function(rangeCheck.i,parVect,border.coefs,lpdensEps){
	cat(sprintf('\r%4i %100s',rangeCheck.i,names(parVect[rangeCheck.i])))
	parVectMinCheck.i <- parVect
	parVectMinCheck.i[rangeCheck.i] <- border.coefs[rangeCheck.i]
	lLike <- -negLLike(parVectMinCheck.i)
	if(is.null(lLike)){
		lLike<-NA
	}
	borderLogLikeError <- lLike-lpdensEps
	if(abs(lLike-lpdensEps) >= rangeTol*10){
		cat(sprintf('\r%4i %100s %+12.6f\n',rangeCheck.i,names(parVect[rangeCheck.i]),lLike-lpdensEps))
	}
	return(borderLogLikeError)
}

# requires baseNegLL in the global env
orderOfMagNegLLErrorFun <- function(delta,par.i){
	jParVect.i <- jParVect
	jParVect.i[par.i] <- jParVect[par.i] + delta
	return(abs(baseNegLL-jnegLLikelihood.f(jParVect.i))-1)
}
funFindParScale <- function(par.i,niter=100,useOrdersOfMagGuesses=F){
	if(!useOrdersOfMagGuesses|length(ordersOfMagGuesses)<par.i){
		# the fallback sweep, bounded per parameter
		minOrderOfMag <- ordersOfMagLimits[par.i,'min']
		maxOrderOfMag <- ordersOfMagLimits[par.i,'max']
	} else {
		minOrderOfMag <- ordersOfMagGuesses[par.i] -2
		maxOrderOfMag <- ordersOfMagGuesses[par.i] +1
	}
	# A parameter with no order of magnitude to work from, a zero width range or a
	# zero variance, has no scale to find.
	if(!is.finite(minOrderOfMag)||!is.finite(maxOrderOfMag)||
		 maxOrderOfMag<minOrderOfMag){
		cat(sprintf('\r%4i %-50s ... %+e                     \n',
								par.i,substr(names(jParVect)[par.i],1,50),NA))
		return(NA)
	}
	ordersOfMag <- minOrderOfMag:maxOrderOfMag
	cat(sprintf('%4i %-50s ... magscale:     ',
							par.i,substr(names(jParVect)[par.i],1,50)))
	ordersOfMagDeltRes <- c()
	ordersOfMagNegLLResp <- c()
	# secant starts every sweep below from x0=0, and the value there is known
	# without asking the model: orderOfMagNegLLErrorFun(0,par.i) puts jParVect back
	# exactly as it was, so jnegLLikelihood.f returns baseNegLL and the whole
	# expression is abs(baseNegLL-baseNegLL)-1.
	negLLErrorAtZero <- -1
	for(ord.i in 1:length(ordersOfMag)){
		cat(sprintf('\b\b\b\b\b\b\b\b\b\b%+10.1e',10^ordersOfMag[ord.i]))
		secantRes <- secant(orderOfMagNegLLErrorFun,x0=0,
																				x1=10^ordersOfMag[ord.i],
																				f0=negLLErrorAtZero,
																				niter=niter,
																				doWarn=F,tol=1e-2,par.i=par.i,
																				hasToBePositive=T)
		ordersOfMagDeltRes[ord.i] <- secantRes
		# secant reports the value at the point it hands back; the returns that never
		# evaluated that point leave it unset, and those still have to be measured.
		respFromSecant <- attr(secantRes,'fval')
		ordersOfMagNegLLResp[ord.i] <- if(is.null(respFromSecant)){
			orderOfMagNegLLErrorFun(ordersOfMagDeltRes[ord.i],par.i)
		} else {
			respFromSecant
		}
		if(!is.na(ordersOfMagNegLLResp[ord.i])&&!is.nan(ordersOfMagNegLLResp[ord.i])&&
			 abs(ordersOfMagNegLLResp[ord.i])<responseTolerance){
			break
		}
	}
	if(sum(is.na(ordersOfMagDeltRes))==length(ordersOfMagDeltRes) ||
		 length(ordersOfMagDeltRes)==0){
		cat(sprintf('\r%4i %-50s ... %+e                     \n',
								par.i,substr(names(jParVect)[par.i],1,50),NA))
		return(NA)
	} else {
		retScale <- ordersOfMagDeltRes[which.min(abs(ordersOfMagNegLLResp))]
		cat(sprintf('\r%4i %-50s ... %+e                     \n',
								par.i,substr(names(jParVect)[par.i],1,50),retScale))
		return(ordersOfMagDeltRes[which.min(abs(ordersOfMagNegLLResp))])
	}
}

# verify cached parscales
funParScaleHolds <- function(par.i,parscale){
	response <- orderOfMagNegLLErrorFun(parscale[[par.i]],par.i)
	return(is.finite(response)&&abs(response)<responseTolerance)
}

funOrderOfMagnitude <- function(x){
	return(floor(log10(abs(x))))
}

# ranges and cached determinations ####
# The Min and Max their authors gave each parameter, in the order
# sampleParms has them. Matching by name so a parameter that is not in
# frida_info comes out as NA rather than shifting following rows.
funParBoundsForSampleParms <- function(sampleParms,frida_info){
	parBounds <- frida_info[match(sampleParms$Variable,frida_info$Variable),c('Min','Max')]
	rownames(parBounds) <- sampleParms$Variable
	colnames(parBounds) <- c('Min','Max')
	return(parBounds)
}

# The parscales a previous run determined, as a list of the values, a logical
# saying which of them are remembered failures, and whether that run sampled the
# same parameters as this one.
#
# Entries are matched by name. The cached vector need not hold the same
# parameters, nor hold them in the same order.
funReadCachedParscale <- function(file,jParVectNames,redoFailedParscales=FALSE,
																	currentKey=NULL){
	parscale <- rep(NA_real_,length(jParVectNames))
	names(parscale) <- jParVectNames
	notDetermined <- rep(FALSE,length(jParVectNames))
	names(notDetermined) <- jParVectNames
	if(!file.exists(file)){
		return(list(parscale=parscale,notDetermined=notDetermined,
								parNamesChanged=FALSE))
	}
	cached <- readRDS(file)
	if(is.list(cached)&&!is.null(currentKey)){
		# Anything but the parameter list having moved means these values were
		# computed from something else. The parameter list alone leaves the values
		# both runs share, matched by name below and checked by the caller.
		mismatch <- funDeterminationKeyMismatch(cached$key,currentKey)
		if(length(mismatch)>0){
			cat(sprintf('Cached parscales in %s cannot be used, %s has changed. Redetermining.\n',
									basename(file),paste(mismatch,collapse='; ')))
			return(list(parscale=parscale,notDetermined=notDetermined,
									parNamesChanged=FALSE))
		}
	}
	if(is.list(cached)){
		parscale.old <- cached$parscale
		status.old <- cached$status
	} else {
		parscale.old <- cached
		status.old <- rep(NA_character_,length(parscale.old))
		names(status.old) <- names(parscale.old)
		cat(sprintf('%s predates the status record, failed determinations in it will be redone\n',
								basename(file)))
	}
	matches <- names(parscale)[names(parscale)%in%names(parscale.old)]
	parscale[matches] <- parscale.old[matches]
	if(!redoFailedParscales){
		notDetermined[matches] <- status.old[matches]%in%'notDetermined'
	}
	return(list(parscale=parscale,notDetermined=notDetermined,
							parNamesChanged=!setequal(names(parscale.old),jParVectNames)))
}

# What a cached determination was computed from. A changed model, changed
# calibration data or changed likelihood settings must invalidate the cached results.
#
# The model is the checksum of the frida directory without the files in exclude,
# the ones the analysis writes into it: the scenario and the export spec. A
# scenario that changes the likelihood shows in baseNegLL.
funDeterminationKey <- function(location.frida,location.frida.info,name.frida_info,
																calDat,resSigma,parNames,settings,
																baseNegLL=NULL,exclude=c()){
	md5Of <- function(path){
		if(!file.exists(path)){
			return('missing')
		}
		return(unname(tools::md5sum(path)))
	}
	digestOf <- function(x){
		if(is.null(x)){
			return('absent')
		}
		return(digest::digest(x))
	}
	list(fridaModel=funFridaFilesChecksum(location.frida,exclude=exclude)$checksum,
			 fridaInfo=md5Of(file.path(location.frida.info,name.frida_info)),
			 calDat=digestOf(calDat),
			 resSigma=digestOf(resSigma),
			 parNames=parNames,
			 settings=settings,
			 baseNegLL=baseNegLL,
			 writtenAt=Sys.time())
}

# Which parts of the key have moved since the cache was written, named so the log
# says why a determination is being redone. The parameter list is reported
# separately from everything else because it is the one mismatch that does not
# invalidate the parameters both runs have in common.
funDeterminationKeyMismatch <- function(cached,current){
	if(is.null(cached)||!is.list(cached)){
		return('no key was recorded with it')
	}
	mismatch <- character(0)
	described <- c(fridaModel='the FRIDA model files',
								 fridaInfo='frida_info',
								 calDat='the calibration data',
								 resSigma='the residual covariance')
	for(field in names(described)){
		if(!identical(cached[[field]],current[[field]])){
			mismatch <- c(mismatch,described[[field]])
		}
	}
	if(!is.null(cached$baseNegLL)&&!is.null(current$baseNegLL)&&
		 !isTRUE(all.equal(cached$baseNegLL,current$baseNegLL))){
		mismatch <- c(mismatch,
									sprintf('the likelihood at the starting parameters (%.10g, was %.10g)',
													current$baseNegLL,cached$baseNegLL))
	}
	for(setting in names(current$settings)){
		if(!identical(cached$settings[[setting]],current$settings[[setting]])){
			mismatch <- c(mismatch,sprintf('%s (%s, was %s)',setting,
																		 format(current$settings[[setting]]),
																		 if(is.null(cached$settings[[setting]])){
																		 	'not recorded'
																		 } else {
																		 	format(cached$settings[[setting]])
																		 }))
		}
	}
	return(mismatch)
}

# TRUE when the two runs do not sample the same parameters. Kept apart from the
# mismatch above: everything there invalidates the whole determination, while this
# only means some parameters are new.
funDeterminationParNameMismatch <- function(cached,current){
	if(is.null(cached)||!is.list(cached)){
		return(TRUE)
	}
	return(!identical(cached$parNames,current$parNames))
}

# The sampleParms a previous run left behind after determining ranges, or NULL
# when that file cannot stand in for a determination. Everything downstream of
# the determination rebuilds itself from these columns, so a file lacking any of
# them is redetermined rather than half used. A column that is there but holds
# nothing but NA is missing too.
funReadCachedRangedSampleParms <- function(file,currentKey=NULL,
																					 keyFile=paste0(tools::file_path_sans_ext(file),
																					 							 '.key.RDS')){
	required <- c('Variable','Value','Min','Max','MinAfterDet','MaxAfterDet',
								'MinNotDeterminedBorder','MaxNotDeterminedBorder',
								'parscale','parscaleStatus')
	# The key lives beside the file rather than in it, because this one is written
	# as a plain sampleParms and read back as one in several places.
	if(!is.null(currentKey)){
		cachedKey <- if(file.exists(keyFile)){readRDS(keyFile)}else{NULL}
		mismatch <- c(funDeterminationKeyMismatch(cachedKey,currentKey),
									if(funDeterminationParNameMismatch(cachedKey,currentKey)){
										'the set of sampled parameters'
									})
		if(length(mismatch)>0){
			cat(sprintf('Cached ranges in %s cannot be used, %s has changed. Redetermining.\n',
									basename(file),paste(mismatch,collapse='; ')))
			return(NULL)
		}
	}
	sampleParms <- readRDS(file)
	absent <- required[!required%in%colnames(sampleParms)]
	empty <- character(0)
	if(length(absent)==0){
		empty <- required[sapply(sampleParms[,required],function(x){all(is.na(x))})]
	}
	if(length(c(absent,empty))>0){
		cat(sprintf('Cached ranges in %s cannot be used, %s. Redetermining.\n',
								basename(file),
								paste(c(if(length(absent)>0){sprintf('no %s column',paste(absent,collapse=', '))},
												if(length(empty)>0){sprintf('nothing but NA in %s',paste(empty,collapse=', '))}),
											collapse='; ')))
		return(NULL)
	}
	return(sampleParms)
}

# sampling the parameter space ####
# Make the sampled range symmetric around the parameter value, and say what that
# did. Two kinds of range are left alone unless asked for:
#
#   external ranges, from frida_external_ranges.csv, are a deliberate statement
#   of the range to sample and are used as given, and
#
#   ranges that fell back to the author bound because a border could not be
#   determined, or was skipped for want of a parscale. Symmetrifying these is
#   what collapses ranges to zero width: the fallback bound is often exactly the
#   parameter value, so the smaller half width is zero and both sides snap onto
#   the value.
#
# The fallback test is per parameter, not per direction. Both borders determined,
# it symmetrifies as normal; neither, symmetrifyFallbackAuthorRanges decides;
# one of each, the two would disagree, so nothing is done and each side keeps
# what it was given. With the flag set it symmetrifies regardless, taking a
# fallback bound as equally valid as a determined one.
funSymmetrifyRanges <- function(sampleParms,parBounds,notDeterminedBorders,
																externalRangeParmNames=character(0),
																symmetricRanges='Min',
																allowAssymetricToAvoidZeroRanges=FALSE,
																symmetricRangesBoundByAuthors=TRUE,
																symmetrifyExternalRanges=FALSE,
																symmetrifyFallbackAuthorRanges=FALSE){
	if(!symmetricRanges%in%c('Min','Max')){
		cat(sprintf('Not symmetrifying parameter ranges (symmetricRanges is \'%s\')\n',
								as.character(symmetricRanges)[1]))
		return(sampleParms)
	}
	cat(sprintf('Symmetrifying ranges using procedure %s\n',symmetricRanges))
	widthBefore <- sampleParms$Max-sampleParms$Min
	# which parameters keep the range they came in with
	excludedExternal <- !symmetrifyExternalRanges &
		sampleParms$Variable%in%externalRangeParmNames
	excludedFallback <- !symmetrifyFallbackAuthorRanges &
		(notDeterminedBorders[,'Min']|notDeterminedBorders[,'Max'])
	excluded <- excludedExternal|excludedFallback
	if(symmetricRanges=='Max'){
		sampleParms$distance <- pmax(sampleParms$Value-sampleParms$Min,
																 sampleParms$Max-sampleParms$Value)
		if(symmetricRangesBoundByAuthors){
			sampleParms$distance <- pmin(sampleParms$distance,
																	 pmin(sampleParms$Value-parBounds[,'Min'],
																	 		 parBounds[,'Max']-sampleParms$Value))
		}
	} else {
		sampleParms$distance <- pmin(sampleParms$Value-sampleParms$Min,
																 sampleParms$Max-sampleParms$Value)
	}
	toApply <- !excluded
	# those that would have a distance of zero, we do not reassign
	if(allowAssymetricToAvoidZeroRanges){
		toApply <- toApply & sampleParms$distance!=0
	}
	toApply[is.na(toApply)] <- FALSE
	sampleParms$Max[toApply] <- sampleParms$Value[toApply]+sampleParms$distance[toApply]
	sampleParms$Min[toApply] <- sampleParms$Value[toApply]-sampleParms$distance[toApply]
	# report. The counts below partition sampleParms, so they add up to its rows
	widthAfter <- sampleParms$Max-sampleParms$Min
	tally <- function(x){sum(x&!excluded,na.rm=TRUE)}
	collapsed <- tally(widthAfter==0&widthBefore>0)
	narrowed <- tally(widthAfter<widthBefore&widthAfter>0)
	widened <- tally(widthAfter>widthBefore)
	unchanged <- tally(widthAfter==widthBefore)
	if(narrowed>0){cat(sprintf('  %5i parameter ranges narrowed\n',narrowed))}
	if(widened>0){cat(sprintf('  %5i parameter ranges widened\n',widened))}
	cat(sprintf('  %5i parameter ranges collapsed\n',collapsed))
	if(sum(excluded)>0){
		cat(sprintf('  %5i parameter ranges not symmetrified (%i author range fallback, %i external range, %i both)\n',
								sum(excluded),
								sum(excludedFallback&!excludedExternal),
								sum(excludedExternal&!excludedFallback),
								sum(excludedFallback&excludedExternal)))
	}
	if(unchanged>0){cat(sprintf('  %5i parameter ranges unchanged\n',unchanged))}
	return(sampleParms)
}

generateSobolSequenceForSampleParms <- function(sampleParms,numSample,
																								restretchSamplePoints=F,
																								ignoreExistingResults=F,
																								integerParms=NULL,
																								nullProb=0){
	if(!ignoreExistingResults && file.exists(file.path(location.output,'samplePoints.RDS'))){
		cat('Reading sampling points...')
		samplePoints <- readRDS(file.path(location.output,'samplePoints.RDS'))
		cat('done\n')
	} else {
		if(nrow(sampleParms)==1){
			samplePoints.base <- array(seq(0,1,length.out=numSample),dim=c(numSample,1))
			colnames(samplePoints.base) <- sampleParms[1,1]
			rownames(samplePoints.base) <- 1:numSample
		} else {
			cat('Generate sampling points using sobol sequence...')
			# sobolSequence.points generates points on the unit interval for each var
			# transformed, so vars are in rows samples in cols, makes the next steps easier
			samplePoints.base <- sobolSequence.points(nrow(sampleParms),31,numSample) 
			if(sum(duplicated(samplePoints.base))>0){
				stop('Not enough unique sample points. Check the sobol generation\n')
			}
		}
		if(nullProb>0){
			samplePoints.base[samplePoints.base<=nullProb] <- NA
			samplePoints.base <- samplePoints.base[!duplicated(samplePoints.base),]
			samplePoints.base[!is.na(samplePoints.base)] <- 
				(samplePoints.base[!is.na(samplePoints.base)]-nullProb)/(1-nullProb)
		}
		if(!is.null(integerParms)){
			sampleParms[sampleParms$Variable %in% integerParms$Variable,c('Max')] <- 
				sampleParms[sampleParms$Variable %in% integerParms$Variable,c('Max')] + 1
		}
		samplePoints <- funStretchSamplePoints(samplePoints.base,sampleParms,
																					 restretchSamplePoints)
		colnames(samplePoints) <- sampleParms$Variable
		if(!is.null(integerParms) && nrow(integerParms)>0){
			cat('rounding integer parms...')
			for(p.i in 1:nrow(integerParms)){
				if(integerParms$Variable[p.i]%in%sampleParms$Variable){
					samplePoints[,as.character(integerParms$Variable[p.i])] <-
						 floor(samplePoints[,as.character(integerParms$Variable[p.i])])
					samplePoints[samplePoints[,as.character(integerParms$Variable[p.i])]==sampleParms$Max[p.i],
											 as.character(integerParms$Variable[p.i])] <- sampleParms$Max[p.i]-1
				}
			}
			samplePoints <- samplePoints[!duplicated(samplePoints),]
		}
		if(nrow(sampleParms)==1){
			samplePoints <- array(samplePoints,dim=c(numSample,1))
			colnames(samplePoints) <- sampleParms[1,1]
			rownames(samplePoints) <- 1:numSample
		}
		if(sum(duplicated(samplePoints))>0){
			stop('Not enough possible combinations in specified parameters to satisfy numSample\n')
		}
		saveRDS(samplePoints,file.path(location.output,'samplePoints.RDS'))
		cat('done\n')
	}
	return(samplePoints)
}


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
