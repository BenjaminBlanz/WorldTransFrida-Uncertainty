# The following functions work together to 
# find the borders of the density where density is equal to pdensEps
# repeat two steps, find root in current parm
# maximize density for root of cuurrent parm using other parms
# This function relies on the negLLike function being present in the global env.
# This allows the user to specify the type of likelihood function.
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
	if(sum(colnames(calDat)==colnames(runDat))<ncol(calDat)){
		stop('not all calDat columns in runDat\n')
	}
	resDat <- calDat-runDat[1:nrow(calDat),colnames(calDat)]
	lLikelihood <- funLogLikelihood(resDat,resSigma)
	# If the logLike is not NA but the run did not complete assign 
	# lowest value. We use this when narrowing the parms space
	if(is.na(runDat[[1]][nrow(runDat)])){
		lLikelihood <- -.Machine$double.xmax+(sum(!is.na(runDat[[1]]))*.Machine$double.eps)
	}
	return(-lLikelihood)
}
negLLike <- function(parVect){
	runDat <- runFRIDASpecParms(parVect)
	if(sum(colnames(calDat)==colnames(runDat))<ncol(calDat)){
		stop('not all calDat columns in runDat\n')
	}
	resDat <- calDat-runDat[1:nrow(calDat),colnames(calDat)]
	lLikelihood <- funLogLikelihood(resDat,resSigma)
	# If the logLike is not NA but the run did not complete assign 
	# lowest value. We use this when narrowing the parms space
	if(is.na(runDat[[1]][nrow(runDat)])){
		lLikelihood <- -.Machine$double.xmax+(sum(!is.na(runDat[[1]]))*.Machine$double.eps)
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
densMaxGivenParFun <- function(otherPars,parVect,parIdx,idcToMod, ...){
	parVect[idcToMod[-parIdx]] <- otherPars
	return(negLLike(parVect, ...))
}
# if ceterisParibusPars is TRUE, the densValBorder is found for the selecteed par Idx,
# keeping all other pars at the values in parVect, i.e. this will not account
# for rotatet elipsoid parameter distributions, but just for the slice through
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
															workerStagger=FALSE,
															...){
	if(workerStagger){
		Sys.sleep(workerID*0.04)
	}
	if(length(parIdx)>1){
		stop('only one parIdx at a time\n')
	}
	if(is.list(idcToMod)){
		idcToMod <- idcToMod[[parIdx]]
	}
	idcToMod.base <- idcToMod
	for(idcsToMod.i in 2:length(idcToMod.base)){
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
			#if there is no sign change between the endpoints of root.range, use secant's
			#method otherwise use uniroot
			if(((likeGoalDiffFun(root.range[1],parVect,parIdx,lpdensEps,...)>0)-
					(likeGoalDiffFun(root.range[2],parVect,parIdx,lpdensEps,...)>0))==0){
				if(max){
					root.range <- c(par.val,par.val+parscale[parIdx])
				} else {
					root.range <- c(par.val-parscale[parIdx],par.val)
				}
				par.val.old <-par.val
				par.val <- secant(likeGoalDiffFun,
													root.range[1],root.range[1]*1.001,
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
																						parVect=parVect,
																						parIdx=parIdx,
																						lpdensEps=lpdensEps,
																						tol = 1e-16, maxiter = niter,...)$root)
			}
			if(ceterisParibusPars){
				return(par.val)
			} else {
				parVect[parIdx] <- par.val
				if(trace>0){
					cat('iter ',iter,' ',par.val,' : ')
				}
				# maximize density at root using other parms
				otherIdx <- idcToMod[-parIdx]
				otherPars <- parVect[otherIdx]
				res <- suppressWarnings(optimx(otherPars,densMaxGivenParFun,
																			 method = 'Nelder-Mead',
																			 control=list(dowarn = F,
																			 						 parscale=parscale[otherIdx]),
																			 parVect=parVect,
																			 parIdx=parIdx, 
																			 idcToMod=idcToMod,...))
				parVect[otherIdx] <- unlist(res[1:length(otherIdx)])
				likeAtMaxOld <- likeAtMax
				likeAtMax <- -res$value
				if(trace>0){
					cat(parVect,' ',likeAtMax-lpdensEps,'\n')
				}
				if(likeAtMax>likeAtMaxOld){
					if(trace>0){
						cat('Likelihood Imporovement After Step Reoptimizing\n')
						res <- suppressWarnings(optimx(parVect,negLLike,
																					 method = 'Nelder-Mead',
																					 control=list(dowarn = F,
																					 						 parscale=parscale,...)))#,trace=99)))
						parVect <- unlist(res[1:length(parVect)])
						par.val <- parVect[parIdx]
						likeAtMax <- -res$value
						cat(parVect,' ',likeAtMax-lpdensEps,'\n')
					}
				}
			}
			iter <- iter+1
		} 
	}
	return(par.val)
}


secant <- function(fun, x0, x1, tol=1e-07, niter=1e4, doWarn=T, trace=0,
									 bound=NULL,hasToBePositive=FALSE,...){
	if(is.null(bound)){
		bound <- sign(x1-x0)*Inf
	}
	for ( i in 1:niter ){
		# cat(sprintf('x0=%10.2e x1=%10.2e',x0,x1))
		f0 <- fun(x0,...)
		f1 <- fun(x1,...)
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
		if(abs(fun(x2,...)) < tol || abs(x2)>abs(bound)){
			return(x2)
		}
		if(x0==x2){
			if(doWarn){
				warning("In secant cycle detected\n")
			}
			return(x2)
		}
		# cat(sprintf(' x2=%10.2e\n',x2))
		x0 <- x1
		x1 <- x2
	}
	if(doWarn){
		warning("In secant exceeded allowed number of iterations\n")
	}
	return(x2)
}


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
		minOrderOfMag <- min(ordersOfMagLimits)
		maxOrderOfMag <- max(ordersOfMagLimits)
	} else {
		minOrderOfMag <- ordersOfMagGuesses[par.i] -2
		maxOrderOfMag <- ordersOfMagGuesses[par.i] +1
	}
	ordersOfMag <- minOrderOfMag:maxOrderOfMag
	cat(sprintf('%4i %-50s ... magscale:     ',
							par.i,substr(names(jParVect)[par.i],1,50)))
	ordersOfMagDeltRes <- c()
	ordersOfMagNegLLResp <- c()
	for(ord.i in 1:length(ordersOfMag)){
		cat(sprintf('\b\b\b\b\b\b\b\b\b\b%+10.1e',10^ordersOfMag[ord.i]))
		ordersOfMagDeltRes[ord.i] <- secant(orderOfMagNegLLErrorFun,x0=0,
																				x1=10^ordersOfMag[ord.i],
																				niter=niter,
																				doWarn=F,tol=1e-2,par.i=par.i,
																				hasToBePositive=T)
		ordersOfMagNegLLResp[ord.i] <- orderOfMagNegLLErrorFun(ordersOfMagDeltRes[ord.i],par.i)
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

funOrderOfMagnitude <- function(x){
	return(floor(log10(abs(x))))
}

# sobol sequence ####
generateSobolSequenceForSampleParms <- function(sampleParms,numSample,
																								restretchSamplePoints=F,
																								ignoreExistingResults=F,
																								integerParms=NULL){
	if(!ignoreExistingResults && file.exists(file.path(location.output,'samplePoints.RDS'))){
		cat('Reading sampling points...')
		samplePoints <- readRDS(file.path(location.output,'samplePoints.RDS'))
		samplePoints.base <- readRDS(file.path(location.output,'samplePointsBase.RDS'))
		cat('done\n')
	} else {
		cat('Generate sampling points using sobol sequence...')
		# sobolSequence.points generates points on the unit interval for each var
		# transformed, so vars are in rows samples in cols, makes the next steps easier
		samplePoints.base <- sobolSequence.points(nrow(sampleParms),31,numSample) 
		if(sum(duplicated(samplePoints.base))>0){
			stop('Not enough unique sample points. Check the sobol generation\n')
		}
		samplePoints <- funStretchSamplePoints(samplePoints.base,sampleParms,restretchSamplePoints)
		# samplePoints <- rbind(samplePoints, t(sampleParms$Value))
		if(!is.null(integerParms)){
			cat('rounding integer parms...')
			for(p.i in 1:nrow(integerParms)){
				if(integerParms$Variable[p.i]%in%sampleParms$Variable){
					samplePoints[,integerParms$Variable[p.i]] <- round(samplePoints[,integerParms$Variable[p.i]])
				}
			}
		}
		# if('Climate Units.selected climate case'%in%sampleParms$Variable){
		# 	samplePoints[,'Climate Units.selected climate case'] <- round(samplePoints[,'Climate Units.selected climate case'])
		# }
		saveRDS(samplePoints,file.path(location.output,'samplePoints.RDS'))
		saveRDS(samplePoints.base,file.path(location.output,'samplePointsBase.RDS'))
		cat('done\n')
	}
	return(samplePoints)
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
