# testFullReoptBranch.R ####
#
# The branch of findDensValBorder that reoptimises the other parameters at each
# step (ceterisParibusPars = FALSE) is a supported configuration, but nothing
# exercises it: the production config sets treatVarsAsIndep, which makes
# ceterisParibusPars TRUE and returns before any of it runs. Four defects had
# accumulated in it by the time anyone looked.
#
# This is the standing guard against that happening again. findDensValBorder takes
# its objective from negLLike in the global environment, which is what makes the
# fixture possible: a synthetic gaussian likelihood stands in for the stella run,
# so the branch is exercised in seconds on a handful of parameters with no model,
# no cluster and no calibration data.
#
# Run from the repository root:  Rscript developmentTools/testFullReoptBranch.R

suppressPackageStartupMessages(library(optimx))
source('funParmSpace.R')

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){
		cat(sprintf('  ok   %s\n',label)); ok <<- ok+1
	} else {
		cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),'')))
		fail <<- fail+1
	}
}

# ---- the stand-in model ####
# A gaussian log density, so every border this fixture asks for has a closed form.
nPar <- 6
mu <- c(1,2,3,4,5,6)
sd <- c(1,1,2,1,0.5,1)
names(mu) <- names(sd) <- paste0('p',1:nPar)
corMat <- diag(nPar)

negLLikeCalls <- NULL
recordCalls <- FALSE
makeNegLLike <- function(corMat){
	sigma <- diag(sd)%*%corMat%*%diag(sd)
	sigmaInv <- solve(sigma)
	function(parVect,...){
		if(recordCalls){
			negLLikeCalls[[length(negLLikeCalls)+1]] <<- parVect
		}
		d <- as.numeric(parVect)-mu
		0.5*as.numeric(t(d)%*%sigmaInv%*%d)
	}
}
negLLike <- makeNegLLike(corMat)

parVect <- mu
# the contour to find, in log density units below the maximum
lpdensEps <- -2
parscale <- sd
bounds <- cbind(mu-20*sd,mu+20*sd)

# ---- 15a: the ranged parameter must not be handed to the inner optimiser ####
# densMaxGivenParFun is where that happens. idcToMod is deliberately not 1:n and
# not 1 based, which is the case the old positional negative index got wrong: for
# idcToMod = c(2,5,8) and parIdx = 5, idcToMod[-5] is out of range and R hands
# back the vector unchanged, so parIdx stays in the set the optimiser moves.
cat('15a  the ranged parameter is held fixed\n')
idcNonContig <- c(2,5,8)
parIdx.a <- 5
parVect.a <- setNames(rep(0,8),paste0('q',1:8))
negLLikeCalls <- list()
recordCalls <- TRUE
negLLike.orig <- negLLike
negLLike <- function(parVect,...){
	negLLikeCalls[[length(negLLikeCalls)+1]] <<- parVect
	0
}
invisible(densMaxGivenParFun(otherPars=c(11,88),parVect=parVect.a,
														 parIdx=parIdx.a,idcToMod=idcNonContig))
negLLike <- negLLike.orig
recordCalls <- FALSE
built <- negLLikeCalls[[1]]
check('the ranged parameter keeps its incoming value',
			built[parIdx.a]==parVect.a[parIdx.a],
			sprintf('parVect[%d] is %s, was %s',parIdx.a,
							format(built[parIdx.a]),format(parVect.a[parIdx.a])))
check('the other parameters in idcToMod take the optimiser values',
			all(built[setdiff(idcNonContig,parIdx.a)]==c(11,88)))
check('parameters outside idcToMod are untouched',
			all(built[setdiff(1:8,idcNonContig)]==0))

# ---- 15d: a single element idcToMod must not index past the end ####
cat('15d  a one element idcToMod is handled\n')
res <- try(findDensValBorder(1,parVect=parVect,lpdensEps=lpdensEps,
														 ceterisParibusPars=FALSE,
														 idcToMod=list(1),
														 parscale=parscale,bounds=bounds,
														 tol=1e-4,niter=50,max=TRUE),silent=TRUE)
check('no error on a one element idcToMod',!inherits(res,'try-error'),
			if(inherits(res,'try-error')){as.character(res)}else{''})

# ---- 15b: tracing must not change the answer ####
# The corrective reoptimisation after a likelihood improvement used to sit inside
# if(trace>0), so a traced run and an untraced run computed different borders.
cat('15b  tracing does not change the border\n')
runBorder <- function(trace,parIdx=3,max=TRUE,corMat=diag(nPar)){
	negLLike <<- makeNegLLike(corMat)
	out <- capture.output(
		r <- suppressWarnings(findDensValBorder(parIdx,parVect=parVect,lpdensEps=lpdensEps,
																						ceterisParibusPars=FALSE,
																						idcToMod=1:nPar,
																						parscale=parscale,bounds=bounds,
																						tol=1e-4,niter=50,trace=trace,max=max)))
	as.numeric(r)
}
b.untraced <- runBorder(0)
b.traced <- runBorder(1)
check('untraced and traced borders agree',
			isTRUE(all.equal(b.untraced,b.traced,tolerance=1e-6)),
			sprintf('trace=0 gives %.10g, trace=1 gives %.10g',b.untraced,b.traced))

# ---- the branch finds the right contour ####
# With independent parameters the others reoptimise back to mu, so the profile
# border and the conditional border coincide and both have a closed form.
cat('     the border sits on the requested contour\n')
analytic <- mu[3]+sd[3]*sqrt(-2*lpdensEps)
check('independent case matches the closed form',
			isTRUE(all.equal(b.untraced,as.numeric(analytic),tolerance=1e-3)),
			sprintf('found %.10g, closed form %.10g',b.untraced,analytic))

# ---- the reoptimisation actually does something ####
# Under correlation the profile border, which lets the other parameters move, lies
# strictly outside the conditional border that holds them at mu. If the full
# branch ever silently degrades into the ceteris paribus one, this is what notices.
cat('     the reoptimisation widens the border under correlation\n')
corCorrelated <- diag(nPar)
corCorrelated[3,4] <- corCorrelated[4,3] <- 0.9
negLLike <- makeNegLLike(corCorrelated)
b.profile <- suppressWarnings(as.numeric(
	findDensValBorder(3,parVect=parVect,lpdensEps=lpdensEps,
										ceterisParibusPars=FALSE,idcToMod=1:nPar,
										parscale=parscale,bounds=bounds,
										tol=1e-4,niter=50,max=TRUE)))
b.conditional <- suppressWarnings(as.numeric(
	findDensValBorder(3,parVect=parVect,lpdensEps=lpdensEps,
										ceterisParibusPars=TRUE,idcToMod=1:nPar,
										parscale=parscale,bounds=bounds,
										tol=1e-4,niter=50,max=TRUE)))
check('profile border lies outside the conditional one',
			b.profile > b.conditional*(1+1e-6),
			sprintf('profile %.10g, conditional %.10g',b.profile,b.conditional))

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('the full reoptimisation branch of findDensValBorder is not behaving\n')
}
