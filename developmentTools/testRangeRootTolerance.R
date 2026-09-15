# testRangeRootTolerance.R ####
#
# What rangeRootTol does in the border search, against a stand-in objective. It
# is the one knob that can move a border.
#
#   Where the likelihood is smooth near the border, the tolerance changes
#   nothing. Brent's method converges in about sixteen evaluations whatever
#   tolerance it is given, and never comes near the iteration limit: the
#   bisection fallback terminates once the bracket stops shrinking.
#
#   Where the border falls on a discontinuity it costs. A parameter pushed far
#   enough that the model stops completing gets logLike.failedRun, a cliff rather
#   than a crossing, and there a tight tolerance keeps subdividing an interval
#   whose root is not a root. That is the case the roughly hundred not-determined
#   borders per direction are in, and the case rangeRootTol is for.
#
# Run from the repository root:
#   Rscript developmentTools/testRangeRootTolerance.R

source('funParmSpace.R')

mu <- 3
sdev <- 2
lpdensEps <- -2
outputDigits <- 6   # stella writes its output at finite precision
analytic <- mu+sdev*sqrt(-2*lpdensEps)

nCalls <- 0
# cliff stands in for the point past which the model no longer completes, whose
# runs come back with the failed-run marker instead of a likelihood
makeNegLLike <- function(cliff){
	function(parVect,...){
		nCalls <<- nCalls+1
		x <- as.numeric(parVect[1])
		if(x>cliff){
			return(1e7)
		}
		signif(0.5*((x-mu)/sdev)^2,outputDigits)
	}
}

runSearch <- function(rootTolFactor,rootMaxIter,parscale=sdev){
	nCalls <<- 0
	r <- suppressWarnings(findDensValBorder(
		1,parVect=c(p1=mu),lpdensEps=lpdensEps,
		ceterisParibusPars=TRUE,idcToMod=1,parscale=parscale,
		bounds=cbind(mu-40*sdev,mu+40*sdev),
		tol=1e-15,niter=1000,max=TRUE,
		rootTolFactor=rootTolFactor,rootMaxIter=rootMaxIter))
	list(border=as.numeric(r),calls=nCalls)
}

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),'')))
		fail <<- fail+1}
}

cat(sprintf('contour to find: %.10g\n\n',analytic))
cat(sprintf('%-34s %12s %11s %12s %8s\n',
						'','border','model runs','error','saving'))

# ---- smooth: the border is well clear of the cliff
negLLike <- makeNegLLike(Inf)
old.smooth <- runSearch(NA,1e3)
new.smooth <- runSearch(1e-4,60)
cat(sprintf('%-34s %12.8g %11d %12.2e\n','smooth, tol=1e-16 maxiter=1000',
						old.smooth$border,old.smooth$calls,abs(old.smooth$border-analytic)))
cat(sprintf('%-34s %12.8g %11d %12.2e %7.2fx\n','smooth, tol=parscale*1e-4 mi=60',
						new.smooth$border,new.smooth$calls,abs(new.smooth$border-analytic),
						old.smooth$calls/max(new.smooth$calls,1)))

# ---- the model stops completing before the contour is reached
negLLike <- makeNegLLike(6.5)
old.cliff <- runSearch(NA,1e3)
new.cliff <- runSearch(1e-4,60)
cat(sprintf('%-34s %12.8g %11d %12s\n','cliff,  tol=1e-16 maxiter=1000',
						old.cliff$border,old.cliff$calls,'n/a'))
cat(sprintf('%-34s %12.8g %11d %12s %7.2fx\n','cliff,  tol=parscale*1e-4 mi=60',
						new.cliff$border,new.cliff$calls,'n/a',
						old.cliff$calls/max(new.cliff$calls,1)))

cat('\n')
check('a smooth border costs the same either way, no thrashing to maxiter',
			new.smooth$calls==old.smooth$calls&&old.smooth$calls<50,
			sprintf('old %d, new %d',old.smooth$calls,new.smooth$calls))
check('a smooth border stays within one part in 1e3 of the contour',
			abs(new.smooth$border-analytic)/abs(analytic-mu) < 1e-3,
			sprintf('relative error %.3e',abs(new.smooth$border-analytic)/abs(analytic-mu)))
check('a border on a discontinuity costs materially less',
			old.cliff$calls/max(new.cliff$calls,1) > 2,
			sprintf('old %d, new %d',old.cliff$calls,new.cliff$calls))
check('and lands in the same place',
			abs(new.cliff$border-old.cliff$border)/abs(old.cliff$border) < 1e-3,
			sprintf('old %.10g, new %.10g',old.cliff$border,new.cliff$border))
check('NA restores the old absolute tolerance exactly',
			identical(runSearch(NA,1e3)$border,old.cliff$border))
check('a parameter with no usable parscale falls back to the old tolerance',
			identical(runSearch(1e-4,60,parscale=0)$border,runSearch(NA,60)$border))

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('the border root tolerance does not behave as intended\n')
}
cat('\nMeasured on a quantised gaussian, not on FRIDA. What this cannot say is how\n')
cat('many real borders sit on a discontinuity; only a counted run can.\n')
