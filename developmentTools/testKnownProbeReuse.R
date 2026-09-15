# testKnownProbeReuse.R ####
#
# Two probes in these searches are values the caller already has:
#
#   orderOfMagNegLLErrorFun(0,par.i) puts jParVect back exactly as it was, so
#   jnegLLikelihood.f returns baseNegLL and the expression is
#   abs(baseNegLL-baseNegLL)-1. Constant -1.
#
#   The border search's probe at the starting point is -negLLike(parVect) minus
#   lpdensEps, and the caller measured -negLLike(parVect) to derive lpdensEps in
#   the first place. Same value for every parameter and both directions.
#
# Substituting a value for a measurement is only safe if it is the same value, so
# that is what this checks, alongside the saving.
#
# Run from the repository root:
#   Rscript developmentTools/testKnownProbeReuse.R

suppressPackageStartupMessages(library(optimx))
source('funParmSpace.R')

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),'')))
		fail <<- fail+1}
}

# ---- the zero delta probe ####
cat('orderOfMagNegLLErrorFun at delta zero\n')
jParVect <- c(a=1.5,b=-2.25,c=1e6)
sampleParms <- data.frame(Variable=c('a','b'),Value=c(1.5,-2.25))
treatVarsAsIndep <- TRUE
calDat <- data.frame(x=1:3)
jnegLLikelihood.f <- function(v){sum(as.numeric(v)^2)+0.5}
baseNegLL <- jnegLLikelihood.f(jParVect)
measured <- sapply(seq_along(jParVect),function(i){orderOfMagNegLLErrorFun(0,i)})
check('is exactly -1 for every parameter',all(measured==-1),
			paste(format(measured),collapse=', '))
check('and funFindParScale now uses that instead of measuring it',
			any(grepl('negLLErrorAtZero <- -1',readLines('funParmSpace.R'),fixed=TRUE)))

# ---- the starting point probe ####
cat('the border search probe at the starting parameters\n')
mu <- 3
sdev <- 2
likeCutoffRatio <- 1000
nCalls <- 0
negLLike <- function(parVect,...){
	nCalls <<- nCalls+1
	0.5*((as.numeric(parVect[1])-mu)/sdev)^2
}
parVect <- c(p1=mu)
lpdensAtParVect <- -negLLike(parVect)
lpdensEps <- lpdensAtParVect-log(likeCutoffRatio)

check('the substituted value equals what a measurement would give',
			isTRUE(all.equal(lpdensAtParVect-lpdensEps,
											 likeGoalDiffFun(parVect[[1]],parVect,1,lpdensEps))),
			sprintf('substituted %.12g, measured %.12g',lpdensAtParVect-lpdensEps,
							likeGoalDiffFun(parVect[[1]],parVect,1,lpdensEps)))
check('and to log(likeCutoffRatio), which is what makes it constant',
			isTRUE(all.equal(lpdensAtParVect-lpdensEps,log(likeCutoffRatio))))

runSearch <- function(useKnown,max=TRUE){
	nCalls <<- 0
	r <- suppressWarnings(findDensValBorder(
		1,parVect=parVect,lpdensEps=lpdensEps,
		ceterisParibusPars=TRUE,idcToMod=1,parscale=sdev,
		bounds=cbind(mu-40*sdev,mu+40*sdev),
		tol=1e-15,niter=1000,max=max,
		lpdensAtParVect=if(useKnown){lpdensAtParVect}else{NULL}))
	list(border=as.numeric(r),calls=nCalls)
}

for(direction in c(TRUE,FALSE)){
	without <- runSearch(FALSE,max=direction)
	with <- runSearch(TRUE,max=direction)
	cat(sprintf('  %s: %d runs without, %d with, border %.12g vs %.12g\n',
							ifelse(direction,'max','min'),without$calls,with$calls,
							without$border,with$border))
	check(sprintf('%s border is unchanged',ifelse(direction,'max','min')),
				isTRUE(all.equal(without$border,with$border)),
				sprintf('%.15g vs %.15g',without$border,with$border))
	check(sprintf('%s costs one model run fewer',ifelse(direction,'max','min')),
				with$calls==without$calls-1,
				sprintf('%d vs %d',with$calls,without$calls))
}

# ---- the substitution must stop once the parameters move ####
# In the full reoptimisation branch parVect changes between iterations, and the
# known value stops being the value at that point.
cat('the substitution is dropped once parVect has moved\n')
negLLike <- function(parVect,...){
	nCalls <<- nCalls+1
	0.5*sum((as.numeric(parVect)-c(mu,mu))^2/sdev^2)
}
parVect2 <- c(p1=mu,p2=mu)
lpdensAtParVect2 <- -negLLike(parVect2)
lpdensEps2 <- lpdensAtParVect2-log(likeCutoffRatio)
b.known <- suppressWarnings(as.numeric(findDensValBorder(
	1,parVect=parVect2,lpdensEps=lpdensEps2,ceterisParibusPars=FALSE,
	idcToMod=1:2,parscale=c(sdev,sdev),bounds=cbind(c(mu-40*sdev,mu-40*sdev),
																									c(mu+40*sdev,mu+40*sdev)),
	tol=1e-4,niter=50,max=TRUE,lpdensAtParVect=lpdensAtParVect2)))
b.plain <- suppressWarnings(as.numeric(findDensValBorder(
	1,parVect=parVect2,lpdensEps=lpdensEps2,ceterisParibusPars=FALSE,
	idcToMod=1:2,parscale=c(sdev,sdev),bounds=cbind(c(mu-40*sdev,mu-40*sdev),
																									c(mu+40*sdev,mu+40*sdev)),
	tol=1e-4,niter=50,max=TRUE)))
check('the reoptimising branch gets the same border either way',
			isTRUE(all.equal(b.known,b.plain,tolerance=1e-9)),
			sprintf('%.12g vs %.12g',b.known,b.plain))

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('a substituted probe does not match the measurement it replaces\n')
}
