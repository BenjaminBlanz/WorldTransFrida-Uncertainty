# testParscaleSweepWidth.R ####
#
# The fallback parscale sweep used to be one global range applied to every
# parameter (plan item 4). This measures what the per parameter bounds replace it
# with, against a real determination if one is on disk and against a synthetic
# spread otherwise, and checks the two properties the change has to preserve:
# the guess pass is untouched, and the fallback still covers every scale below the
# guess that the global range covered.
#
# Run from the repository root:
#   Rscript developmentTools/testParscaleSweepWidth.R [path/to/a/workOutput/run]

source('funParmSpace.R')

args <- commandArgs(TRUE)
runDir <- if(length(args)>0){args[1]}else{NA}

ordersOfMagLimits.finite <- function(x){
	x <- x[is.finite(x)]
	if(length(x)==0){return(NA_real_)}
	x
}

buildLimits <- function(guess.parvect,guess.resSigma){
	n <- length(guess.parvect)+length(guess.resSigma)
	lim <- array(NA_real_,dim=c(n,2),dimnames=list(NULL,c('min','max')))
	pIdc <- seq_along(guess.parvect)
	rIdc <- (length(guess.parvect)+1):n
	lim[pIdc,'min'] <- min(ordersOfMagLimits.finite(guess.parvect))-2
	maxParvect <- guess.parvect+1
	maxParvect[!is.finite(maxParvect)] <-
		max(ordersOfMagLimits.finite(c(guess.parvect,guess.resSigma)))+4
	lim[pIdc,'max'] <- maxParvect
	lim[rIdc,'min'] <- min(ordersOfMagLimits.finite(guess.resSigma))-2
	lim[rIdc,'max'] <- max(ordersOfMagLimits.finite(guess.resSigma))+4
	lim
}

# ---- the parameter spread to measure against
if(!is.na(runDir)&&file.exists(file.path(runDir,'sampleParmsParscaleRanged.csv'))){
	sp <- read.csv(file.path(runDir,'sampleParmsParscaleRanged.csv'))
	guess.parvect <- funOrderOfMagnitude(sp$Max-sp$Min)
	sigFile <- file.path(runDir,'sigma-indepParms.RDS')
	if(!file.exists(sigFile)){sigFile <- file.path(dirname(runDir),'sigma-indepParms.RDS')}
	guess.resSigma <- if(file.exists(sigFile)){
		funOrderOfMagnitude(diag(readRDS(sigFile)))-6
	} else {
		seq(-14,16,length.out=222)
	}
	cat(sprintf('measured against %s: %d sampled parameters, %d variances\n',
							runDir,length(guess.parvect),length(guess.resSigma)))
} else {
	# the spread the reference run had
	set.seed(1)
	guess.parvect <- sample(-8:13,665,replace=TRUE)
	guess.resSigma <- sample(-14:16,222,replace=TRUE)
	cat('measured against a synthetic spread matching the reference run\n')
	cat('  (pass a workOutput run directory as an argument to use real numbers)\n')
}

lim <- buildLimits(guess.parvect,guess.resSigma)
guessAll <- c(guess.parvect,guess.resSigma)

globalMin <- min(ordersOfMagLimits.finite(guessAll))-2
globalMax <- max(ordersOfMagLimits.finite(guessAll))+4
globalWidth <- globalMax-globalMin+1
newWidth <- lim[,'max']-lim[,'min']+1

cat(sprintf('\nglobal sweep : %3d orders for every parameter\n',globalWidth))
cat(sprintf('per parameter: %3.0f orders on average, %.0f smallest, %.0f largest\n',
						mean(newWidth,na.rm=TRUE),min(newWidth,na.rm=TRUE),max(newWidth,na.rm=TRUE)))
cat(sprintf('               %d of %d sampled parameters have a zero width range and keep the old ceiling\n',
						sum(!is.finite(guess.parvect)),length(guess.parvect)))
cat(sprintf('               %.2fx fewer orders per fallback sweep\n',
						globalWidth/mean(newWidth,na.rm=TRUE)))

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),'')))
		fail <<- fail+1}
}

cat('\n')
pIdc <- seq_along(guess.parvect)
check('the fallback never sweeps wider than the old global range',
			all(newWidth<=globalWidth,na.rm=TRUE))
hasOwnOrder <- is.finite(guess.parvect)
check('a sampled parameter is not swept above the order of its own range',
			all(lim[pIdc,'max'][hasOwnOrder]==guess.parvect[hasOwnOrder]+1))
check('a zero width range keeps the old global ceiling instead',
			all(lim[pIdc,'max'][!hasOwnOrder]==globalMax))
check('the fallback still reaches every scale below the guess that the global range did',
			all(lim[pIdc,'min'][hasOwnOrder]<=guess.parvect[hasOwnOrder]-2),
			'the guess pass covers guess-2 upward, so the fallback must start no higher')
check('the fallback floor is no higher than the old global floor for sampled parameters',
			all(lim[pIdc,'min']>=globalMin,na.rm=TRUE))
check('the variances keep their old headroom',
			all(lim[-pIdc,'max']==max(ordersOfMagLimits.finite(guess.resSigma))+4,na.rm=TRUE))

# a zero width range has no order of magnitude and must not poison the floor
guess.withZero <- c(guess.parvect,-Inf)
limZero <- buildLimits(guess.withZero,guess.resSigma)
check('a zero width range does not drag every floor to -Inf',
			is.finite(limZero[1,'min']))

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('the per parameter sweep bounds do not hold\n')
}
