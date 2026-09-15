suppressPackageStartupMessages({
	library(Rmpfr,quietly=T,warn.conflicts = F) # arbitrary precision math used to calculate the likelihood from loglikelihood
	library(optimx,quietly=T,warn.conflicts = F) # interface to various optimizers
	library(tictoc,quietly=T,warn.conflicts = F) # simple timing measurements
	library(SobolSequence,quietly=T,warn.conflicts = F) # generates multidimensional sobol sequences
	library(lubridate,quietly=T,warn.conflicts = F) # deals with times
	library(cNORM,quietly=T,warn.conflicts = F) # wighted quantiles
	library(spatstat.explore,quietly=T,warn.conflicts = F) # used for the quantile.density function
	library(caret,quietly=T,warn.conflicts = F) # to find linear combinations and remove them in the calib dat
	library(matrixcalc,quietly=T,warn.conflicts = F) # to test positive definitnes of cov matrix
	#library(imputeTS,quietly=T,warn.conflicts = F) # used for interpolating missing values # Only needed for interpolation after MLE, fails for R v4.4 on Levante
	library(data.table,quietly=T,warn.conflicts = F) # fast csv reading/writing of the per var files
	library(parallel)
})
# we parallelise over processes (workers/variables), so data.table must not
# additionally spawn OpenMP threads, that would oversubscribe the node
data.table::setDTthreads(1)

# logLike.failedRun ####
# Marker for a run that failed or produced no usable output; anything greater
# counts as a complete run. It has to survive the 15 significant digit csv round
# trip of the per variable files, which -.Machine$double.xmax does not, so the
# loop finds the largest magnitude negative fixed point of that round trip.
logLike.failedRun <- (function(){
	x <- .Machine$double.xmax
	for(i in 1:64){
		y <- as.numeric(sprintf('%.15g',x))
		if(is.finite(y) && identical(as.numeric(sprintf('%.15g',y)),y)){
			return(-y)
		}
		x <- x*(1-.Machine$double.eps)
	}
	stop('could not determine a csv safe marker for failed runs\n')
})()

# logLike.quasiEps ####
# One of these is added to logLike.failedRun per year of output a partial run did
# produce, so it can be told apart from one that produced nothing. Too small an
# increment is lost: near the marker's 1.8e308 neighbouring doubles are some
# 2e292 apart. One step of the 15 digit decimal grid there, 1e294, survives both
# the addition and the csv round trip.
logLike.quasiEps <- 10^(floor(log10(abs(logLike.failedRun)))-14)

# logLike.failedRun.max ####
# Anything at or below this is one of the failed run markers above, anything
# greater is a real log likelihood. The headroom takes any number of years a
# partial run might report, with orders of magnitude to spare, and is still
# negligible next to the marker itself.
logLike.failedRun.max <- logLike.failedRun+1e4*logLike.quasiEps

source('naturalsort.R') # used to sort the chunked per var files before reading

source('funRunFRIDA.R')
source('funPlot.R')
source('funParmSpace.R')

wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {
	options(width=as.integer(howWide))
}

