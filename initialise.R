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
# The log likelihood that a run which failed or produced no usable output is
# marked with. Anything greater than this counts as a complete run.
# The marker has to survive a csv round trip, because the per variable files
# hold 15 significant digits unless perVarFullPrecision is on.
# -.Machine$double.xmax does not survive it: written with 15 digits it rounds up
# past the largest representable double and reads back as -Inf. So the marker is
# the largest magnitude negative value that is a fixed point of that round trip,
# -1.79769313486231e+308, 29 representable doubles below -.Machine$double.xmax.
# It is derived rather than written out so that it stays correct if the printing
# ever changes.
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
# A run that failed part way through is marked with logLike.failedRun plus one
# of these per year of output it did manage to produce, so that a partial run
# can be told apart from one that produced nothing at all.
# The obvious increment, .Machine$double.eps, cannot do that. eps is the spacing
# of the doubles near 1; near the magnitude of the marker, 1.8e308, two
# neighbouring doubles are about 2e292 apart, so adding 2.2e-16 to 1.8e308
# returns 1.8e308 unchanged.
# So add a quasi eps instead: one step of the 15 significant digit decimal grid
# at the magnitude of the marker, 1e294. That is far coarser than the spacing of
# the doubles there, so the increment survives the addition, and it is one full
# step of the grid a csv written without perVarFullPrecision snaps to, so it
# survives being written out and read back in again as well.
# In relative terms it is still nothing: even a thousand of them move the marker
# by 6e-11 of its own magnitude. The result stays around -1.8e308, three hundred
# orders of magnitude below any log likelihood a run could produce.
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

