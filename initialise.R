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
})

source('funRunFRIDA.R')
source('funPlot.R')
source('funParmSpace.R')

wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {
	options(width=as.integer(howWide))
}

