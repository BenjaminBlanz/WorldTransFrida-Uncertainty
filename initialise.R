suppressPackageStartupMessages({
	library(caret,quietly=T,warn.conflicts = F) # to find linear combinations and remove them in the calib dat
	library(matrixcalc,quietly=T,warn.conflicts = F) # to test positive definitnes of cov matrix
	library(imputeTS,quietly=T,warn.conflicts = F) # used for interpolating missing values
	library(Rmpfr,quietly=T,warn.conflicts = F) # use to calculate the likelihood from loglikelihood
	library(parallel,quietly=T,warn.conflicts = F)
	library(optimx,quietly=T,warn.conflicts = F)
	library(MASS,quietly=T,warn.conflicts = F)
	library(spatstat.explore,quietly=T,warn.conflicts = F) # used for the quantile.density function
	library(tictoc)
})

source('config.R')
source('funRunFRIDA.R')
source('funPlot.R')
source('config.R')
source('funParmSpace.R')


