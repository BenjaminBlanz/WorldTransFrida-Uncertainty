suppressPackageStartupMessages({
	library(Rmpfr,quietly=T,warn.conflicts = F) # arbitrary precision math used to calculate the likelihood from loglikelihood
	library(optimx,quietly=T,warn.conflicts = F) # interface to various optimizers
	library(tictoc,quietly=T,warn.conflicts = F) # simple timing measurements
	library(SobolSequence,quietly=T,warn.conflicts = F) # generates multidimensional sobol sequences
	library(lubridate,quietly=T,warn.conflicts = F) # deals with times
	library(cNORM,quietly=T,warn.conflicts = F) # wighted quantiles
	library(spatstat.explore,quietly=T,warn.conflicts = F) # used for the quantile.density function
	library(parallel,quietly = T,warn.conflicts = F) # for parallel runs
	library(paletteer,quietly = T,warn.conflicts = F) # for color scales
	library(optparse,quietly = T,warn.conflicts = F) # for passing input arguments to scripts
})

source('funRunFRIDA.R')
source('funPlot.R')
source('funParmSpace.R')

quietgc <- function(){
	sink(file='/dev/null')
	gc()
	sink()	
}

filterSpecsAreEqual <- function(filterSpec1,filterSpec2){
	if(length(filterSpec1)!=length(filterSpec2)){
		return(F)
	}
	for(entry in names(filterSpec1)){
		if(length(filterSpec1[[entry]])!=length(filterSpec2[[entry]])){
			return(F)
		}
		for(subentry.i in 1:length(filterSpec1[[entry]])){
			if(filterSpec1[[entry]][subentry.i]!=filterSpec2[[entry]][subentry.i]){
				return(F)
			}
		}
	}
	return(T)
}
