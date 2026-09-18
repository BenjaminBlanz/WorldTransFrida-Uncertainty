# packages ####
packages.attach <- c(
	'Rmpfr', # arbitrary precision math used to calculate the likelihood from loglikelihood
	'optimx', # interface to various optimizers
	'tictoc', # simple timing measurements
	'SobolSequence', # generates multidimensional sobol sequences
	'lubridate', # deals with times
	'cNORM', # wighted quantiles
	'spatstat.explore', # used for the quantile.density function
	'caret', # to find linear combinations and remove them in the calib dat
	'matrixcalc', # to test positive definitnes of cov matrix
	#'imputeTS', # used for interpolating missing values # Only needed for interpolation after MLE, fails for R v4.4 on Levante
	'data.table', # fast csv reading/writing of the per var files
	'parallel') # for running things in parallel
# only used through :: or by other packages. Attaching R.utils would mask base
# functions such as cat, load and save.
packages.namespace <- c(
	'R.utils', # data.table::fread needs it to read gz files
	'mvtnorm', # multivariate normal likelihood
	'digest', # hashes the inputs a cached parscale determination was computed from
	'spatstat.univar') # weighted quantiles in the plots

# ensurePackages ####
# Installs the missing ones of pkgs from CRAN into the first writable library
# and stops if any remain missing.
ensurePackages <- function(pkgs){
	isInstalled <- function(p){suppressPackageStartupMessages(requireNamespace(p,quietly=T))}
	missing <- pkgs[!vapply(pkgs,isInstalled,logical(1))]
	if(length(missing)==0){
		return(invisible(TRUE))
	}
	cat(sprintf('Installing missing R package%s: %s\n',
							ifelse(length(missing)==1,'','s'),paste(missing,collapse=', ')))
	repos <- getOption('repos')
	repos[repos=='@CRAN@'] <- 'https://cloud.r-project.org'
	if(length(repos)==0){
		repos <- 'https://cloud.r-project.org'
	}
	lib <- .libPaths()[file.access(.libPaths(),2)==0][1]
	if(is.na(lib)){
		lib <- path.expand(strsplit(Sys.getenv('R_LIBS_USER'),.Platform$path.sep)[[1]][1])
		dir.create(lib,recursive=T,showWarnings=F)
		.libPaths(c(lib,.libPaths()))
	}
	try(install.packages(missing,lib=lib,repos=repos))
	missing <- missing[!vapply(missing,isInstalled,logical(1))]
	if(length(missing)>0){
		# cat before stop, so this ends up in the log file
		cat(sprintf('\nMissing R package%s: %s\nInstall with the R used for the runs, on a node with internet access (e.g. a login node):\n  install.packages(c(%s))\n\n',
								ifelse(length(missing)==1,'','s'),paste(missing,collapse=', '),
								paste0("'",missing,"'",collapse=',')))
		stop(sprintf('missing R packages: %s\n',paste(missing,collapse=', ')))
	}
	invisible(TRUE)
}
ensurePackages(c(packages.attach,packages.namespace))
suppressPackageStartupMessages(
	invisible(lapply(packages.attach,library,character.only=T,quietly=T,warn.conflicts=F)))
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

