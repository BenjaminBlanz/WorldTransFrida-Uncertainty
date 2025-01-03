# redo all calculations instead of using stored values
redoAllCalc <- F

# parallel things
if(!exists('numWorkers')){
	numWorkers <- parallel::detectCores()
}
# How large the chunks of work are, smaller means more frequent pauses to write out
# itermediate results (and update the diagnostic output).
chunkSizePerWorker <- 100
# tyoe of cluster. PSOCK allows connections across a network
# FORK forks the currently running process, but with copy on write memory
# sharing
clusterType <- 'psock'
# tmpfs location for the worker directories to not churn the hard drive
# and be faster
# alternative path to /dev/shm would be /run/user/####/ where #### is the uid
# /dev/shm is not executable on some distros use /run/user
tmpfsDir <- paste0('/run/user/',system('id -u',intern = T),'/rwork')

#plotting related things
plotWhileRunning <- F
plotDatWhileRunning <- F
plotDatPerChunWhileRunning <- F
whatToPlot <- tolower('GDP_Real_GDP_in_2021c')
# padding for data plots y axis in share of the data range
yaxPad <- 0.4
# pretty plots
location.plots <- 'figures'
yearsToPlot.names <- c('allYears','1980-2023')
uncertaintiesToPlot <- c('all uncertainty','noise uncertainty','fit uncertainty')
alsoPlotMean.vals <- c(TRUE,FALSE)
mean.lty <- 'solid'
mean.lwd <- 2
mean.col <- 'blue'
alsoPlotDefaultRun.vals <- c(TRUE,FALSE)
def.lty <- 'solid'
def.lwd <- 2
def.col <- 'green'
plotWidth <- 20
plotHeight <- 20
plotUnit <- 'cm'
plotRes <- 150
plotWeightTypes <- c('equaly')#,'linearly','logCutoff')#,'likelihood') #options are equal, linear, logCutoff, likelihood
CIsToPlot <- c(0,0.5,0.95)
CIsToPlot.lty <- c('solid','longdash','dotted')#,'dotdash','dotted')
CIsToPlot.lwd <- c(3,1,1)
CIsToPlot.lcol <- c(1,1,1)
CIsToPlot.col <- c(NA,gray(0.7,0.5),gray(0.8,0.5))

calDat.col <- 'red'

# sampling the parameters ####

# number of samples for the sobol sequence across all dimensions
numSample <- 1439
# by default sobol sequence covers the entire range between min and max with 
# equal density.
# However we might want to ensure that there are similar number of points above and 
# below the Value in our baseline calibration, our prior.
restretchSamplePoints <- F
# For the likelihood we require a positive definite covariance matrix of the residuals.
# A greater number of minObs increase the change of having more complete cases to 
# work with, increasing our ods of a good cov mat.
minObsForLike <- 5
# linear combinations in the residuals will make the matrix singular, i.e. not 
# positive definite. So we remove then
removeLinearCombinations <- F
# To increase the complete cases we can impute missing observations of individual 
# vars (by linear interpolation, smarter later maybe)
imputeMissingVars <- F
# In addition we can try to extrapolate from the calibration data we have to cover more
# years. However this comes with the risk of producing nonsense, check the diagnostic plots!
# So far only the 'n' and 'f' options are implemented.
# 'n'      do not extrapolate
# 'f'      fill in the last good value for all missing values
# 'l##'    linear extrapolation using the first/last ##% of observations
# 'q##'    quadratic extrapolation using the first/last ##% of observations
extrapolateMissingVarMethod <- 'n'

# calculating the likelihoods ####

# do we assume or pretend we assume that all residuals are independent.
# I.e. the cov matrix is a diagonal withe the per variable variance on the diagonal
treatVarsAsIndep <- T
# for changing the parm space where should our threshold be.
# The threshold is a ratio between the maximum likelihood parms and the least likely
# parms.
# The parm range will be either increased or decreased to make this happen in each
# parameter.
likeCutoffRatio <- 200
# tolerance for the search of the likelihood border
rangeTol <- 1e-15
kickParmsErrorRangeDet <- FALSE
kickParmsErrorRangeDet.tolerance <- 1e-2
# further overrides
ignoreParBounds <- FALSE
forceParBounds <- FALSE
# should we make the parameter range be symmetric
# can specify 'Min' or 'Max' to decide if the larger or smaller of the ranges of 
# the parameter should be used to decide the new range. Any other value deactivates this
# feature.
symmetricRanges <- 'Min'
# should we skip the parameter maximum likelihood estimation and use the default
# frida pars as MLE
if(!exists('skipParMLE')){
	skipParMLE <- T
}


# locations and names ####
# location of frida/stella for running
location.frida <- './FRIDAforUncertaintyAnalysis'
location.stella<- './Stella_Simulator_Linux'
# location frida/stella is stored while the above is located in tmpfs
location.frida.storage <- './FRIDAforUncertaintyAnalysis-store'
location.stella.storage <- './Stella_Simulator_Linux-store'

name.fridaExportVarsFile <- 'varsForExport.txt'
name.fridaInputFile <- 'uncertainty_analysis_paramter_values.csv'
name.fridaOutputFile <- 'uncertainty_analysis_exported_variables.csv'
location.output <- file.path('workOutput',paste0('NumSample-',numSample,
																								 '-chunkSizePerWorker-',chunkSizePerWorker,
																								 '-likeCutoffRatio-',likeCutoffRatio,
																								 '-ignoreParBounds-',ignoreParBounds,
																								 '-forceParBounds-',forceParBounds,
																								 '-kickParmsErrorRangeDet-',kickParmsErrorRangeDet,
																								 '-symmetricRanges-',symmetricRanges))
location.output.base <- location.output
cat(sprintf('Output folder: %s\n',location.output))
if(file.exists(location.output)){
	cat('  exists\n')
} else {
	dir.create(file.path(location.output),recursive = T,showWarnings = F)
	cat('  created\n')
}

# save the config to the output folder
file.copy('config.R',location.output,overwrite = T)

