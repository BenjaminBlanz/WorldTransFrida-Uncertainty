
location.frida <- './FRIDAforUncertaintyAnalysis'
location.stella<- './Stella_Simulator_Linux'
name.fridaExportVarsFile <- 'varsForExport.txt'
name.fridaInputFile <- 'uncertainty_analysis_paramter_values.csv'
name.fridaOutputFile <- 'uncertainty_analysis_exported_variables.csv'

location.output <- file.path('workOutput',paste0('NumSample-',numSample,
																								 '-chunkSizePerWorker-',chunkSizePerWorker))
location.output.base <- location.output
dir.create(file.path(location.output),recursive = T,showWarnings = F)



# parallel things
numWorkers <- parallel::detectCores()
# How large the chunks of work are, smaller means more frequent pauses to write out
# itermediate results (and update the diagnostic output).
chunkSizePerWorker <- 100

#plotting related things
plotWhileRunning <- T
plotDatWhileRunning <- F
whatToPlot <- tolower('GDP_Real_GDP_in_2021c')
# padding for data plots y axis in share of the data range
yaxPad <- 0.4

# sampling the parameters ####

# number of samples for the sobol sequence across all dimensions
numSample <- 5e3
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
likeCutoffRatio <- 100
# should we skip the parameter maximum likelihood estimation and use the default
# frida pars
if(!exists('skipParMLE')){
	skipParMLE <- F
}

# save the config to the output folder
file.copy('config.R',location.output)
