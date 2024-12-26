numWorkers <- parallel::detectCores()
numSample <- 5e3

# How large the chunks of work are, smaller means more frequent pauses to write out
# itermediate results (and update the diagnostic output).
chunkSizePerWorker <- 100

# by default sobol sequence covers the entire range between min and max with 
# equal density.
# However we might want to ensure that there are similar number of points above and 
# below the Value in our baseline calibration, our prior.
restretchSamplePoints <- F

plotWhileRunning <- T
plotDatWhileRunning <- F
whatToPlot <- tolower('GDP_Real_GDP_in_2021c')

location.frida <- './FRIDAforUncertaintyAnalysis'
location.stella<- './Stella_Simulator_Linux'
name.fridaExportVarsFile <- 'varsForExport.txt'
name.fridaInputFile <- 'uncertainty_analysis_paramter_values.csv'
name.fridaOutputFile <- 'uncertainty_analysis_exported_variables.csv'

location.output <- file.path('workOutput',paste0('NumSample-',numSample,
																								 '-chunkSizePerWorker-',chunkSizePerWorker))
location.output.base <- location.output
dir.create(file.path(location.output),recursive = T,showWarnings = F)
# For the likelihood we require a positive definite covariance matrix of the residuals.
# A greater number of minObs increase the change of having more complete cases to 
# work with, increasing our ods of a good cov mat.
minObsForLike <- 40
# linear combinations in the residuals will make the matrix singular, i.e. not 
# positive definite. So we remove then
removeLinearCombinations <- T
# To increase the complete cases we can impute missing observations of individual 
# vars (by linear interpolation, smarter later maybe)
imputeMissingVars <- T
# In addition we can try to extrapolate from the calibration data we have to cover more
# years. However this comes with the risk of producing nonsense, check the diagnostic plots!
# So far only the 'n' and 'f' options are implemented.
# 'n'      do not extrapolate
# 'f'      fill in the last good value for all missing values
# 'l##'    linear extrapolation using the first/last ##% of observations
# 'q##'    quadratic extrapolation using the first/last ##% of observations
extrapolateMissingVarMethod <- 'n'


# for changing the parm space where should our threshold be.
# The threshold is a ratio between the maximum likelihood parms and the least likely
# parms.
# The parm range will be either increased or decreased to make this happen in each
# parameter.
likeThresholdRatio <- 100

# padding for data plots y axis in share of the data range
yaxPad <- 0.4

file.copy('config.R',location.output)
