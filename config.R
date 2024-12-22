numWorkers <- parallel::detectCores()
numSample <- 1e5

# How large the chunks of work are, smaller means more frequent pauses to write out
# itermediate results (and update the diagnostic output).
chunkSizePerWorker <- 100

# by default sobol sequence covers the entire range between min and max with 
# equal density.
# However we might want to ensure that there are similar number of points above and 
# below the Value in our baseline calibration, our prior.
restretchSamplePoints <- F

plotWhileRunning <- T
whatToPlot <- tolower('GDP_Real_GDP_in_2021c')

location.frida <- './FRIDAforUncertaintyAnalysis'
location.stella<- './Stella_Simulator_Linux'
name.fridaExportVarsFile <- 'varsForExport.txt'
name.fridaInputFile <- 'uncertainty_analysis_paramter_values.csv'
name.fridaOutputFile <- 'uncertainty_analysis_exported_variables.csv'

location.output <- file.path('workOutput',paste0('NumSample-',numSample,
																								 '-chunkSizePerWorker-',chunkSizePerWorker))
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
