numWorkers <- 2
numSample <- 4e2

# How large the chunks of work are, smaller means more frequent pauses to write out
# itermediate results (and update the diagnostic output).
chunkSizePerWorker <- 100

# by default sobol sequence covers the entire range between min and max with 
# equal density.
# However we might want to ensure that there are similar number of points above and 
# below the Value in our baseline calibration, our prior.
restretchSamplePoints <- F

plotWhileRunning <- T
whatToPlot <- 'GDP_Real_GDP_in_2021c'

location.frida <- './FRIDAforUncertaintyAnalysis'
location.stella<- './Stella_Simulator_Linux'
name.fridaExportVarsFile <- 'varsForExport.txt'
name.fridaInputFile <- 'uncertainty_analysis_paramter_values.csv'
name.fridaOutputFile <- 'uncertainty_analysis_exported_variables.csv'

location.output <- file.path('workOutput',paste0('NumSample-',numSample,
																								 '-chunkSizePerWorker-',chunkSizePerWorker))

minObsForLike <- 35
removeLinearCombinations <- F
imputeMissingVars <- T
# this can 
# 'n'      do not extrapolate
# 'f'      fill in the last good value for all missing values
# 'l##'    linear extrapolation using the first/last ##% of observations
# 'q##'    quadratic extrapolation using the first/last ##% of observations
extrapolateMissingVarMethod <- 'n'
