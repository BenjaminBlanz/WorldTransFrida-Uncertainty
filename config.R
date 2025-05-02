# redo all calculations instead of using stored values
redoAllCalc <- F

# parallel ####
#if(!exists('numWorkers')){
#	numWorkers <- min(parallel::detectCores(), 120)
#}
numWorkers <- parallel::detectCores()
# How large the chunks of work are, smaller means more frequent pauses to write out
# itermediate results (and update the diagnostic output).
chunkSizePerWorker <- 100
# tyoe of cluster. PSOCK allows connections across a network
# FORK forks the currently running process, but with copy on write memory
# sharing
clusterType <- 'psock'
# should the workers save their output independently or send it back to the
# main thread.
# If true each worker writes its results to disk in a seperate file. This should be
# much faster than handling all output in a single thread.
writePerWorkerFiles <- TRUE
doNotReturnRunDataSavePerWorkerOnly <- FALSE
perVarOutputTypes <- c('RDS','csv')

#plotting ####
#related things
plotWhileRunning <- F
plotDatWhileRunning <- F
plotDatPerChunWhileRunning <- F
whatToPlot <- tolower('GDP_Real_GDP_in_2021c')
# padding for data plots y axis in share of the data range
yaxPad <- 0.4
# pretty plots
location.plots <- 'figures'
yearsToPlot.names <- c('allYears')#,'1980-2023')
uncertaintiesToPlot <- c('fit uncertainty')#,'noise uncertainty','all uncertainty')
alsoPlotMean.vals <- c(FALSE)
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
plotWeightTypes <- c('completeEqually')#,'logLikelihood')#,'linearly','logCutoff')#,'likelihood') #options are equaly, completeEqually, linear, logCutoff, likelihood
CIsToPlot <- c(0,0.67,0.95)
CIsToPlot.lty <- c('solid','longdash','dotted')#,'dotdash','dotted')
CIsToPlot.lwd <- c(3,1,1)
CIsToPlot.lcol <- c(1,1,1)
CIsToPlot.col <- c(NA,gray(0.7,0.5),gray(0.8,0.5))

calDat.col <- 'red'

# sampling ####
# number of samples for the sobol sequence across all dimensions
numSample <- 1e4
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

# parameger ranges ####
# do we assume or pretend we assume that all residuals are independent.
# I.e. the cov matrix is a diagonal withe the per variable variance on the diagonal
treatVarsAsIndep <- T
# for changing the parm space where should our threshold be.
# The threshold is a ratio between the maximum likelihood parms and the least likely
# parms.
# The parm range will be either increased or decreased to make this happen in each
# parameter.
likeCutoffRatio <- 1000
# tolerance for the search of the likelihood border
rangeTol <- 1e-15
checkBorderErrors <- FALSE
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
allowAssymetricToAvoidZeroRanges <- FALSE
symmetricRangesBoundByAuthors <- TRUE
# should we skip the parameter maximum likelihood estimation and use the default
# frida pars as MLE
if(!exists('skipParMLE')){
	skipParMLE <- T
}

# representative subsample ####
subSample.NumSamplePerVar <- 11
subSample.Ps <- seq(0.5/subSample.NumSamplePerVar,1-0.5/subSample.NumSamplePerVar,
									 length.out=subSample.NumSamplePerVar)
subSample.TargetVars <- c('demographics_real_gdp_per_person')
subSample.sampleJointly <- FALSE

# FRIDA config ####
climateFeedbackSpecFile <- 'ClimateFeedback_On.csv'
climateOverrideSpecFile <- 'ClimateSTAOverride_Off.csv'
policyFileName <- 'policy_EMB.csv'#'policy_100DollarCarbonTax.csv' #'policy_EMB.csv'


# locations and names ####
# location of frida/stella for running
baselocation.frida <-location.frida <- './FRIDAforUncertaintyAnalysis'
baselocation.stella <- location.stella<- './Stella_Simulator_Linux'
# location frida/stella is stored while the above is located in tmpfs
location.frida.storage <- './FRIDAforUncertaintyAnalysis-store'
location.stella.storage <- './Stella_Simulator_Linux-store'
# FRIDA config
# location for setting parameters for FRIDA
# e.g. turnig climate feedbacks on or off
# or policy
location.frida.configs <- './FRIDA-configs'
# location for files used to set up parameter ranges,
# variables to use/not use
# export preferences etc
location.frida.info <- './FRIDA-info'
name.frida_external_ranges <-'frida_external_ranges.csv'
name.frida_info <- 'frida_info.csv'
name.frida_integer_parms <- 'frida_integer_parms.csv'
name.frida_parameter_exclusion_list <- 'frida_parameter_exclusion_list.csv'
# list of variables to exclude from the likelihood calculations
name.frida_variable_exclusion_list <- 'frida_variable_exclusion_list.csv'
name.frida_extra_variables_to_export_list <- 'frida_extra_variables_to_export_list.csv'

# names of files written to the FRIDA Data directory for the running and export
name.fridaExportVarsFile <- 'varsForExport.txt'
name.fridaInputFile <- 'uncertainty_analysis_paramter_values.csv'
name.fridaOutputFile <- 'uncertainty_analysis_exported_variables.csv'


# execute config ####
name.output <- 'dummyNameForSubmitSlurmScriptToOverwrite'
# if this was not run by slurm, the above will not be overwritten and so we set a 
# sensible output name. Otherwise name.output will be set by the slurm submit script
if(name.output=='dummyNameForSubmitSlurmScriptToOverwrite'){
	name.output <- paste0('N-',numSample,
												'-ChS-',chunkSizePerWorker,
												'-LCR-',likeCutoffRatio,
												'-IgB-',ignoreParBounds,
												'-FrB-',forceParBounds,
												'-KcE-',kickParmsErrorRangeDet,
												'-Sym-',symmetricRanges,
												'-AAZ-',allowAssymetricToAvoidZeroRanges,
												'-CFB-',strsplit(tools::file_path_sans_ext(climateFeedbackSpecFile),'_')[[1]][2],
												'-Pol-',tools::file_path_sans_ext(policyFileName),
												'-CTO-',strsplit(tools::file_path_sans_ext(climateOverrideSpecFile),'_')[[1]][2])
}
location.output <- file.path('workOutput',name.output)
location.output.base <- location.output
# tmpfs location for the worker directories to not churn the hard drive
# and be faster
# typical options on linux are /dev/shm or /run/user/####/ where #### is the uid
# if both of these are unavailable use notTMPFS or some other arbitrary location on disk
# tmpfsBaseDir <- paste0('/run/user/',system('id -u',intern = T),'/rwork')
tmpfsBaseDir <- paste0('/dev/shm/',system('id -u',intern = T),'/rwork')
# tmpfsBaseDir <- 'notTMPFS'
tmpfsDir <- file.path(tmpfsBaseDir,name.output)

name.workDir <- paste0('workerDirs-',name.output)
name.workerDirBasename <- 'workDir_'


cat(sprintf('Output folder: %s\n',location.output))
if(file.exists(location.output)){
	cat('  exists\n')
} else {
	dir.create(file.path(location.output),recursive = T,showWarnings = F)
	cat('  created\n')
}
# save the config to the output folder
file.copy('config.R',location.output,overwrite = T)

# run setupTMPFS now, so that location.frida points to the one specific for this
# configuration
# but only do this if the executing process is not a worker running in its own work dir
# we can detect this by the file not existing, as workers do not get this file
if(file.exists('setupTMPFS.R')){
	source('setupTMPFS.R')
} else {
	location.frida <- paste0(baselocation.frida,'-',name.output)
	location.stella <- paste0(baselocation.stella,'-',name.output)
}

# copy slected policy file and climate feedbacks config to frida
cat(sprintf('Copying %s, %s, and %s to the frida directory.\n',
						climateFeedbackSpecFile,
						policyFileName,climateOverrideSpecFile))
file.copy(file.path(location.frida.configs,climateFeedbackSpecFile),
					file.path(location.frida,'Data','climateFeedbackSwitches.csv'),T)
file.copy(file.path(location.frida.configs,policyFileName),
					file.path(location.frida,'Data','policyParameters.csv'),T)
file.copy(file.path(location.frida.configs,climateOverrideSpecFile),
					file.path(location.frida,'Data','ClimateSTAOverride.csv'),T)

