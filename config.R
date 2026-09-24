# redo all calculations instead of using stored values
redoAllCalc <- F
# Retry the parscales a previous run could not determine. Those are the most
# expensive parameters in the determination, since failing means having swept
# every order of magnitude, and the answer does not change unless the model or
# the calibration data has. Off, so a rerun keeps the previous verdict.
redoFailedParscales <- F

# parallel ####
#if(!exists('numWorkers')){
#	numWorkers <- min(parallel::detectCores(), 120)
#}
numWorkers <- parallel::detectCores()
# Workers used to merge the per chunk files into one file per variable.
# The merge streams the chunk files into the final file, so for csv only output
# its memory footprint is a fixed buffer per worker and this can be set high.
# If 'RDS' is among the perVarOutputTypes, budget roughly
# numSample * length(outputDataYears) * 8 bytes per worker on top of that,
# as building an RDS requires the whole variable to be in memory once.
numWorkersFileMerge <- numWorkers
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
# When the workers write their own files, the run data is already on disk by the
# time they return. With this set the workers return only the parameter index and
# the log likelihood, which is all the main thread and the resume logic need,
# instead of a second copy of the run data in workUnit-<i>.RDS.
# Set it to FALSE if you need whole runs back, as verificationCases.R does via
# loadClusterRuns. Plotting while running does not need it, that path already
# turns writePerWorkerFiles off.
doNotReturnRunDataSavePerWorkerOnly <- TRUE
# Format(s) of the *final* one file per variable results.
# The per chunk intermediates the workers write are always plain uncompressed
# csv, mergePerVarFiles derives every requested final format from those.
# Outputting csv files only massively reduces the amount of memory needed in the 
# merging step. If enabling RDS files make sure to reduce the number of workers
# used in the merge step.
# Allowed options: c('csv','RDS')
perVarOutputTypes <- c('csv')
# gzip the final per variable csv files. The per chunk intermediates stay
# uncompressed either way, so that the merge can concatenate them byte wise.
compressCsv <- TRUE
# compression of the final per variable RDS files. FALSE is markedly faster to
# write and to read back, at the cost of much larger files.
perVarRdsCompress <- TRUE
# write the doubles of the per chunk csv intermediates with enough digits to
# reproduce them exactly. The merged files, RDS included, are built from those
# intermediates, so with this off they carry a relative error of some 1e-15
# against what the model produced. Off writes the chunks about three times
# faster and 10% smaller. The marker for a failed run (logLike.failedRun in
# initialise.R) survives either setting.
# The eps values marking incomplete years of a run are large enough to survive the
# lower precision, and no result needs more, so off is safe.
perVarFullPrecision <- FALSE

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
rs.lty <- 'solid'
rs.lwd <- '1'
rs.col <- 'black'
calDat.col <- 'red'
alsoPlotRepSample.vals <- c(FALSE,TRUE)
repSample.lty <- 'solid'
repSample.lwd <- 1
repSample.col <- 'black'
shareOfYearsThatYlimShouldAdjustTo <- 0.9

# sampling ####
# file that contains the baseline parametrisation of the model
# which is applied before the sample points. I.e. sample points override the
# baseline parametrisation, but if a parameter is not present in the samplePoints, but
# is present in the baseline parms it will be used.
# If this is left empty, no baseline parms are applied (the ones within the model
# files will be used).
name.baselineParmFile <- NA
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
# uniroot's tolerance in the border search, as a multiple of the parameter's own
# parscale. The likelihood comes from a stella run read back out of a csv with
# some loss of precision, so a tolerance far below what the model resolves is one
# it cannot answer to. NA selects an absolute 1e-16.
rangeRootTol <- 1e-4
# iteration limit for that same search
rangeRootMaxIter <- 60
# Should we drop parameters for which we can not determine the parameter scale?
# This is likely because they do not affect the run.
# However we often run EMB first, where policy related parameters have no effect, 
# but then use those sample points to run scenarios with policies. Dropping these
# parameters would undersample the uncertainty of the policies as their uncertain 
# parameters would have been dropped in emb. So only set this to true if you are 
# certain the samplePoints won't be reused for experiments with other specified 
# policies.
# FALSE keeps them and samples them over the ranges their authors gave them in
# frida_info, the same fallback a failed min or max border determination gets.
# Their border determination is skipped, it cannot succeed without a parscale.
# TRUE drops them and writes the names to frida_parameter_exclusion_list.csv.
kickParmsParScaleDet <-FALSE
# Should we check for errors in determining the likelihood border
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
# Should parameters whose range comes from an external override in
# frida_external_ranges.csv be symmetrified? The override is a deliberate
# statement of the range to sample, so by default it is used as given.
symmetrifyExternalRanges <- FALSE
# Should parameters that fell back to their author range be symmetrified?
# Applies to a parameter with a not determined border in either direction,
# whether the border search failed or was skipped for want of a parscale.
# Symmetrifying these is what collapses ranges to zero width, because the
# parameter value often sits exactly on the author bound it fell back to.
symmetrifyFallbackAuthorRanges <- TRUE
# should we skip the parameter maximum likelihood estimation and use the default
# frida pars as MLE
if(!exists('skipParMLE')){
	skipParMLE <- T
}

# representative subsample ####
subSample.NumSamplePerVar <- 11
subSample.Ps <- seq(0.5/subSample.NumSamplePerVar,1-0.5/subSample.NumSamplePerVar,
										length.out=subSample.NumSamplePerVar)
subSample.TargetVars <- c('energy_balance_model_surface_temperature_anomaly')


# FRIDA config ####
climateFeedbackSpecFile <- 'ClimateFeedback_On.csv'
climateOverrideSpecFile <- 'ClimateSTAOverride_Off.csv'
climateOverrideSpecFileTS <- 'ClimateSTAOverrideTS_none.csv'
policyFileName <- 'policy_EMB.csv'#'policy_100DollarCarbonTax.csv' #'policy_EMB.csv'


# locations and names ####
# branch of frida to use
name.frida_branch <- 'main'
# location of frida/stella for running
baselocation.frida <-location.frida <- './FRIDAforUncertaintyAnalysis'
baselocation.stella <- location.stella <- './Stella_Simulator_Linux'
# git checkout the model files come from. FRIDAforUncertaintyAnalysis itself is an
# rsync copy without the .git directory, so the version of frida in use has to be
# read from here. Maintained by uncertainity_update_frida.sh.
location.frida.git <- paste0(baselocation.frida,'Git')
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

# file containing names of uncertainty parameters and ranges set by the authors
# pre v3.1 this is a file that has to be manually updated.
# post v3.1 this is a file that is generated as part of the frida files and can be 
# found in the FRIDAforUncertaintyAnalysis folder. This step of the config tests
# if that file exists and sets it for use, otherwise falling back to the manually
# updated file.
suppressWarnings(rm('name.frida_info'))
# If you want to use a user specified frida_info file, e.g. for certain uncertainty analysis,
# uncomment the following lines to specify name.frida_info and place the file in location.frida.info
# 
# name.frida_info <- 'frida_info_override.csv'
if(exists('name.frida_info')){
	frida_info_type <- 'user'
} else {
	if(file.exists(file.path(baselocation.frida,'Parameter Info.csv'))){
		name.frida_info <- 'link_to_frida_info_from_model_repo.csv'
		# every run started from this directory reads the link, so it is never
		# removed, only replaced by renaming a new link over it
		local({
			link <- file.path(location.frida.info,name.frida_info)
			target <- sprintf('../%s/Parameter Info.csv',baselocation.frida)
			if(!identical(unname(Sys.readlink(link)),target)){
				tmpLink <- paste0(link,'.',Sys.info()[['nodename']],'.',Sys.getpid())
				unlink(tmpLink)
				file.symlink(target,tmpLink)
				file.rename(tmpLink,link)
			}
		})
		frida_info_type <- 'StellaExport'
	} else {
		name.frida_info <- 'frida_info_preV3.csv'
		frida_info_type <- 'OldStyleFromBilly'
	}
}
name.frida_integer_parms <- 'frida_integer_parms.csv'
name.frida_parameter_exclusion_list <- 'frida_parameter_exclusion_list.csv'
# list of variables to exclude from the likelihood calculations
name.frida_variable_exclusion_list <- 'frida_variable_exclusion_list.csv'
name.frida_extra_variables_to_export_list <- 'frida_extra_variables_to_export_list.csv'

# names of files written to the FRIDA Data directory for the running and export
name.fridaExportVarsFile <- 'varsForExport.txt'
name.fridaInputFile <- 'uncertainty_analysis_paramter_values.csv'
name.fridaOutputFile <- 'uncertainty_analysis_exported_variables.csv'
# everything this analysis writes into the FRIDA Data directory: the scenario of
# the run and the export spec. The rest of the frida directory is the model.
name.fridaAnalysisDataFiles <- c('climateFeedbackSwitches.csv','policyParameters.csv',
																 'ClimateSTAOverride.csv','ClimateSTAOverrideTS.csv',
																 name.fridaExportVarsFile,name.fridaInputFile,name.fridaOutputFile)


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
# copies of the run scripts, the config and the input files of the run
location.output.runScripts <- file.path(location.output,'runScriptsAndConfiguration')
# tmpfs location for the worker directories to not churn the hard drive
# and be faster
# typical options on linux are /dev/shm or /run/user/####/ where #### is the uid
# if both of these are unavailable use notTMPFS or some other arbitrary location on disk
# tmpfsBaseDir <- paste0('/run/user/',system('id -u',intern = T),'/rwork')
tmpfsBaseDir <- paste0('/dev/shm/',system('id -u',intern = T),'/rwork')
# tmpfsBaseDir <- 'notTMPFS'
origTmpfsDir <- tmpfsDir <- file.path(tmpfsBaseDir,name.output)

origName.workDir <- name.workDir <- paste0('workerDirs-',name.output)
name.workerDirBasename <- 'workDir_'


cat(sprintf('Output folder: %s\n',location.output))
if(file.exists(location.output)){
	cat('  exists\n')
} else {
	dir.create(file.path(location.output),recursive = T,showWarnings = F)
	cat('  created\n')
}
# save the config to the output folder
# A run that only post processes an existing ensemble sets recordRunProvenance to
# FALSE and leaves the record of the ensemble that produced the data alone.
if(!exists('recordRunProvenance')){
	recordRunProvenance <- TRUE
}
if(recordRunProvenance){
	dir.create(location.output.runScripts,recursive = T,showWarnings = F)
	file.copy('config.R',location.output.runScripts,overwrite = T)
}
# which config file this run is using. The submit script rewrites every mention
# of config.R in its copy of this file, so this ends up naming the copy, which is
# what the run metadata needs to diff against the default config.
name.configFile <- 'config.R'

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
file.copy(file.path(location.frida.configs,climateOverrideSpecFileTS),
					file.path(location.frida,'Data','ClimateSTAOverrideTS.csv'),T)

# record the run metadata ####
# Which version of frida, which version of these scripts and which config
# settings produced a set of results is not visible from the output folder name,
# so write it into the folder itself. The completion of the ensemble is appended
# to the same file once it has run, by funAppendRunCompletionSummary.
# Only the process that runs the analysis does this. Workers re-source this config
# from their own work dirs, where writing the file would neither be correct nor
# useful. They are the ones without setupTMPFS.R, same test as above.
if(recordRunProvenance&&file.exists('setupTMPFS.R')&&exists('funWriteRunMetadataFile',mode='function')){
	fridaVersion <- funWriteRunMetadataFile(
		location.output,location.frida.git,location.frida,name.output,
		configFile=name.configFile,
		exclude=name.fridaAnalysisDataFiles,
		specFiles=c(policyFileName=policyFileName,
								climateFeedbackSpecFile=climateFeedbackSpecFile,
								climateOverrideSpecFile=climateOverrideSpecFile,
								climateOverrideSpecFileTS=climateOverrideSpecFileTS),
		baselineParmFile=name.baselineParmFile,
		location.frida.configs=location.frida.configs,
		inputFiles=c(file.path(location.frida.configs,
													 c(policyFileName,climateFeedbackSpecFile,climateOverrideSpecFile,
													 	climateOverrideSpecFileTS,name.baselineParmFile)),
								 file.path(location.frida.info,
								 					c(name.frida_info,name.frida_integer_parms,name.frida_external_ranges,
								 						name.frida_parameter_exclusion_list,name.frida_variable_exclusion_list,
								 						name.frida_extra_variables_to_export_list))),
		location.inputs=file.path(location.output.runScripts,'input'))
	if(fridaVersion['commit']=='noGit'){
		cat('FRIDA version: could not be determined, see runMetadata.txt\n')
	} else {
		cat(sprintf('FRIDA version: %s (%s, %s)\n',
								substr(fridaVersion['commit'],1,7),
								fridaVersion['branch'],fridaVersion['date']))
	}
}
