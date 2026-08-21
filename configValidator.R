#
# Validates that everything specified in config.R is actually present and
# usable, before an expensive run is started.
#
# Source this after initialise.R and config.R, e.g.
#   source('initialise.R')
#   source('config.R')
#   source('configValidator.R')
#
# Why this exists:
# Most of the config is applied either by copying files into the frida
# directory (file.copy) or by reading them lazily, much later in the run.
# Both fail silently. A missing policy file for example leaves whatever
# policyParameters.csv was in the frida directory before in place, and the run
# then happily produces results for the wrong policy, or for no policy at all,
# without ever complaining.
# This script turns all of those silent failures into a single hard error that
# lists everything that is wrong with the current configuration.
#

cat('Validating config...\n')

# environment this validator was sourced into. Usually the global env, but
# using this rather than globalenv() keeps the validator usable from within
# local() or a function.
cfgVal.env <- environment()

cfgVal.errors <- c()
cfgVal.warnings <- c()

cfgVal.error <- function(...){
	cfgVal.errors <<- c(cfgVal.errors,sprintf(...))
	invisible(FALSE)
}
cfgVal.warning <- function(...){
	cfgVal.warnings <<- c(cfgVal.warnings,sprintf(...))
	invisible(FALSE)
}

# report and abort ####
# prints everything found and stops if there were any errors.
cfgVal.report <- function(){
	if(length(cfgVal.warnings)>0){
		cat(sprintf('\nConfig validation warnings (%i):\n',length(cfgVal.warnings)))
		cat(paste0('  - ',cfgVal.warnings,collapse='\n'))
		cat('\n')
	}
	if(length(cfgVal.errors)>0){
		cat(sprintf('\nConfig validation FAILED with %i error%s:\n',
								length(cfgVal.errors),ifelse(length(cfgVal.errors)==1,'','s')))
		cat(paste0('  - ',cfgVal.errors,collapse='\n'))
		cat('\n\n')
		# cat before stop, so the above ends up in the log file, as the error
		# message itself goes to stderr, which is not sunk.
		stop(sprintf('config validation failed with %i error(s), see above.\n',
								 length(cfgVal.errors)))
	}
	cat('Config validation passed.\n')
	invisible(TRUE)
}

# helpers ####
cfgVal.has <- function(name){
	exists(name,envir=cfgVal.env)
}
cfgVal.get <- function(name,default=NULL){
	if(cfgVal.has(name)){
		get(name,envir=cfgVal.env)
	} else {
		default
	}
}
# is the value of a config variable a usable, non empty file name
cfgVal.isSet <- function(name){
	val <- cfgVal.get(name)
	!is.null(val) && length(val)==1 && !is.na(val) && nzchar(as.character(val))
}

cfgVal.checkDir <- function(path,what,writable=FALSE){
	if(is.null(path)||length(path)!=1||is.na(path)||!nzchar(path)){
		return(cfgVal.error('%s: no location specified in the config',what))
	}
	if(!dir.exists(path)){
		return(cfgVal.error('%s: directory not found: %s',what,path))
	}
	if(writable&&file.access(path,2)!=0){
		return(cfgVal.error('%s: directory is not writable: %s',what,path))
	}
	invisible(TRUE)
}

# note that file.exists is FALSE for broken symlinks, which is exactly what we
# want here, name.frida_info may be a link into the frida model repo.
cfgVal.checkFile <- function(path,what,allowEmpty=FALSE){
	if(is.null(path)||length(path)!=1||is.na(path)||!nzchar(path)){
		return(cfgVal.error('%s: no file specified in the config',what))
	}
	if(!file.exists(path)){
		return(cfgVal.error('%s: file not found: %s',what,path))
	}
	if(dir.exists(path)){
		return(cfgVal.error('%s: is a directory, not a file: %s',what,path))
	}
	if(file.access(path,4)!=0){
		return(cfgVal.error('%s: file is not readable: %s',what,path))
	}
	if(!allowEmpty&&file.size(path)==0){
		return(cfgVal.error('%s: file is empty: %s',what,path))
	}
	invisible(TRUE)
}

# check that a csv has the columns the code reading it relies on
cfgVal.checkCsvColumns <- function(path,what,requiredCols){
	if(!file.exists(path)||file.size(path)==0){
		return(invisible(FALSE))
	}
	header <- try(read.csv(path,nrows=1),silent=TRUE)
	if(inherits(header,'try-error')){
		return(cfgVal.error('%s: can not be read as csv: %s',what,path))
	}
	missingCols <- requiredCols[!requiredCols%in%colnames(header)]
	if(length(missingCols)>0){
		return(cfgVal.error('%s: missing column%s %s in %s',what,
												ifelse(length(missingCols)==1,'','s'),
												paste0("'",missingCols,"'",collapse=', '),path))
	}
	invisible(TRUE)
}

cfgVal.checkFlag <- function(name){
	if(!cfgVal.has(name)){
		return(cfgVal.error('config variable %s is not defined',name))
	}
	val <- cfgVal.get(name)
	if(!is.logical(val)||length(val)!=1||is.na(val)){
		return(cfgVal.error("%s must be a single TRUE/FALSE, is '%s'",
												name,paste0(format(val),collapse=',')))
	}
	invisible(TRUE)
}

cfgVal.checkNumber <- function(name,min=-Inf,max=Inf,integer=FALSE){
	if(!cfgVal.has(name)){
		return(cfgVal.error('config variable %s is not defined',name))
	}
	val <- cfgVal.get(name)
	if(!is.numeric(val)||length(val)!=1||is.na(val)||!is.finite(val)){
		return(cfgVal.error("%s must be a single finite number, is '%s'",
												name,paste0(format(val),collapse=',')))
	}
	if(val<min||val>max){
		return(cfgVal.error('%s must be in [%s,%s], is %s',name,
												format(min),format(max),format(val)))
	}
	if(integer&&val!=round(val)){
		return(cfgVal.error('%s must be a whole number, is %s',name,format(val)))
	}
	invisible(TRUE)
}

# prerequisites ####
# without these there is nothing sensible to validate, so stop right away
# instead of reporting hundreds of follow up errors.
if(!cfgVal.has('location.output')){
	stop('configValidator.R: run config.R before this script\n')
}
if(!exists('cleanNames',mode='function')||!exists('prepareSampleParms',mode='function')){
	stop('configValidator.R: run initialise.R before this script\n')
}

cfgVal.requiredVars <- c(
	# parallel and sampling
	'numWorkers','chunkSizePerWorker','clusterType','numSample','likeCutoffRatio',
	'perVarOutputTypes','symmetricRanges','minObsForLike',
	'extrapolateMissingVarMethod',
	# frida config files
	'climateFeedbackSpecFile','policyFileName',
	'climateOverrideSpecFile','climateOverrideSpecFileTS',
	# locations
	'baselocation.frida','baselocation.stella','location.frida','location.stella',
	'location.frida.configs','location.frida.info',
	'name.output','location.output',
	# frida info files
	'frida_info_type','name.frida_info','name.frida_integer_parms',
	'name.frida_external_ranges','name.frida_parameter_exclusion_list',
	'name.frida_variable_exclusion_list',
	'name.frida_extra_variables_to_export_list',
	# files written into the frida Data directory
	'name.fridaExportVarsFile','name.fridaInputFile','name.fridaOutputFile',
	# representative subsample
	'subSample.NumSamplePerVar','subSample.TargetVars')
cfgVal.missingVars <- cfgVal.requiredVars[!sapply(cfgVal.requiredVars,cfgVal.has)]
if(length(cfgVal.missingVars)>0){
	for(cfgVal.v in cfgVal.missingVars){
		cfgVal.error('config variable %s is not defined',cfgVal.v)
	}
	# everything below depends on these, so do not even try
	cfgVal.report()
}

# values ####
cfgVal.checkNumber('numWorkers',min=1,integer=TRUE)
if(cfgVal.has('numWorkersFileMerge')){
	cfgVal.checkNumber('numWorkersFileMerge',min=1,integer=TRUE)
}
cfgVal.checkNumber('chunkSizePerWorker',min=1,integer=TRUE)
cfgVal.checkNumber('numSample',min=1,integer=TRUE)
if(isTRUE(cfgVal.checkNumber('likeCutoffRatio',min=0))&&likeCutoffRatio<=1){
	cfgVal.error('likeCutoffRatio must be greater than 1, is %s',format(likeCutoffRatio))
}
if(cfgVal.has('rangeTol')){
	cfgVal.checkNumber('rangeTol',min=0)
}
cfgVal.checkNumber('minObsForLike',min=1,integer=TRUE)
cfgVal.checkNumber('subSample.NumSamplePerVar',min=1,integer=TRUE)
for(cfgVal.v in c('redoAllCalc','writePerWorkerFiles','doNotReturnRunDataSavePerWorkerOnly',
									'compressCsv','perVarRdsCompress','perVarFullPrecision',
									'plotWhileRunning','plotDatWhileRunning','plotDatPerChunWhileRunning',
									'restretchSamplePoints','imputeMissingVars','removeLinearCombinations',
									'treatVarsAsIndep','checkBorderErrors','kickParmsErrorRangeDet',
									'ignoreParBounds','forceParBounds',
									'allowAssymetricToAvoidZeroRanges','symmetricRangesBoundByAuthors',
									'skipParMLE')){
	if(cfgVal.has(cfgVal.v)){
		cfgVal.checkFlag(cfgVal.v)
	}
}
if(length(clusterType)!=1||!clusterType%in%c('psock','fork')){
	cfgVal.error("clusterType must be 'psock' or 'fork', is '%s'",
							 paste0(as.character(clusterType),collapse=','))
}
if(length(perVarOutputTypes)<1||!all(perVarOutputTypes%in%c('RDS','csv'))){
	cfgVal.error("perVarOutputTypes must be a non empty subset of c('RDS','csv'), is c(%s)",
							 paste0("'",perVarOutputTypes,"'",collapse=','))
}
if(!is.character(symmetricRanges)||length(symmetricRanges)!=1||is.na(symmetricRanges)){
	cfgVal.error('symmetricRanges must be a single string')
} else if(!symmetricRanges%in%c('Min','Max','none','None','FALSE','off')){
	cfgVal.warning("symmetricRanges is '%s', which is neither 'Min' nor 'Max', ranges will not be symmetrified",
								 symmetricRanges)
}
if(length(extrapolateMissingVarMethod)!=1||
	 !grepl('^([nf]|[lq][0-9]+)$',extrapolateMissingVarMethod)){
	cfgVal.error("extrapolateMissingVarMethod must be 'n', 'f', 'l##' or 'q##', is '%s'",
							 paste0(as.character(extrapolateMissingVarMethod),collapse=','))
} else if(!extrapolateMissingVarMethod%in%c('n','f')){
	cfgVal.warning("extrapolateMissingVarMethod '%s' is not implemented, only 'n' and 'f' are",
								 extrapolateMissingVarMethod)
}
if(isTRUE(ignoreParBounds)&&isTRUE(forceParBounds)){
	cfgVal.error('ignoreParBounds and forceParBounds are both TRUE, these contradict each other')
}
if(isTRUE(imputeMissingVars)){
	# initialise.R does not load imputeTS anymore, na_interpolation would not be found
	if(!exists('na_interpolation',mode='function')){
		cfgVal.error('imputeMissingVars is TRUE but na_interpolation (imputeTS) is not available')
	}
}

# locations ####
cfgVal.checkDir(location.frida.configs,'location.frida.configs')
cfgVal.checkDir(location.frida.info,'location.frida.info',writable=TRUE) # exclusion list is written here
cfgVal.checkDir(location.output,'location.output',writable=TRUE)
cfgVal.checkDir(baselocation.frida,'baselocation.frida')
cfgVal.checkDir(baselocation.stella,'baselocation.stella')
# location.frida/location.stella are the per run copies/links set up by
# setupTMPFS.R, they are what is actually run
cfgVal.fridaOK <- cfgVal.checkDir(location.frida,'location.frida (run copy of frida)')
cfgVal.checkDir(location.stella,'location.stella (run copy of the stella simulator)')
if(cfgVal.has('name.workDir')&&!file.exists(name.workDir)){
	cfgVal.warning('worker directory %s does not exist, it will be created by setupTMPFS.R',
								 name.workDir)
}

# frida model ####
if(isTRUE(cfgVal.fridaOK)){
	cfgVal.checkFile(file.path(location.frida,'FRIDA.stmx'),'frida model')
	cfgVal.checkDir(file.path(location.frida,'Data'),'frida Data directory',writable=TRUE)
	cfgVal.checkFile(file.path(location.frida,'Data','Calibration Data.csv'),'calibration data')
	cfgVal.checkFile(file.path(location.frida,'Data','frida_input_data.csv'),'frida input data')
	# the names of the files we write for and read back from stella have to be
	# the ones the model actually imports/exports, otherwise stella silently
	# runs on whatever was in the Data directory before.
	if(file.exists(file.path(location.frida,'FRIDA.stmx'))){
		cfgVal.stmx <- paste0(readLines(file.path(location.frida,'FRIDA.stmx'),warn=FALSE),collapse='\n')
		for(cfgVal.v in c('name.fridaExportVarsFile','name.fridaInputFile','name.fridaOutputFile')){
			if(!grepl(cfgVal.get(cfgVal.v),cfgVal.stmx,fixed=TRUE)){
				cfgVal.error("%s <- '%s' is not referenced by %s, the model does not use that file",
										 cfgVal.v,cfgVal.get(cfgVal.v),file.path(location.frida,'FRIDA.stmx'))
			}
		}
	}
}
cfgVal.stellaBin <- file.path(location.stella,'stella_simulator')
if(!file.exists(cfgVal.stellaBin)){
	cfgVal.error('stella simulator not found: %s',cfgVal.stellaBin)
} else if(file.access(cfgVal.stellaBin,1)!=0){
	cfgVal.error('stella simulator is not executable: %s',cfgVal.stellaBin)
}

# frida version ####
# The version of frida used is recorded in fridaVersion.txt in the output folder.
# That record is only as good as the git checkout it is read from.
if(cfgVal.has('location.frida.git')){
	if(!file.exists(file.path(location.frida.git,'.git'))){
		cfgVal.warning('no git checkout at %s, the frida version can not be determined and will be recorded as noGit',
									 location.frida.git)
	} else if(isTRUE(cfgVal.fridaOK)&&exists('funFridaCheckoutDiff',mode='function')){
		# the model directory is an rsync copy of the checkout, if the two have
		# drifted apart the recorded commit does not describe the model being run.
		# The edits uncertainity_update_frida.sh makes to FRIDA.stmx do not count,
		# funFridaCheckoutDiff sorts that out.
		cfgVal.checkoutDiff <- funFridaCheckoutDiff(
			location.frida,location.frida.git,
			exclude=c('climateFeedbackSwitches.csv','policyParameters.csv',
								'ClimateSTAOverride.csv','ClimateSTAOverrideTS.csv',
								name.fridaExportVarsFile,name.fridaInputFile,name.fridaOutputFile))
		if(cfgVal.checkoutDiff$verdict=='differs'){
			cfgVal.warning('the model in %s does not match the checkout in %s, so the recorded frida version will not describe the model being run: %s',
										 location.frida,location.frida.git,
										 funFridaCheckoutDiffNote(cfgVal.checkoutDiff))
		}
	}
}

# frida config files ####
# These are copied into the frida Data directory by config.R. file.copy fails
# silently, so check both that the source exists and that the copy actually
# arrived with the same content. Otherwise a run inherits the settings of
# whatever ran here before.
# An empty spec file is valid and means stella runs without any changes from
# it, e.g. ClimateSTAOverrideTS_none.csv for no override at all. It is only
# reported as a warning, so that a file that is empty by accident is still
# visible in the log.
cfgVal.configFiles <- list(
	list(var='climateFeedbackSpecFile', dest='climateFeedbackSwitches.csv'),
	list(var='policyFileName',          dest='policyParameters.csv'),
	list(var='climateOverrideSpecFile', dest='ClimateSTAOverride.csv'),
	list(var='climateOverrideSpecFileTS',dest='ClimateSTAOverrideTS.csv'))
for(cfgVal.cf in cfgVal.configFiles){
	cfgVal.what <- sprintf('%s (copied to Data/%s)',cfgVal.cf$var,cfgVal.cf$dest)
	if(!cfgVal.isSet(cfgVal.cf$var)){
		# config.R copies these unconditionally, so an unset name means the
		# previous Data/%s stays in place. For no policy/no override name an
		# empty spec file instead.
		cfgVal.error('%s is not set, config.R copies it to Data/%s unconditionally, name a file in %s (an empty spec file if you want no %s at all)',
								 cfgVal.cf$var,cfgVal.cf$dest,location.frida.configs,cfgVal.cf$var)
		next
	}
	cfgVal.src <- file.path(location.frida.configs,cfgVal.get(cfgVal.cf$var))
	if(!isTRUE(cfgVal.checkFile(cfgVal.src,cfgVal.what,allowEmpty=TRUE))){
		next
	}
	if(file.size(cfgVal.src)==0){
		cfgVal.warning('%s is empty, the model runs without any settings from it',cfgVal.src)
	}
	if(!isTRUE(cfgVal.fridaOK)){
		next
	}
	cfgVal.dst <- file.path(location.frida,'Data',cfgVal.cf$dest)
	if(!file.exists(cfgVal.dst)){
		cfgVal.error('%s was not copied to the frida directory, %s is missing',
								 cfgVal.get(cfgVal.cf$var),cfgVal.dst)
	} else if(!identical(unname(tools::md5sum(cfgVal.src)),unname(tools::md5sum(cfgVal.dst)))){
		cfgVal.error('%s differs from %s, the copy in config.R did not succeed, the run would use stale settings',
								 cfgVal.dst,cfgVal.src)
	}
}

# baseline parms are optional, but if specified they have to be there and be a
# single set of parameters (runMLEandParmSpace.R stops on more than one row)
if(cfgVal.has('name.baselineParmFile')&&cfgVal.isSet('name.baselineParmFile')){
	cfgVal.baselineFile <- file.path(location.frida.configs,name.baselineParmFile)
	if(isTRUE(cfgVal.checkFile(cfgVal.baselineFile,'name.baselineParmFile'))){
		cfgVal.baselineParms <- try(read.csv(cfgVal.baselineFile),silent=TRUE)
		if(inherits(cfgVal.baselineParms,'try-error')){
			cfgVal.error('name.baselineParmFile: can not be read as csv: %s',cfgVal.baselineFile)
		} else if(nrow(cfgVal.baselineParms)!=1){
			cfgVal.error('name.baselineParmFile: must contain exactly one set (line) of parameters, %s has %i',
									 cfgVal.baselineFile,nrow(cfgVal.baselineParms))
		}
	}
}

# frida info files ####
# name.frida_info is special, its expected columns depend on where it came from
cfgVal.fridaInfoFile <- file.path(location.frida.info,name.frida_info)
if(!frida_info_type%in%c('user','StellaExport','OldStyleFromBilly')){
	cfgVal.error("frida_info_type '%s' is unknown, prepareSampleParms would stop",frida_info_type)
}
if(isTRUE(cfgVal.checkFile(cfgVal.fridaInfoFile,'name.frida_info'))){
	cfgVal.checkCsvColumns(cfgVal.fridaInfoFile,'name.frida_info',c('Variable','Value','Min','Max'))
	cfgVal.fridaInfo <- try(read.csv(cfgVal.fridaInfoFile),silent=TRUE)
	if(inherits(cfgVal.fridaInfo,'try-error')){
		cfgVal.error('name.frida_info: can not be read as csv: %s',cfgVal.fridaInfoFile)
	} else {
		if(frida_info_type=='StellaExport'){
			cfgVal.checkCsvColumns(cfgVal.fridaInfoFile,'name.frida_info (stella export)',
														 c('Assumption','Calibrated.Output','Climate.Impact.Parameter',
														 	'Discrete.Outflow','Explanatory.Variable','External.Parameter',
														 	'Internal.Calibration.Parameter','No.Sensi','Output',
														 	'Parameter','Policy','Unit'))
		} else if(frida_info_type=='OldStyleFromBilly'&&ncol(cfgVal.fridaInfo)<12){
			cfgVal.error('name.frida_info: %s has %i columns, at least 12 are needed for frida_info_type %s',
									 cfgVal.fridaInfoFile,ncol(cfgVal.fridaInfo),frida_info_type)
		}
		if(all(c('Value','Min','Max')%in%colnames(cfgVal.fridaInfo))){
			cfgVal.numUsable <- sum(!is.na(cfgVal.fridaInfo$Min)&!is.na(cfgVal.fridaInfo$Max)&
																!is.na(cfgVal.fridaInfo$Value)&
																cfgVal.fridaInfo$Max>cfgVal.fridaInfo$Min&
																cfgVal.fridaInfo$Min<=cfgVal.fridaInfo$Value&
																cfgVal.fridaInfo$Value<=cfgVal.fridaInfo$Max)
			if(cfgVal.numUsable==0){
				cfgVal.error('name.frida_info: %s contains no parameter with a valid Min<=Value<=Max range, nothing would be sampled',
										 cfgVal.fridaInfoFile)
			}
		}
	}
}

# the remaining info files, allowEmpty for the ones the code guards with
# file.size()>0, generated for the ones that are written by the run itself
cfgVal.infoFiles <- list(
	list(var='name.frida_integer_parms',
			 cols=c('Variable','Value','Min','Max'),allowEmpty=FALSE,generated=FALSE),
	list(var='name.frida_external_ranges',
			 cols=c('Variable','Min','Max'),allowEmpty=FALSE,generated=FALSE),
	list(var='name.frida_extra_variables_to_export_list',
			 cols=c('FRIDA.FQN'),allowEmpty=FALSE,generated=FALSE),
	list(var='name.frida_variable_exclusion_list',
			 cols=c('Variable'),allowEmpty=TRUE,generated=FALSE),
	list(var='name.frida_parameter_exclusion_list',
			 cols=c('excludedName'),allowEmpty=TRUE,generated=TRUE))
for(cfgVal.inf in cfgVal.infoFiles){
	if(!cfgVal.isSet(cfgVal.inf$var)){
		cfgVal.error('%s is not set',cfgVal.inf$var)
		next
	}
	cfgVal.path <- file.path(location.frida.info,cfgVal.get(cfgVal.inf$var))
	if(!file.exists(cfgVal.path)&&cfgVal.inf$generated){
		cfgVal.warning('%s: %s does not exist yet, it will be created during the run',
									 cfgVal.inf$var,cfgVal.path)
		next
	}
	if(isTRUE(cfgVal.checkFile(cfgVal.path,cfgVal.inf$var,allowEmpty=cfgVal.inf$allowEmpty))){
		cfgVal.checkCsvColumns(cfgVal.path,cfgVal.inf$var,cfgVal.inf$cols)
	}
}

# variables to be exported ####
# runMLEandParmSpace.R writes the export spec from
# name.frida_extra_variables_to_export_list alone, so anything that is wanted
# from the sample runs has to be listed there. Variables that are only in the
# calibration data are used for the likelihood, but are not written to the per
# variable output files.
# Without this check a typo in subSample.TargetVars only shows up in
# runDetermineRepresentativeSamples.R, long after the expensive part of the
# run is done.
cfgVal.exportedNames <- c()
cfgVal.extraExportFile <- file.path(location.frida.info,name.frida_extra_variables_to_export_list)
if(file.exists(cfgVal.extraExportFile)&&file.size(cfgVal.extraExportFile)>0){
	cfgVal.extraExport <- try(read.csv(cfgVal.extraExportFile)$FRIDA.FQN,silent=TRUE)
	if(!inherits(cfgVal.extraExport,'try-error')&&!is.null(cfgVal.extraExport)){
		cfgVal.extraExport <- cfgVal.extraExport[nchar(cfgVal.extraExport)>4]
		cfgVal.exportedNames <- unique(cleanNames(cfgVal.extraExport))
	}
}
cfgVal.exportedNames <- cfgVal.exportedNames[nzchar(cfgVal.exportedNames)]
if(length(cfgVal.exportedNames)>0){
	cfgVal.missingTargets <- subSample.TargetVars[!subSample.TargetVars%in%cfgVal.exportedNames]
	if(length(cfgVal.missingTargets)>0){
		cfgVal.error('subSample.TargetVars %s not in the list of variables to export, add them to %s',
								 paste0("'",cfgVal.missingTargets,"'",collapse=', '),cfgVal.extraExportFile)
	}
	if(cfgVal.has('whatToPlot')&&!all(whatToPlot%in%cfgVal.exportedNames)){
		cfgVal.warning("whatToPlot '%s' is not in the list of variables to export, plotting while running will fail",
									 paste0(whatToPlot[!whatToPlot%in%cfgVal.exportedNames],collapse=', '))
	}
} else {
	cfgVal.warning('could not read the list of variables to export from %s, skipping the subSample.TargetVars check',
								 cfgVal.extraExportFile)
}

# scripts ####
# Only scripts whose absence would go unnoticed are checked here. Everything
# else (funRunFRIDA.R, funPlot.R, cleanup.R, ...) has already been sourced by
# initialise.R and config.R by the time this runs, and source() on a missing
# file is a loud error anyway.
# Note that the scripts the submit script modifies (config.R,
# runInitialiseData.R, clusterHelp.R, ...) are sourced under their expID
# prefixed per run names, so checking the templates here would say nothing
# about the run.
if(!file.exists('workerFileMergeScript.R')){
	# only used by mergePerVarFiles(parStrat=3), which starts it with system()
	cfgVal.error('workerFileMergeScript.R not found in the working directory %s, merging the per variable files with parStrat 3 would fail',
							 getwd())
}
if(!file.exists('setupTMPFS.R')){
	cfgVal.warning('setupTMPFS.R not found, frida is run from %s without a tmpfs',location.frida)
}

# plotting ####
# mismatched lengths here silently recycle and produce wrong looking plots
if(cfgVal.has('CIsToPlot')){
	for(cfgVal.v in c('CIsToPlot.lty','CIsToPlot.lwd','CIsToPlot.lcol','CIsToPlot.col')){
		if(cfgVal.has(cfgVal.v)&&length(cfgVal.get(cfgVal.v))!=length(CIsToPlot)){
			cfgVal.warning('%s has %i entries, but CIsToPlot has %i',
										 cfgVal.v,length(cfgVal.get(cfgVal.v)),length(CIsToPlot))
		}
	}
}

# report ####
cfgVal.report()

# clean up, so nothing of the validator lingers in the run environment
rm(list=ls(envir=cfgVal.env)[grepl('^cfgVal\\.',ls(envir=cfgVal.env))],envir=cfgVal.env)
