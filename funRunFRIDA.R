#
# Functions to run frida
#


# prepareSampleParms ####
prepareSampleParms <- function(excludeNames=c(),sampleParms=NULL,integerParms=NULL){
	cat('Specify sampling parameters...')
	if(is.null(sampleParms)){
		cat('reading frida_info...')
		# read in the parameters in frida that have ranges defined
		frida_info <- read.csv(file.path(location.frida.info,name.frida_info))
		if(frida_info_type == 'StellaExport'){ 
			# this file is from the new export in stella 4.11
			columnsThatAreFlags <- c('Assumption', 'Calibrated.Output', 'Climate.Impact.Parameter', 
															 'Discrete.Outflow', 'Explanatory.Variable', 'External.Parameter', 
															 'Internal.Calibration.Parameter', 'No.Sensi', 'Output', 
															 'Parameter', 'Policy', 'Unit')
			# skip lines that are not parameters with ranges
			frida_info <- frida_info[!is.na(frida_info$Min)&!is.na(frida_info$Max),]
			# write zeroes for NAs in columns that are flags
			temp <- unlist(frida_info[,columnsThatAreFlags])
			temp[is.na(temp)] <- 0
			frida_info[,columnsThatAreFlags] <- temp
		} else if (frida_info_type == 'OldStyleFromBilly'){ 
			# this is a frida_info file proided by billy pre stella 4.11
			columnsThatAreFlags <- c(2,3,4,5,6,7,8,9,10,11,12)
		}	else if (frida_info_type =='user' && sum(c('Variable','Value','Min','Max')%in%colnames(frida_info))==4){ 
			# this is a simple parm file with Variable, Min, Max
			columnsThatAreFlags <- NULL
		} else {
			stop("Unkown frida_info type check config. 
If using a user supplied frida_info ensure the columns 'Variable','Value','Min','Max' are present.\n")
		}
		# select the parameters to be sampled
		if(!is.null(columnsThatAreFlags)){
			sampleParms <- frida_info[rowSums(frida_info[,columnsThatAreFlags])>0 &
																	frida_info$No.Sensi==0 &
																	frida_info$Policy==0 &
																	frida_info$Unit==0,
																-which(!colnames(frida_info)%in%c('Variable','Value','Min','Max'))]
		} else {
			sampleParms <- frida_info[,c('Variable','Value','Min','Max')]
		}
		invalidLines <- which(!((sampleParms$Max-sampleParms$Min)>0 &
															sampleParms$Min <= sampleParms$Value &
															sampleParms$Value <= sampleParms$Max))
		suppressWarnings(file.remove(file.path(location.output,'frida_info_errorCases.csv')))
		if(length(invalidLines)>0){
			cat('invalid lines detected, see frida_info_errorCases.csv...')
		}
		write.csv(sampleParms[invalidLines,],file.path(location.output,'frida_info_errorCases.csv'))
	} else {
		invalidLines <- c()
	}
	# deal with manually excluded items
	if(!redoAllCalc &&
		 file.exists(file.path(location.frida.info,name.frida_parameter_exclusion_list)) &&
		 file.size(file.path(location.frida.info,name.frida_parameter_exclusion_list))>0){
		manParExclusionList <- read.csv(file.path(location.frida.info,name.frida_parameter_exclusion_list))
		excludeNames <- unique(c(excludeNames,manParExclusionList$excludedName))
	}
	# deal with excluded Names
	excludedIdc <- which(sampleParms$Variable %in% excludeNames)
	if(length(c(invalidLines,excludedIdc))>0){
		cat('excluding invalid and excluded parms...')
		excludedIdc <- unique(c(invalidLines,excludedIdc))
		sampleParms <- sampleParms[-excludedIdc,]
	}
	# add the integer vars e.g. climate case
	sampleParms$isInteger <- rep(FALSE,nrow(sampleParms))
	if(!is.null(integerParms) && nrow(integerParms)>0){
		cat('adding integer parms...')
		for(p.i in 1:nrow(integerParms)){
			if(!integerParms$Variable[p.i]%in%excludeNames &&
				 !integerParms$Variable[p.i]%in%sampleParms$Variable){
				newIdx <- nrow(sampleParms)+1
				sampleParms[newIdx,c('Variable')] <- integerParms$Variable[p.i] # e.g. 'Climate Units.selected climate case'
				sampleParms[newIdx,c('Value','Min','Max')] <- 
					c(integerParms$Value[p.i],integerParms$Min[p.i],integerParms$Max[p.i])
				sampleParms[newIdx,c('isInteger')] <- TRUE
				sampleParmsRowNames <- rownames(sampleParms)
				if(nrow(sampleParms)>1){
					sampleParmsRowNames[newIdx] <- as.character(as.numeric(sampleParmsRowNames[newIdx-1])+1)
				} else {
					sampleParmsRowNames[newIdx] <- 1
				}
				rownames(sampleParms) <- sampleParmsRowNames
			}
		}
	}
	cat('done\n')
	return(sampleParms)
}

# write firda export vars ####
writeFRIDAExportSpec <- function(varsForExport.fridaNames,location.frida){
	varsForExport.cleanNames <- cleanNames(varsForExport.fridaNames)
	dupe.lst <- split(seq_along(varsForExport.cleanNames), varsForExport.cleanNames)
	nonDupeIdc <- c()
	for(i in 1:length(dupe.lst)){
		nonDupeIdc[i] <- dupe.lst[[i]][1]
	}
	# nonDupeIdc <- sort(nonDupeIdc)
	varsForExport.fridaNames <- varsForExport.fridaNames[nonDupeIdc]
	sink(file=file.path(location.frida,'Data',name.fridaExportVarsFile))
	cat(paste0(varsForExport.fridaNames,collapse='\n'))
	sink()
}

# write frida input ####
# uses location.frida and name.fridaInputFile from the global env.
writeFRIDAInput <- function(variables,values,policyMode=F){
	if(policyMode){
		sink(file.path(location.frida,'Data',name.fridaInputFile))
		for(domID in names(values)){
			if(!is.na(values[domID])){
				dplID <- values[domID]
				sdmID <- jointPolicies$sdmID[jointPolicies$domID==domID & 
																		 	jointPolicies$dplID==dplID]
				pfID <- pdpMeta$pfID[pdpMeta$domID==domID & pdpMeta$sdmID==sdmID]
				sdpID <- jointPolicies$sdpID[jointPolicies$domID==domID & jointPolicies$dplID==dplID]
				pf <- pdp.lst[[pfID]]
				cat(paste0(pf[which(pf$polID==sdpID),c('policyString')],'\n',collapse = ''))
			}
		}
		sink()
	} else {
		parmValues <- data.frame(Variable=unlist(variables),Value=unlist(values))
		write.table(parmValues,file = file.path(location.frida,'Data',name.fridaInputFile),
								row.names = F,col.names = F,sep=',')
	}
}

# check for free space
disk.free <- function(path = getwd()) {
	if(length(system("which df", intern = TRUE))) {
		cmd <- sprintf("df %s", path)
		exec <- system(cmd, intern = TRUE)
		exec <- strsplit(exec[length(exec)], "[ ]+")[[1]]
		exec <- as.numeric(exec[4])
		structure(exec, names = c("available"))
	} else {
		stop("'df' command not found")
	}
}

# funGitInfo ####
# The state of a git checkout: which commit it is on and whether it has been
# modified since. Used both for the frida model checkout and for the checkout of
# these analysis scripts themselves.
# Returns a named character vector. commit is 'noGit' when the directory is not
# a checkout, so callers always have something to report.
funGitInfo <- function(location.git){
	info <- c(commit='noGit',branch=NA,date=NA,author=NA,subject=NA,origin=NA,
						dirty=NA,dirtyFiles=NA)
	if(is.null(location.git)||!file.exists(file.path(location.git,'.git'))){
		return(info)
	}
	gitOut <- function(args,all=FALSE){
		out <- suppressWarnings(
			try(system(sprintf('git -C "%s" %s 2>/dev/null',location.git,args),
								 intern=TRUE),silent=TRUE))
		if(inherits(out,'try-error')||length(out)<1||is.na(out[1])||!nzchar(out[1])){
			return(if(all){character(0)}else{NA})
		}
		if(all){out}else{out[1]}
	}
	commit <- gitOut('rev-parse HEAD')
	if(is.na(commit)||!grepl('^[0-9a-f]+$',commit)){
		return(info)
	}
	info['commit'] <- commit
	info['branch'] <- gitOut('rev-parse --abbrev-ref HEAD')
	info['date'] <- gitOut('log -1 --format=%cd --date=short')
	info['author'] <- gitOut('log -1 --format=%an')
	info['subject'] <- gitOut('log -1 --format=%s')
	info['origin'] <- gitOut('config --get remote.origin.url')
	# uncommitted changes mean the commit above does not fully describe what ran
	status <- gitOut('status --porcelain',all=TRUE)
	info['dirty'] <- as.character(length(status))
	if(length(status)>0){
		names <- sub('^...','',status)
		if(length(names)>10){
			names <- c(names[1:10],sprintf('and %i more',length(names)-10))
		}
		info['dirtyFiles'] <- paste(names,collapse=', ')
	}
	return(info)
}

# funFridaVersionInfo ####
# Version of the frida model, read from the git checkout the model files were
# rsynced from. FRIDAforUncertaintyAnalysis itself has no .git, it is an rsync
# copy made with --exclude=".*".
funFridaVersionInfo <- function(location.frida.git){
	return(funGitInfo(location.frida.git))
}

# funParseConfigAssignments ####
# The top level name <- value assignments of a config file, as text, keyed by
# name. Reading the file rather than sourcing it is deliberate: sourcing
# config.R copies files into the frida directory, creates the output folder and
# writes the run metadata file, none of which may happen just to inspect it.
funParseConfigAssignments <- function(lines){
	lines <- sub('#.*$','',lines)
	m <- regmatches(lines,regexec('^([A-Za-z.][A-Za-z0-9._]*)[ \t]*<-[ \t]*(.*[^ \t])[ \t]*$',lines))
	assignments <- c()
	for(mm in m){
		if(length(mm)==3){
			assignments[mm[2]] <- mm[3]
		}
	}
	return(assignments)
}

# funConfigDiffToDefault ####
# Which config settings this run used that the default config does not.
# The default is the committed config.R, which is the template the submit script
# copies and seds per experiment, falling back to the config.R on disk when
# there is no git to read it out of.
# Returns a data frame of name, default and used, or an empty one when the run
# used the default unchanged.
funConfigDiffToDefault <- function(configFile='config.R',location.git='.',
																	 defaultRef='HEAD:config.R',defaultFile='config.R'){
	res <- data.frame(name=character(0),default=character(0),used=character(0))
	if(is.null(configFile)||!file.exists(configFile)){
		return(res)
	}
	defaultLines <- suppressWarnings(
		try(system(sprintf('git -C "%s" show %s 2>/dev/null',location.git,defaultRef),
							 intern=TRUE),silent=TRUE))
	if(inherits(defaultLines,'try-error')||length(defaultLines)==0){
		if(!file.exists(defaultFile)){
			return(res)
		}
		defaultLines <- readLines(defaultFile,warn=FALSE)
	}
	default <- funParseConfigAssignments(defaultLines)
	used <- funParseConfigAssignments(readLines(configFile,warn=FALSE))
	# bookkeeping the submit script rewrites, not a setting anybody chose
	names.all <- setdiff(union(names(default),names(used)),'name.configFile')
	for(n in names.all){
		d <- if(n%in%names(default)){default[[n]]}else{NA}
		u <- if(n%in%names(used)){used[[n]]}else{NA}
		if(!identical(d,u)){
			res[nrow(res)+1,] <- list(n,d,u)
		}
	}
	return(res)
}

# funFridaFilesChecksum ####
# Checksums the model files as they are actually run: FRIDA.stmx, the module
# .itmx files and the input data. The files this analysis writes into the Data
# directory are excluded, they say nothing about the model version and change
# from run to run.
# Returns the aggregate checksum and the md5 of each file, in path order.
funFridaFilesChecksum <- function(location.frida,exclude=c()){
	files <- list.files(location.frida,recursive=TRUE)
	files <- sort(files[!basename(files)%in%exclude])
	if(length(files)==0){
		return(list(checksum=NA,files=character(0)))
	}
	md5s <- unname(tools::md5sum(file.path(location.frida,files)))
	names(md5s) <- files
	# base R has no hash for strings, so the aggregate is the md5 of the per file
	# lines written out to a temporary file
	aggFile <- tempfile()
	writeLines(paste0(md5s,'  ',files),aggFile)
	checksum <- unname(tools::md5sum(aggFile))
	unlink(aggFile)
	return(list(checksum=checksum,files=md5s))
}

# funFridaCheckoutDiff ####
# Compares the model files that are actually run against the git checkout they
# came from, to catch a checkout that has moved on since the last
# uncertainity_update_frida.sh, in which case the recorded commit does not
# describe the model being run.
# Only files present in both are compared, the checkout has no Data files of ours
# and we have no .git.
# uncertainity_update_frida.sh edits FRIDA.stmx before it rsyncs the model out
# (the sed on sim_specs plus FRIDAforAnalysis.patch). Those edits must not count
# as a mismatch, so if FRIDA.stmx is the only file that differs and git reports
# it as unmodified in the checkout, the difference is exactly the update script's
# doing and the verdict stays clean.
# verdict is one of
#   inSync                   nothing differs
#   updateScriptChangesOnly  only FRIDA.stmx differs, and only by those edits
#   differs                  a real mismatch, name the files
#   unknown                  nothing to compare against
funFridaCheckoutDiff <- function(location.frida,location.frida.git,exclude=c()){
	res <- list(compared=character(0),differing=character(0),
							stmxPatched=NA,verdict='unknown')
	if(is.null(location.frida.git)||!dir.exists(location.frida.git)){
		return(res)
	}
	files <- list.files(location.frida,recursive=TRUE)
	files <- sort(files[!basename(files)%in%exclude])
	files <- files[file.exists(file.path(location.frida.git,files))]
	if(length(files)==0){
		return(res)
	}
	res$compared <- files
	same <- unname(tools::md5sum(file.path(location.frida,files)))==
		unname(tools::md5sum(file.path(location.frida.git,files)))
	# unreadable files give NA, treat those as differing rather than as equal
	res$differing <- files[!same%in%TRUE]
	# has the update script applied its changes to the stmx in the checkout?
	if(file.exists(file.path(location.frida.git,'.git'))){
		status <- suppressWarnings(
			try(system(sprintf('git -C "%s" status --porcelain -- FRIDA.stmx 2>/dev/null',
												 location.frida.git),intern=TRUE),silent=TRUE))
		if(!inherits(status,'try-error')){
			res$stmxPatched <- length(status)>0&&any(grepl('^.?M',status))
		}
	}
	if(length(res$differing)==0){
		res$verdict <- 'inSync'
	} else if(identical(res$differing,'FRIDA.stmx')&&isFALSE(res$stmxPatched)){
		res$verdict <- 'updateScriptChangesOnly'
	} else {
		res$verdict <- 'differs'
	}
	return(res)
}

# funFridaCheckoutDiffNote ####
# One line summary of funFridaCheckoutDiff, for the version file and the validator.
funFridaCheckoutDiffNote <- function(cmp,maxNames=5){
	if(cmp$verdict=='inSync'){
		sprintf('all %i files shared with the checkout are identical',length(cmp$compared))
	} else if(cmp$verdict=='updateScriptChangesOnly'){
		sprintf('the %i other files shared with the checkout are identical, FRIDA.stmx differs only by the changes uncertainity_update_frida.sh applies to it',
						length(cmp$compared)-1)
	} else if(cmp$verdict=='differs'){
		names <- cmp$differing
		if(length(names)>maxNames){
			names <- c(names[1:maxNames],sprintf('and %i more',length(names)-maxNames))
		}
		sprintf('%i of %i files shared with the checkout DIFFER (%s), re-run uncertainity_update_frida.sh',
						length(cmp$differing),length(cmp$compared),paste(names,collapse=', '))
	} else {
		'could not be compared to the checkout'
	}
}

# funWriteRunMetadataFile ####
# Writes a human readable record of the model version into the output folder, so
# that a result folder can still be traced back to a model version later.
# The folder names are unchanged by this, all of the information lives in the file.
# funRunMetadataJobKey ####
# What identifies the job a run metadata file belongs to. The scripts that run
# after the ensemble source the config again, sometimes from a second R session
# (the failure cleanup in the .run templates does), so a session local marker
# cannot tell "later sourcing within this run" from "a new run in this folder".
# The batch job id can, and both the real sbatch and the localSlurm stand in
# export it. Outside a job there is no id, and then the process is the run.
funRunMetadataJobKey <- function(){
	key <- Sys.getenv('SLURM_JOB_ID',unset='')
	if(!nzchar(key)){key <- paste0('pid-',Sys.getpid())}
	key
}

# funRunMetadataField ####
# The value of one 'label   value' field of a run metadata file, NA when the
# file does not carry that field.
funRunMetadataField <- function(lines,label){
	hit <- grep(sprintf('^%s +',label),lines)
	if(length(hit)==0){return(NA_character_)}
	trimws(sub(sprintf('^%s +',label),'',lines[hit[1]]))
}

# funStripRunCompletionBlock ####
# A run metadata file without the run completion summary funAppendRunCompletionSummary
# splices into it. The summary is appended after the metadata was written, so it
# has to come out again before the two can be compared, and it has to come out
# before a new one goes in.
funStripRunCompletionBlock <- function(lines){
	start <- which(lines=='run completion')
	if(length(start)==0){return(lines)}
	start <- start[1]
	ends <- which(lines=='model file checksums')
	end <- if(any(ends>start)){min(ends[ends>start])-1}else{length(lines)}
	lines[-(start:end)]
}

funWriteRunMetadataFile <- function(location.output,location.frida.git,location.frida,
																		name.output=NULL,exclude=c(),
																		fileName='runMetadata.txt',
																		configFile='config.R',location.analysis.git='.'){
	jobKey <- funRunMetadataJobKey()
	info <- funFridaVersionInfo(location.frida.git)
	analysisInfo <- funGitInfo(location.analysis.git)
	configDiff <- funConfigDiffToDefault(configFile=configFile,
																			 location.git=location.analysis.git)
	checksums <- funFridaFilesChecksum(location.frida,exclude=exclude)
	orUnknown <- function(x){
		if(is.null(x)||length(x)!=1||is.na(x)||!nzchar(as.character(x))){'unknown'}else{as.character(x)}
	}
	field <- function(label,value){sprintf('%-12s %s',label,orUnknown(value))}
	# is the model directory still what the checkout says it is?
	stmx <- file.path(location.frida,'FRIDA.stmx')
	stmx.md5 <- if(file.exists(stmx)){unname(tools::md5sum(stmx))}else{NA}
	cmp <- funFridaCheckoutDiff(location.frida,location.frida.git,exclude=exclude)
	dirtyNote <- function(gitInfo){
		if(is.na(gitInfo['dirty'])){
			'unknown'
		} else if(gitInfo['dirty']=='0'){
			'clean, the commit above describes what ran'
		} else {
			sprintf('%s file(s) UNCOMMITTED, the commit above does not fully describe what ran (%s)',
							gitInfo['dirty'],orUnknown(gitInfo['dirtyFiles']))
		}
	}
	lines <- c('FRIDA model version used for this run',
						 '=====================================',
						 '',
						 field('commit',info['commit']),
						 field('branch',info['branch']),
						 field('date',info['date']),
						 field('author',info['author']),
						 field('subject',info['subject']),
						 field('origin',info['origin']),
						 '',
						 field('model dir',location.frida),
						 field('git dir',location.frida.git),
						 sprintf('%-12s %s  (md5 over %i files, listed below)',
						 				'model files',orUnknown(checksums$checksum),length(checksums$files)),
						 field('FRIDA.stmx',stmx.md5),
						 field('vs checkout',funFridaCheckoutDiffNote(cmp)),
						 '',
						 'The model files are the checkout above after uncertainity_update_frida.sh',
						 'has applied its cleanup and FRIDAforAnalysis.patch, so they do not match',
						 'the commit byte for byte. The checksums below cover the model as it ran,',
						 'excluding the files this analysis writes into the Data directory.',
						 '',
						 field('run',name.output),
						 field('job',jobKey),
						 field('written',format(Sys.time(),'%Y-%m-%d %H:%M:%S')))
	# the analysis scripts are as much a part of what produced these results as
	# the model is. The submit script copies them to <expID>_*.R before running
	# them, so the state of this checkout is the only record of what they were.
	lines <- c(lines,'',
						 'Analysis scripts used for this run',
						 '==================================',
						 '',
						 field('commit',analysisInfo['commit']),
						 field('branch',analysisInfo['branch']),
						 field('date',analysisInfo['date']),
						 field('author',analysisInfo['author']),
						 field('subject',analysisInfo['subject']),
						 field('origin',analysisInfo['origin']),
						 field('script dir',normalizePath(location.analysis.git,mustWork=FALSE)),
						 field('working tree',dirtyNote(analysisInfo)))
	# which settings this run does not share with the default config
	lines <- c(lines,'',
						 'Config settings differing from the default',
						 '==========================================',
						 '',
						 field('config used',configFile),
						 field('default',sprintf('%s as committed',orUnknown('config.R'))),
						 '')
	if(nrow(configDiff)==0){
		lines <- c(lines,'none, this run used the default config unchanged')
	} else {
		width <- max(nchar(configDiff$name))
		lines <- c(lines,
							 sprintf('%-*s  %s',width,'setting','default -> used'),
							 sprintf('%-*s  %s -> %s',width,configDiff$name,
							 				ifelse(is.na(configDiff$default),'(not set)',configDiff$default),
							 				ifelse(is.na(configDiff$used),'(not set)',configDiff$used)))
	}
	# the run completion summary is appended here by funAppendRunCompletionSummary
	# once the ensemble has run
	lines <- c(lines,'',
						 'model file checksums',
						 '--------------------',
						 paste0(ifelse(is.na(checksums$files),strrep(' ',32),checksums$files),
						 			 '  ',names(checksums$files)))
	# Write once per job, verify every time after that. The scripts that run after
	# the ensemble source the config again, and rewriting the file there would drop
	# the run completion summary funAppendRunCompletionSummary appended to it in
	# between, which is how the summary kept going missing. So a later sourcing
	# within the same job only checks: what it would write now has to match what
	# the first sourcing wrote. If it does not, the config, the model or the
	# analysis scripts changed on disk while the run was going, and the recorded
	# metadata no longer describes every stage of it. There is no version of that
	# worth continuing with, so it is an error rather than a warning.
	file <- file.path(location.output,fileName)
	# the timestamp is the one field that legitimately differs between the write
	# and the checks that follow it
	substantive <- function(l){l[!startsWith(l,'written')]}
	previous <- if(file.exists(file)){readLines(file,warn=FALSE)}else{character(0)}
	if(length(previous)>0&&identical(funRunMetadataField(previous,'job'),jobKey)){
		recorded <- substantive(funStripRunCompletionBlock(previous))
		current <- substantive(lines)
		if(!identical(recorded,current)){
			only <- function(x,y){setdiff(x,y)}
			gone <- only(recorded,current)
			came <- only(current,recorded)
			detail <- c(if(length(gone)>0){c('  recorded when this run started:',paste0('    ',utils::head(gone,12)))},
									if(length(came)>0){c('  would be recorded now:',paste0('    ',utils::head(came,12)))})
			stop(sprintf(paste0('The run metadata in %s was written by an earlier stage of job %s,\n',
													'and sourcing the config again now would record something different.\n',
													'Something this metadata describes, the config, the model files or the\n',
													'analysis scripts, changed on disk while the run was going, so the run\n',
													'is no longer the one the metadata claims.\n%s\n'),
									 file,jobKey,paste(detail,collapse='\n')),call.=FALSE)
		}
		return(invisible(info))
	}
	writeLines(lines,file)
	invisible(info)
}

# run completion status ####
# The log likelihood doubles as a record of whether a run completed: a run that
# failed carries logLike.failedRun plus one logLike.quasiEps per year of output
# it produced (see initialise.R). That is precise but unreadable, so every run
# also reports its completion directly, in the runStatus pseudo variable.
#   completed     1 when the run produced output for the final year
#   failYear      the first year without output, NA when completed
#   likelihoodOK  1 when logLike is a real value rather than one of the markers,
#                 NA in policy mode, where there is no calibration likelihood to
#                 compute and logLike is unconditionally a marker
# completed is about the model run reaching the end, likelihoodOK about the
# likelihood being computable. Both together are what the log likelihood test
# logLike > logLike.failedRun.max used to say on its own.

# funRunReachedFinalYear ####
# Whether one FRIDA execution produced usable output for the whole model
# horizon. The same test funRunStatusOfRun applies per state of the world, over
# all of the columns at once, for the callers that only need the yes or no.
# Stella writes a short output file when a run stops early rather than a full
# length one padded with NAs, so testing the last row of runDat for an NA, as
# the callers of this used to, finds nothing wrong with a run that stopped in
# its first year. The horizon therefore has to come from outputDataYears, the
# horizon of the default run, and never from nrow(runDat), which is itself
# short for a truncated run.
funRunReachedFinalYear <- function(runDat,
																	 years=if(exists('outputDataYears')){outputDataYears}else{rownames(runDat)}){
	nrow(runDat)>=length(years)&&!any(is.na(as.matrix(runDat)))
}

# funLikelihoodOK ####
# Whether a log likelihood is a real value rather than one of the failed run
# markers. NA in policy mode, where there is no calibration likelihood to
# compute and the log likelihood is unconditionally a marker.
funLikelihoodOK <- function(logLike,policyMode=F){
	if(policyMode){
		NA_integer_
	} else {
		as.integer(length(logLike)==1&&is.finite(logLike)&&logLike>logLike.failedRun.max)
	}
}

# funRunStatusOfRun ####
# The completion status of a single FRIDA execution, one row per state of the
# world. In policy mode one execution carries numSOW states of the world at
# once, and they can stop at different times, so the status is determined per
# SOW rather than per execution. Outside policy mode there is one SOW and this
# returns a single row.
# The failure year is the first year in which any variable of that SOW is in a
# failed state, so the run is reported as having failed the moment anything in
# the model goes bad, not once the bulk of it has.
# Two things a failed state looks like, and both have to be caught: stella
# writes fewer rows than the model horizon has years, and it writes rows whose
# values are NA. So the candidates are the first NA of every column and, when
# the output is short, the first year that has no row at all, and the failure
# year is the earliest of them. Whether the run completed is that measured
# against outputDataYears, the horizon of the default run, never against
# nrow(runDat), which is itself short for a truncated run.
funRunStatusOfRun <- function(runDat,origColNames,logLike=NULL,policyMode=F,
															years=if(exists('outputDataYears')){outputDataYears}else{rownames(runDat)}){
	numYears <- length(years)
	# the state of the world each column belongs to, 1 for columns without a
	# [n] suffix
	sowOfCol <- suppressWarnings(as.integer(sub('^.*\\[(\\d+)\\].*$','\\1',origColNames)))
	sowOfCol[is.na(sowOfCol)] <- 1L
	if(length(sowOfCol)!=ncol(runDat)){
		# nothing to group by, treat the run as a single state of the world
		sowOfCol <- rep(1L,ncol(runDat))
	}
	sowIDs <- sort(unique(sowOfCol))
	# a run whose output is short has already failed in the first year it has no
	# row for, whatever its columns say
	firstFailIdx <- rep(nrow(runDat)+1L,length(sowIDs))
	for(s.i in seq_along(sowIDs)){
		cols <- which(sowOfCol==sowIDs[s.i])
		isNA <- is.na(as.matrix(runDat[,cols,drop=FALSE]))
		# the first year each column of this SOW is in a failed state
		perColFirstNA <- apply(isNA,2,function(col){
			failed <- which(col)
			if(length(failed)==0){NA_integer_}else{min(failed)}
		})
		perColFirstNA <- perColFirstNA[!is.na(perColFirstNA)]
		if(length(perColFirstNA)>0){
			firstFailIdx[s.i] <- min(firstFailIdx[s.i],min(perColFirstNA))
		}
	}
	completed <- as.integer(firstFailIdx>numYears)
	failYear <- ifelse(completed==1,NA_real_,
										 suppressWarnings(as.numeric(years[firstFailIdx])))
	likelihoodOK <- funLikelihoodOK(logLike,policyMode=policyMode)
	return(data.frame(sowID=seq_along(sowIDs),
										completed=completed,
										failYear=failYear,
										likelihoodOK=rep(likelihoodOK,length(sowIDs))))
}

# funDecodeLogLikeRunStatus ####
# Reconstructs the run status from the log likelihood markers, for output that
# predates the runStatus file and for work unit files resumed from such a run.
# The decode is the inverse of the marker arithmetic in runFridaParmsByIndex:
# the number of years of output a failed run produced is how many
# logLike.quasiEps it sits above logLike.failedRun.
# This encoding has no per SOW resolution, so the result is always one row per
# run, with an id key.
funDecodeLogLikeRunStatus <- function(logLike,ids=seq_along(logLike),
																			years=outputDataYears,policyMode=F){
	logLike <- as.numeric(logLike)
	numYears <- length(years)
	isMarker <- !is.na(logLike)&logLike<=logLike.failedRun.max
	nYears <- rep(numYears,length(logLike))
	nYears[isMarker] <- round((logLike[isMarker]-logLike.failedRun)/logLike.quasiEps)
	nYears[is.na(logLike)] <- NA
	completed <- as.integer(nYears>=numYears)
	failYear <- rep(NA_real_,length(logLike))
	incomplete <- which(completed%in%0)
	failYear[incomplete] <- suppressWarnings(as.numeric(years[nYears[incomplete]+1]))
	likelihoodOK <- if(policyMode){
		rep(NA_integer_,length(logLike))
	} else {
		as.integer(!isMarker)
	}
	likelihoodOK[is.na(logLike)] <- NA
	return(data.frame(id=ids,completed=completed,failYear=failYear,
										likelihoodOK=likelihoodOK))
}


# funRunStatusFolders ####
# The two folders the run status is written to. location.output is sometimes
# handed in already pointing at detectedParmSpace (runMLEandParmSpace.R does
# that), so collapse a doubled detectedParmSpace the same way mergePerVarFiles
# does, and derive the run folder by stripping it off again.
funRunStatusFolders <- function(location.output,baseWD=NULL){
	folder <- location.output
	if(!is.null(baseWD)&&!grepl('^/',folder)){
		folder <- file.path(baseWD,folder)
	}
	while(grepl('/detectedParmSpace/detectedParmSpace',folder)){
		folder <- gsub('/detectedParmSpace/detectedParmSpace','/detectedParmSpace',folder)
	}
	folder <- sub('/detectedParmSpace/?$','',folder)
	return(list(run=folder,detectedParmSpace=file.path(folder,'detectedParmSpace')))
}

# funWriteRunStatusPlainCsv ####
# A plain uncompressed copy of the merged run status at the top of the run
# folder, so that the completion of an ensemble can be read without unpacking
# anything. na='' is what makes failYear come out empty for completed runs.
funWriteRunStatusPlainCsv <- function(runStatus,location.output,baseWD=NULL,
																			fileName='runStatus.csv'){
	folders <- funRunStatusFolders(location.output,baseWD)
	file <- file.path(folders$run,fileName)
	write.csv(runStatus,file,row.names=FALSE,na='',quote=FALSE)
	return(invisible(file))
}

# funOutputDataYearsFromPerVarFiles ####
# The years of the model output, read off the column names of any merged per
# variable file. outputDataYears is only in scope for scripts that ran
# runInitialiseData.R or clusterHelp.R, and decoding a failure year out of the
# log likelihood markers needs the years whatever the caller sourced.
funOutputDataYearsFromPerVarFiles <- function(location.runFiles,outputType='RDS'){
	files <- list.files(location.runFiles)
	varNames <- setdiff(gsub('\\.RDS$|\\.csv\\.gz$|\\.csv$','',files),
											c('logLike','runStatus'))
	for(varName in varNames){
		varData <- try(readPerVarFile(file.path(location.runFiles,varName),
																	outputType=outputType),silent=TRUE)
		if(inherits(varData,'try-error')){next}
		years <- suppressWarnings(as.numeric(colnames(varData)))
		years <- colnames(varData)[!is.na(years)]
		if(length(years)>0){return(years)}
	}
	return(NULL)
}

# funReadRunStatus ####
# The run status of an ensemble, from the merged runStatus file when there is
# one and decoded out of the merged log likelihood markers when there is not,
# which is what output produced before the runStatus file looks like.
# With numSample given the result covers every expected id, so that a run that
# never produced any output at all shows up as NA rather than being absent.
funReadRunStatus <- function(location.output,outputType=perVarOutputTypes[1],
														 numSample=NULL,baseWD=NULL,policyMode=F){
	folders <- funRunStatusFolders(location.output,baseWD)
	location.runFiles <- file.path(folders$detectedParmSpace,paste0('PerVarFiles-',outputType))
	file.runStatus <- file.path(location.runFiles,'runStatus')
	file.logLike <- file.path(location.runFiles,'logLike')
	exts <- c('.RDS','.csv','.csv.gz')
	if(any(file.exists(paste0(file.runStatus,exts)))){
		runStatus <- readPerVarFile(file.runStatus,outputType=outputType)
	} else if(any(file.exists(paste0(file.logLike,exts)))){
		logLike.perVar <- readPerVarFile(file.logLike,outputType=outputType)
		years <- if(exists('outputDataYears')){outputDataYears}else{NULL}
		if(is.null(years)){
			years <- funOutputDataYearsFromPerVarFiles(location.runFiles,outputType=outputType)
		}
		if(is.null(years)){
			stop(sprintf(paste0('no merged runStatus file in %s, and the years of the model output\n',
													'could not be determined to decode the failure years out of the log\n',
													'likelihood markers instead\n'),location.runFiles))
		}
		runStatus <- funDecodeLogLikeRunStatus(logLike.perVar$logLike,
																					 ids=logLike.perVar[[1]],
																					 years=years,
																					 policyMode=policyMode)
	} else {
		stop(sprintf(paste0('no merged runStatus and no merged logLike file in %s\n',
												'Have you run runMLEandParmSpace, and did the per variable files get merged?\n'),
								 location.runFiles))
	}
	# the merged file comes back as doubles (coercePerVarTypes), the decoder
	# produces integers. Normalise, so that callers see the same thing either way.
	for(cn in colnames(runStatus)){
		runStatus[[cn]] <- as.double(runStatus[[cn]])
	}
	# pad out to every expected id, so missing runs are visible
	if(!is.null(numSample)&&!'polID'%in%colnames(runStatus)){
		missing <- setdiff(1:numSample,runStatus$id)
		if(length(missing)>0){
			pad <- runStatus[rep(NA_integer_,length(missing)),]
			pad$id <- missing
			runStatus <- rbind(runStatus,pad)
		}
		runStatus <- runStatus[order(runStatus$id),]
		rownames(runStatus) <- NULL
	}
	return(runStatus)
}

# funRunCompletionSummary ####
# The completion of an ensemble as lines of text, for the run metadata file and
# for the terminal.
funRunCompletionSummary <- function(runStatus,numSample=NULL,maxYears=50){
	isPolicy <- 'polID'%in%colnames(runStatus)
	numRuns <- nrow(runStatus)
	present <- !is.na(runStatus$completed)
	complete <- sum(runStatus$completed%in%1)
	incomplete <- sum(runStatus$completed%in%0)
	noLike <- sum(runStatus$completed%in%1&runStatus$likelihoodOK%in%0)
	missing <- sum(!present)
	if(!is.null(numSample)&&!isPolicy){
		missing <- missing+max(0,numSample-numRuns)
		numRuns <- max(numRuns,numSample)
	}
	pct <- function(n){if(numRuns>0){sprintf(' (%.2f%%)',100*n/numRuns)}else{''}}
	field <- function(label,n,withPct=TRUE){
		sprintf('%-16s %8i%s',label,n,if(withPct){pct(n)}else{''})
	}
	lines <- c('run completion',
						 '--------------',
						 field(if(isPolicy){'policy runs'}else{'runs'},numRuns,FALSE),
						 field('complete',complete),
						 field('incomplete',incomplete),
						 field('no likelihood',noLike),
						 field('missing',missing))
	if(isPolicy){
		lines <- c(lines,
							 sprintf('%-16s %8i','policies',length(unique(runStatus$polID))),
							 sprintf('%-16s %8i','states of world',length(unique(runStatus$sowID))),
							 '',
							 'One row per policy and state of the world. There is no calibration',
							 'likelihood in policy mode, so likelihoodOK is not applicable.')
	}
	failYears <- runStatus$failYear[!is.na(runStatus$failYear)]
	if(length(failYears)>0){
		tab <- table(failYears)
		years <- as.numeric(names(tab))
		ord <- order(years)
		tab <- tab[ord]
		years <- years[ord]
		shown <- min(length(tab),maxYears)
		lines <- c(lines,'',
							 sprintf('%-16s %8s','failure year','runs'),
							 sprintf('%-16.0f %8i',years[1:shown],as.integer(tab)[1:shown]))
		if(length(tab)>shown){
			lines <- c(lines,sprintf('and %i further year(s), see runStatus.csv',
															 length(tab)-shown))
		}
	} else if(incomplete>0){
		lines <- c(lines,'','no failure year could be determined for the incomplete runs')
	}
	return(lines)
}

# funAppendRunCompletionSummary ####
# Puts the completion summary into the run metadata file and onto the terminal.
# The metadata file is written by funWriteRunMetadataFile at config time, long
# before the ensemble has run, and by the time it has the tmpfs frida directory
# may be gone, so this appends to the file rather than rebuilding it. The block
# goes in front of the model file checksums, so that the long checksum list
# stays at the end where it does not get in the way.
funAppendRunCompletionSummary <- function(location.output,runStatus,numSample=NULL,
																					baseWD=NULL,fileName='runMetadata.txt',
																					printToTerminal=TRUE){
	lines <- funRunCompletionSummary(runStatus,numSample=numSample)
	if(printToTerminal){
		cat(paste0(paste(lines,collapse='\n'),'\n'))
	}
	folders <- funRunStatusFolders(location.output,baseWD)
	file <- file.path(folders$run,fileName)
	if(!file.exists(file)){
		warning(sprintf('no %s to write the run completion summary into',file),
						call.=FALSE,immediate.=TRUE)
		return(invisible(NULL))
	}
	# an earlier summary of the same run would otherwise accumulate
	old <- funStripRunCompletionBlock(readLines(file,warn=FALSE))
	at <- which(old=='model file checksums')
	block <- c(lines,'')
	if(length(at)>0){
		at <- at[1]
		new <- c(old[1:(at-1)],block,old[at:length(old)])
	} else {
		new <- c(old,'',block)
	}
	writeLines(new,file)
	return(invisible(file))
}


# runFridaParmsByIndex ####
# Uses from global env:
#   sampleParms,samplePoints,location.frida, and name.fridaInputFile
# If retNegLogLike also uses from global env:
# 	calDat,resSigma
runFridaParmsByIndex <- function(runid,silent=T,policyMode=F,testStellaGood=F){
	retlist <- vector(mode = "list", length = length(runid))
	cat('\n')
	for(i in runid){
		cat(paste('\r',i))
		if(i <= nrow(samplePoints)){
			sink(file.path(location.frida,'lastRun.txt'))
			cat(row.names(samplePoints)[i],'\n')
			sink()
			writeFRIDAInput(colnames(samplePoints),samplePoints[i,],policyMode=policyMode)
			# see runFRIDASpecParms, same reasoning. Here a stale read would not even fail,
			# it would score this sample point with the results of the previous one.
			outputFile <- file.path(location.frida,'Data',name.fridaOutputFile)
			outputFile.mtimeBefore <- file.mtime(outputFile)
			stella_simulator_exit_status <- system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',
									 file.path(location.frida,'FRIDA.stmx')),
						 ignore.stdout = silent,ignore.stderr = silent,wait = T)
			if(testStellaGood&&stella_simulator_exit_status!=0){
				cat('Something wrong with stella simulator.\n')
				system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',#'-s', #to output isdb
										 file.path(location.frida,'FRIDA.stmx')),
							 ignore.stdout = F,ignore.stderr = F,wait = T)
				stop('Something wrong with stella simulator.\n')
			}
			outputFile.mtimeAfter <- file.mtime(outputFile)
			if(is.na(outputFile.mtimeAfter)||
				 (!is.na(outputFile.mtimeBefore)&&outputFile.mtimeAfter<=outputFile.mtimeBefore)){
				stop(sprintf(paste0('The stella simulator wrote no output for sample point %s.\n',
														'  exit status: %i\n',
														'  output file: %s\n',
														'  %s\n',
														'Reading it would have returned the result of an earlier run.\n'),
										 row.names(samplePoints)[i],stella_simulator_exit_status,outputFile,
										 ifelse(is.na(outputFile.mtimeAfter),'does not exist',
										 			 sprintf('unchanged since %s',
										 			 				format(outputFile.mtimeAfter,'%Y-%m-%d %H:%M:%OS3')))))
			}
			runDat <- read.csv(outputFile)
			origColNames <- unname(unlist(read.table(outputFile,
															 sep=',')[1,]))[-1]
			colnames(runDat) <- cleanNames(colnames(runDat))
			# catch failed runs causing NAs in year variable and crash in the rownames assignment
			runDat$year <- seq(runDat$year[1],length.out=nrow(runDat))
			rownames(runDat) <- runDat$year
			runDat <- runDat[,-1]
			# whether the run completed decides whether it gets a real log likelihood
			# or one of the failed run markers, so the status is determined first and
			# the likelihood follows from it. The completion cannot be read off the
			# last row of runDat: a run that stopped early leaves a short output file,
			# whose last row is a perfectly good year.
			runStatus <- funRunStatusOfRun(runDat,origColNames,policyMode=policyMode)
			runComplete <- all(runStatus$completed%in%1)
			yearsOfOutput <- sum(!is.na(runDat[[1]]))
			if(!policyMode){
				calDatInRunDat <- which(colnames(calDat)%in%colnames(runDat))
				if(length(calDatInRunDat)>0){
					resDat <- calDat[calDatInRunDat]-runDat[1:nrow(calDat),colnames(calDat)[calDatInRunDat]]
					logLike <- funLogLikelihood(resDat,resSigma)
				} else {
					logLike <- rep(1,ncol(runDat))
				}
				# a run that did not complete, or whose likelihood is not a single real
				# value, carries the marker instead. We use this when narrowing the
				# parms space
				if(!runComplete||!(length(logLike)==1&&is.finite(logLike))){
					logLike <- logLike.failedRun+yearsOfOutput*logLike.quasiEps
				}
			} else {
				logLike <- logLike.failedRun+yearsOfOutput*logLike.quasiEps
			}
			runStatus$likelihoodOK <- funLikelihoodOK(logLike,policyMode=policyMode)
			suppressWarnings(parmsIndex<-as.numeric(row.names(samplePoints)[i]))
			if(is.na(parmsIndex)){
				retlist[[i]] <- (list(parmsName=row.names(samplePoints)[i],
															parmsIndex=i,
															runDat=runDat,
															origColNames=origColNames,
															logLike=logLike,
															runStatus=runStatus))
			} else {
				retlist[[i]] <- (list(parmsIndex=parmsIndex,
															runDat=runDat,
															origColNames=origColNames,
															logLike=logLike,
															runStatus=runStatus))
			}
		}
	}
	return(retlist)
}
# runFridaParmsBySamplePoints ####
# the same as above, but runs all samples in samplePoints for pre allocated
# samplePoints per worker.
runFridaParmsBySamplePoints <- function(policyMode=F){
	retlist <- runFridaParmsByIndex(1:nrow(samplePoints),policyMode=policyMode)
	if(writePerWorkerFiles){
		workerID <- ifelse(exists('workerID'),workerID,0)
		workUnit.i <- ifelse(exists('workUnit.i'),workUnit.i,0)
		perVarRet <- saveParOutputToPerVarFiles(parOutput = retlist,workUnit.i = workUnit.i,
														 workerID = workerID,policyMode = policyMode)
		if(doNotReturnRunDataSavePerWorkerOnly){
			newRetlist <- list()
			for(i in 1:length(retlist)){
				newRetlist[[i]] <- list(parmsIndex=retlist[[i]]$parmsIndex,
																logLike=retlist[[i]]$logLike,
																runStatus=retlist[[i]]$runStatus)
			}
			retlist <- newRetlist
		}
		retlist[['logLike.df']] <- perVarRet$logLike
		retlist[['runStatus.df']] <- perVarRet$runStatus
	}
	return(retlist)
}

# runFridaDefaultParms ####
# Uses location.frida, and name.fridaInputFile
# from the global environment
runFridaDefaultParms <- function(silent=T,testStellaGood=F){
	return(runFRIDASpecParms(c(),silent=silent,testStellaGood = testStellaGood))
}

# runFRIDASpecParms ####
runFRIDASpecParms <- function(parVect,silent=T,testStellaGood=F){
	# Wall clock in the parscale and range determinations is the number of stella
	# runs divided by the number of workers, so that count is the thing to measure
	# when either is being made faster. Inert unless developmentTools has switched
	# it on; see developmentTools/countModelRuns.R, which is where everything else
	# about it lives.
	if(exists('devTools.countModelRuns')&&isTRUE(devTools.countModelRuns)){
		devTools.modelRunCount <<- if(exists('devTools.modelRunCount')){
			devTools.modelRunCount+1
		} else {
			1
		}
	}
	if(is.null(names(parVect))&length(parVect)>0){
		stop('need names in parVect to write FRIDA input\n')
	}
	if(length(parVect)==0){
		writeFRIDAInput(c(),c())
	}else{
		writeFRIDAInput(names(parVect),parVect)
	}
	# the simulator overwrites the output file, so remember the state of the one that is
	# there now. Whatever we read below has to be newer than this, otherwise it is the
	# result of some earlier run and not of this one.
	outputFile <- file.path(location.frida,'Data',name.fridaOutputFile)
	outputFile.mtimeBefore <- file.mtime(outputFile)
	stella_simulator_exit_status <- system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',#'-s', #to output isdb
							 file.path(location.frida,'FRIDA.stmx')),
				 ignore.stdout = silent,ignore.stderr = silent,wait = T)
	if(testStellaGood&&stella_simulator_exit_status!=0){
		cat('Something wrong with stella simulator.\n')
		system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',#'-s', #to output isdb
								 file.path(location.frida,'FRIDA.stmx')),
					 ignore.stdout = F,ignore.stderr = F,wait = T)
		stop('Something wrong with stella simulator.\n')
	}
	# A run that leaves the output file untouched has not produced the result we are
	# about to read. Without this check read.csv silently returns whatever the last run
	# wrote, and the run before that may well have exported a different set of
	# variables. That surfaces much later, as a mismatch between the columns of the
	# calibration data and of the model result data, far away from what caused it.
	# Runs that do not complete are a normal outcome when the optimiser tries extreme
	# parameters. Those still write output, with NAs from the point they stopped on, and
	# the callers detect them by that. So a non zero exit status is only an error when
	# the output file is stale as well, and then it is the reason to report.
	outputFile.mtimeAfter <- file.mtime(outputFile)
	if(is.na(outputFile.mtimeAfter)||
		 (!is.na(outputFile.mtimeBefore)&&outputFile.mtimeAfter<=outputFile.mtimeBefore)){
		stop(sprintf(paste0('The stella simulator wrote no output.\n',
												'  exit status: %i\n',
												'  output file: %s\n',
												'  %s\n',
												'Reading it would have returned the result of an earlier run.\n'),
								 stella_simulator_exit_status,outputFile,
								 ifelse(is.na(outputFile.mtimeAfter),'does not exist',
								 			 sprintf('unchanged since %s',
								 			 				format(outputFile.mtimeAfter,'%Y-%m-%d %H:%M:%OS3')))))
	}
	runDat <- read.csv(outputFile)
	colnames(runDat) <- cleanNames(colnames(runDat))
	if('year' %in% colnames(runDat) &&
		 sum(is.na(runDat$year))<nrow(runDat)){
		runDat <- runDat[!is.na(runDat$year),]
	}
	rownames(runDat) <- runDat$year
	runDat <- runDat[,-1]
	return(runDat)
}

# cleanNames ####
# takes a vector of e.g. column names and brings them into 
# a comparable standard format
# also drops the trailing 1 of the run id which we do not use
# as we only have single run setups of frida
cleanNames <- function(colNames){
	gsub('time','year',
			 gsub('_+$','',
			 		 gsub('_+','_',
			 		 		 gsub(',','_',
			 		 		 		 gsub('\\$','',
			 		 		 		 		 gsub('_1','',
			 		 		 		 		 		 gsub('\\]','_',
			 		 		 		 		 		 		 gsub('\\[\\*','_',
				 		 		 		 		 		 		 gsub('\\[\\d+','_',
				 		 		 		 		 		 		 		 gsub('[. ]','_',
				 		 		 		 		 		 		 		 		 tolower(colNames)))))))))))
}

# idxOfVarName ####
idxOfVarName <- function(varNames,vecOfVarNames){
	varNames <- cleanNames(varNames)
	vecOfVarNames <- cleanNames(vecOfVarNames)
	idx <- c()
	for(i in 1:length(varNames)){
		idx[i] <- which(vecOfVarNames==varNames[i])
	}
	return(idx)
}

# taken from funchir package
# convert plot region sizes to inches
xdev2in <- function (x = 1) 
{
	x * par("pin")[1L]/diff(par("usr")[1L:2L])
}
ydev2in <- function (y = 1) 
{
	y * par("pin")[2L]/diff(par("usr")[3L:4L])
}
xydev2in <-	function (xy = 1) 
{
	u = par("usr")
	xy * par("pin")/c(u[2L] - u[1L], u[4L] - u[3L])
}
dist.f <- function(oi,yi,ys,offsets,keepOutSize){
	return(min(sqrt((ydev2in(ys-yi))^2+(xdev2in(offsets-oi))^2))-keepOutSize)
}

# funValidRange ####
# returns the first and last index in the variable that has a data point
funValidRange <- function(x){
	validRange <- c(1,length(x))
	if(is.na(x[1])){
		validRange[1] <- which(diff(is.na(x))==-1)[1]+1
	}
	if(is.na(x[length(x)])){
		validRange[2] <- max(which(diff(is.na(x))==1))
	}
	return(validRange)
}

# funLogLikelihood ####
# dmvnorm function from the mvtnorm package for reference
# dmvnorm <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE, 
# 										 checkSymmetry = TRUE) 
# {
# 	if (is.vector(x)) 
# 		x <- matrix(x, ncol = length(x))
# 	p <- ncol(x)
# 	if (!missing(mean)) {
# 		if (!is.null(dim(mean))) 
# 			dim(mean) <- NULL
# 		if (length(mean) != p) 
# 			stop("x and mean have non-conforming size")
# 	}
# 	if (!missing(sigma)) {
# 		if (p != ncol(sigma)) 
# 			stop("x and sigma have non-conforming size")
# 		if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
# 																			check.attributes = FALSE)) 
# 			stop("sigma must be a symmetric matrix")
# 	}
# 	dec <- tryCatch(base::chol(sigma), error = function(e) e)
# 	if (inherits(dec, "error")) {
# 		x.is.mu <- colSums(t(x) != mean) == 0
# 		logretval <- rep.int(-Inf, nrow(x))
# 		logretval[x.is.mu] <- Inf
# 	}
# 	else {
# 		tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
# 		rss <- colSums(tmp^2)
# 		logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * 
# 																												pi) - 0.5 * rss
# 	}
# 	names(logretval) <- rownames(x)
# 	if (log) 
# 		logretval
# 	else exp(logretval)
# }
funLogLikelihood <- function(resid,covmat,treatVarsAsIndep=.GlobalEnv$treatVarsAsIndep){
	if(treatVarsAsIndep){
		perVarLogLike <- rep(NA,ncol(resid))
		for(i in 1:ncol(resid)){
			perVarLogLike[i] <- sum(dnorm(resid[!is.na(resid[,i]),i], 0, sqrt(covmat[i,i]),log = T))
			# treat all non observed values as likely as the average observed one
			# perVarLogLike[i] <- perVarLogLike[i]+sum(is.na(resid[,i]))*perVarLogLike[i]/sum(!is.na(resid[,i]))
		}
		return(sum(perVarLogLike))
	} else {
		if(nrow(resid)==0||sum(is.na(resid))>0){
			return(-Inf)
		}	else {
			logLike <- sum(mvtnorm::dmvnorm(resid,rep(0,ncol(resid)),covmat,log=T,checkSymmetry = F))
			if(is.na(logLike)){
				return(-Inf)
			} else {
				return(logLike)
			}
		}
	}
}

# chunk ####
# cuts a vector into n equal parts
chunk <- function(x,n){
	split(x, cut(seq_along(x), n, labels = FALSE)) 
}

# cluster run ####
clusterRunFridaForSamplePoints <- function(samplePoints,chunkSizePerWorker,
																					 location.output,
																					 calDat,resSigma,
																					 redoAllCalc=F,
																					 plotDatWhileRunning=F,
																					 plotDatPerChunWhileRunning=F,
																					 plotPerChunk=T,
																					 yaxPad=0.2,
																					 baseLL=-29567.06,
																					 skipRunJustRead=F,
																					 numWorkers=length(cl),
																					 calDat.impExtrValue=NULL,
																					 plotMinAlpha = 0.01,
																					 plotBaseName = 'runs-'){
	cat('cluster run...\n')
	# If we are not plotting while running we can directly store from the running
	# worker processes. This is MUCH faster.
	if((plotDatWhileRunning|plotDatPerChunWhileRunning)&writePerWorkerFiles){
		cat('Note: Forcing writePerWorkerFiles to false, as plotting is enabled\n')
		writePerWorkerFiles <<- F
	}
	# prevent skipping work if parOutput exists in the global env.
	if(exists('parOutput')){rm(parOutput)}
	numSample <- nrow(samplePoints)
	logLike <- rep(NA,numSample)
	names(logLike) <- 1:numSample
	workUnitBoundaries <- seq(1,numSample,chunkSizePerWorker*numWorkers)
	# in case the chunkSize is not a perfect divisor of the numSample, add numSample as the 
	# final boundary
	if(workUnitBoundaries[length(workUnitBoundaries)]!=numSample){
		workUnitBoundaries <- c(workUnitBoundaries,numSample)
	}
	# add one to the last work unit boundary, as during running we always deduct one from the next boundary
	workUnitBoundaries[length(workUnitBoundaries)] <- numSample+1
	
	# initialise  cluster
	baseWD <- getwd()
	# ensure path is interepreted correctly in case location.output is provided absolute
	if(grepl('^/',location.output) && exists('baseWD')){
		location.output <- system(paste0('realpath --relative-to="',baseWD,'" "',location.output,'"'),intern = T)
	}
	if(!skipRunJustRead){
		clusterExport(cl,list('location.output','baseWD','sampleParms',
													'chunkSizePerWorker','runFridaParmsBySamplePoints',
													'calDat','resSigma',
													'runFridaParmsByIndex','writePerWorkerFiles',
													'treatVarsAsIndep','perVarOutputTypes',
													'doNotReturnRunDataSavePerWorkerOnly'))
	}
	# plot setup 
	if(plotDatWhileRunning & !plotPerChunk){
		if(!(plotCape['X11']|plotCape['aqua'])){
			ncols <- ncol(calDat)
			sqrtNcols <- sqrt(ncols)
			plotCols <- round(sqrtNcols)
			plotRows <- ceiling(sqrtNcols)
			png(file.path(baseWD,location.output,paste0(plotBaseName,'all.png')),
					width=5*plotCols,
					height=5*plotRows+5/4,
					unit='cm',res=150)
		}
		subPlotLocations <- funPlotDat(calDat,calDat.impExtrValue,yaxPad = yaxPad,
																	 shadowIncompleteYears=F)
	}
	# running
	cat(sprintf('  Run of %i runs split up into %i work units of size %i (%i per worker).\n',
							numSample,length(workUnitBoundaries)-1,chunkSizePerWorker*numWorkers,chunkSizePerWorker))
	chunkTimes <- c()
	completeRunsSoFar <- 0
	runStatus.all <- NULL
	i <- 0
	while(i<(length(workUnitBoundaries)-1)){
		i <- i+1
		workUnit.i <- i
		clusterExport(cl,list('workUnit.i'),envir=environment())
		# unlike logLike.df this is deliberately per iteration, so that a work unit
		# without a status cannot silently inherit the one of its predecessor
		if(exists('runStatus.df')){rm(runStatus.df)}
		if(!redoAllCalc && file.exists(file.path(baseWD,location.output,paste0('workUnit-',i,'.RDS')))){
			cat(sprintf('\r(r) Using existing unit %i',i))
			tryCatch({parOutput <- readRDS(file.path(baseWD,location.output,paste0('workUnit-',i,'.RDS')))},
							 error = function(e){},warning=function(w){})
			if(exists('parOutput')){
				lastChunkSize <- length(parOutput)
				if(!is.null(parOutput$logLike.df)){
					logLike.df <- parOutput$logLike.df
					lastChunkSize <- lastChunkSize-1
				}
				if(!is.null(parOutput$runStatus.df)){
					runStatus.df <- parOutput$runStatus.df
					lastChunkSize <- lastChunkSize-1
				}
				if(lastChunkSize>(workUnitBoundaries[i+1]-workUnitBoundaries[i])){
					cat(sprintf(', existing output has different chunkSize (%i rather than %i), resorting remaining work',
											lastChunkSize,
											chunkSizePerWorker*numWorkers))
					workUnitBoundaries[i+1] <- workUnitBoundaries[i]+lastChunkSize
					workUnitBoundaries <- c(workUnitBoundaries[1:i],
																	seq(workUnitBoundaries[i+1],numSample,chunkSizePerWorker*numWorkers))
					if(workUnitBoundaries[length(workUnitBoundaries)]!=numSample){
						workUnitBoundaries <- c(workUnitBoundaries,numSample)
					}
					workUnitBoundaries[length(workUnitBoundaries)] <- numSample+1
				}
			}
		}
		if(!exists('parOutput')){
			if(skipRunJustRead){
				stop('Missing run files.\n')
			}
			cat(sprintf('\r(r) Running unit %i: samples %i to %i. ',
									i, workUnitBoundaries[i],workUnitBoundaries[i+1]-1))
			if((workUnitBoundaries[i]-1)>0){
				cat(sprintf('So far: Complete runs %i (%.1f%%)',
										completeRunsSoFar,100*completeRunsSoFar/(workUnitBoundaries[i]-1)))
			}
			if(length(chunkTimes>1)){
				cat(sprintf(', time per unit %i s (%.2f r/s, %.2f r/s/thread), expect completion in %s',
										round(mean(chunkTimes,na.rm=T)),
										length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
										chunkSizePerWorker/mean(chunkTimes,na.rm=T),
										if(exists('dseconds')){dseconds(round(mean(chunkTimes,na.rm=T))*
																											(length(workUnitBoundaries)-i))} 
										else{paste0(round(mean(chunkTimes,na.rm=T))*(length(workUnitBoundaries)-i),'s')}))
			}
			tic()
			workUnit <- workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)
			workerWorkUnits <- chunk(workUnit,numWorkers)
			# write the samplePoints of the work units to the worker directories
			for(w.i in workers){
				if(w.i <= length(workerWorkUnits) && !is.null(workerWorkUnits[[w.i]])){
					workerSamplePoints <- array(samplePoints[workerWorkUnits[[w.i]],],
																			dim=c(length(workerWorkUnits[[w.i]]),ncol(samplePoints)))
					colnames(workerSamplePoints) <- colnames(samplePoints)
					if(ncol(samplePoints)>1){
						rownames(workerSamplePoints) <- rownames(samplePoints[workerWorkUnits[[w.i]],])
					} else {
						rownames(workerSamplePoints) <- names(samplePoints[workerWorkUnits[[w.i]],])
					}
					saveRDS(workerSamplePoints,
									file.path(name.workDir,paste0(name.workerDirBasename,w.i),'samplePoints.RDS'))
				} else {
					saveRDS(samplePoints[c(),],
									file.path(name.workDir,paste0(name.workerDirBasename,w.i),'samplePoints.RDS'))
				}
			}
			gobble <- clusterEvalQ(cl,{
				samplePoints <- readRDS('samplePoints.RDS')
			})
			parOutput <- clusterEvalQ(cl,runFridaParmsBySamplePoints())
			timing <- toc(quiet=T)
			chunkTimes[i] <- timing$toc-timing$tic
			if(writePerWorkerFiles){
				logLike.df <- data.frame(id=integer(),logLike=double())
				runStatus.df <- NULL
				for(r.i in 1:length(parOutput)){
					logLike.df <- rbind(logLike.df,parOutput[[r.i]]$logLike.df)
					runStatus.df <- rbind(runStatus.df,parOutput[[r.i]]$runStatus.df)
					parOutput[[r.i]]$logLike.df <- NULL
					parOutput[[r.i]]$runStatus.df <- NULL
				}
			}
			cat('\r(s)')
			parOutput <-  unlist(parOutput, recursive = F)
			if(writePerWorkerFiles){
				parOutput$logLike.df <- logLike.df
				parOutput$runStatus.df <- runStatus.df
			}
			saveRDS(parOutput,file.path(baseWD,location.output,paste0('workUnit-',i,'.RDS')))
			if(!writePerWorkerFiles){
				runStatus.df <- saveParOutputToPerVarFiles(parOutput=parOutput,
																									 workUnit.i=i)$runStatus
			}
			cat('\r   ')
		}
		cat('\r(l)')
		
		# if there is a full assembled logLike.df use it
		if(exists('logLike.df')){
			logLike[logLike.df$id] <- logLike.df$logLike
		} else {
			# otherwise assemble the logLikes 
			for(l in 1:length(parOutput)){
				# if the parOutput has logLike.df use it
				if('logLike.df' %in% names(parOutput[[l]])){
					logLike[parOutput[[l]]$logLike.df$id] <- parOutput[[l]]$logLike.df$logLike 
				} else {
					# otherwise collect the logLikes
					logLike[parOutput[[l]]$parmsIndex] <- parOutput[[l]]$logLike 
				}
			}
		}
		# how many runs completed now comes from the run status rather than from
		# the log likelihood markers. A work unit resumed from output that predates
		# the run status gets it decoded back out of those markers.
		if(!exists('runStatus.df')||is.null(runStatus.df)){
			unitIds <- workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)
			runStatus.df <- funDecodeLogLikeRunStatus(logLike[unitIds],ids=unitIds)
		}
		runStatus.all <- rbind(runStatus.all,runStatus.df)
		completeRunsSoFar <- sum(runStatus.all$completed==1,na.rm=T)
		cat('\r   ')
		if(plotDatWhileRunning&&!plotDatPerChunWhileRunning){
			cat('\r(p)')
			for(dat.i in 1:ncol(calDat)){
				par(mfg = which(subPlotLocations==dat.i,arr.ind = T))
				xlims <- c(min(as.numeric(rownames(calDat)))-0.5,
									 max(as.numeric(rownames(calDat)))+0.5)
				yrange <- range(calDat[[dat.i]],na.rm=T)
				ylims <- c(yrange[1]-abs(diff(yrange))*yaxPad,
									 yrange[2]+abs(diff(yrange))*yaxPad)
				plot(rownames(calDat),calDat[[dat.i]],type='n',
						 xaxt='n',yaxt='n',
						 xaxs='i',yaxs='i',
						 xlim=xlims,
						 ylim=ylims)
				for(l in 1:length(parOutput)){
					runLL <- parOutput[[l]]$logLike
					lines(rownames(parOutput[[l]]$runDat),parOutput[[l]]$runDat[[dat.i]],
								col=adjustcolor(i,min(1,max(plotMinAlpha,
																						1/(abs(runLL-baseLL)+1)
								))))
				}
			}
			cat('\r   ')
		}
		if(plotDatPerChunWhileRunning){
			cat('\r(p)')
			ncols <- ncol(calDat)
			sqrtNcols <- sqrt(ncols)
			plotCols <- round(sqrtNcols)
			plotRows <- ceiling(sqrtNcols)
			png(file.path(baseWD,location.output,paste0(plotBaseName,'-Chunk-',i,'.png')),
					width=5*plotCols,
					height=5*plotRows+5/4,
					unit='cm',res=150)
			subPlotLocations.chk <- funPlotDat(calDat,calDat.impExtrValue,yaxPad = yaxPad,
																				 shadowIncompleteYears=F)
			for(dat.i in 1:ncol(calDat)){
				par(mfg = which(subPlotLocations.chk==dat.i,arr.ind = T))
				xlims <- c(min(as.numeric(rownames(calDat)))-0.5,
									 max(as.numeric(rownames(calDat)))+0.5)
				yrange <- range(calDat[[dat.i]],na.rm=T)
				ylims <- c(yrange[1]-abs(diff(yrange))*yaxPad,
									 yrange[2]+abs(diff(yrange))*yaxPad)
				plot(rownames(calDat),calDat[[dat.i]],type='n',
						 xaxt='n',yaxt='n',
						 xaxs='i',yaxs='i',
						 xlim=xlims,
						 ylim=ylims)
				for(l in 1:length(parOutput)){
					runLL <- parOutput[[l]]$logLike
					lines(rownames(parOutput[[l]]$runDat),parOutput[[l]]$runDat[[dat.i]],
								col=adjustcolor(1,min(1,max(plotMinAlpha,
																						1/(abs(runLL-baseLL)+1)
																						))))
				}
			}
			dev.off()
			cat('\r   ')
		}
		rm(parOutput)
	}
	
	if(plotDatWhileRunning&!plotDatPerChunWhileRunning){
		cat(sprintf('  Saving figure...'))
		plotCape <- capabilities()
		if(!(plotCape['X11']|plotCape['aqua'])){
			dev.off()
		} else{
			dev.print(png,width=5*ncol(subPlotLocations),
								height=5*(nrow(subPlotLocations)-1)+5/4,
								unit='cm',res=150,
								file.path(baseWD,location.output,paste0(plotBaseName,'all.png')))
		}
		cat('done\n')
	}
	if(length(chunkTimes)==0){
		cat(sprintf('\r    all runs read, no calculation necessary. Complete runs %i (%.1f%%)                           \n',
								completeRunsSoFar,100*completeRunsSoFar/numSample))
	} else {
		cat(sprintf('\r    complete runs %i (%.2f%%), average chunk time %i sec (%.2f r/s, %.2f r/s/thread), over all run time %s %s\n',
								completeRunsSoFar,100*completeRunsSoFar/numSample,
								round(mean(chunkTimes,na.rm=T)),
								length(cl)*chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								chunkSizePerWorker/mean(chunkTimes,na.rm=T),
								dseconds(round(sum(chunkTimes,na.rm=T))),
								'                                                                             '))
	}
	# merge all the individual per Var files into one complete one
	mergePerVarFiles()
	# the completion of the ensemble, read back from the merged file so that what
	# is reported is what actually made it to disk
	runStatus <- tryCatch(funReadRunStatus(location.output,numSample=numSample,baseWD=baseWD),
												error=function(e){
													warning(sprintf('could not read the merged run status: %s',
																					conditionMessage(e)),call.=FALSE,immediate.=TRUE)
													runStatus.all
												})
	if(!is.null(runStatus)){
		funWriteRunStatusPlainCsv(runStatus,location.output,baseWD=baseWD)
		funAppendRunCompletionSummary(runStatus=runStatus,location.output=location.output,
																	numSample=numSample,baseWD=baseWD)
		attr(logLike,'runStatus') <- runStatus
	}
	return(logLike)
}

# loadClusterRuns ####
# Reads back the whole parOutput of every work unit, run data included. That
# needs the work units to have been run with doNotReturnRunDataSavePerWorkerOnly
# set to FALSE, otherwise the workUnit-<i>.RDS files hold only parameter indices
# and log likelihoods. For the per variable results use readPerVarFile instead.
loadClusterRuns <- function(location.output){
	if(grepl('^/',location.output) && exists('baseWD')){
		location.output <- system(paste0('realpath --relative-to="',baseWD,'" "',location.output,'"'),intern = T)
	}
	runFilesList <- list.files(file.path(baseWD,location.output),pattern = 'workUnit-[0-9]+\\.RDS')
	retList <- c()
	for(f.i in 1:length(runFilesList)){
		parOutput <- readRDS(file.path(baseWD,location.output,paste0('workUnit-',f.i,'.RDS')))
		retList <- c(retList,parOutput)
	}
	return(retList)
}

saveParOutputToPerVarFiles <- function(parOutput, workUnit.i='0', workerID='0',
																			 verbosity=0, policyMode=F){
	# ensure path is interepreted correctly in case location.output is provided absolute
	if(grepl('^/',location.output) && exists('baseWD')){
		location.output <- system(paste0('realpath --relative-to="',baseWD,'" "',location.output,'"'),intern = T)
	}
	varNames <- unique(parOutput[[1]]$origColNames)
	workUnitLength <- length(parOutput)
	perVarData <- list()
	logLike <- data.frame(id=rep(NA,workUnitLength),logLike=rep(NA,workUnitLength))
	varsIdc.lst <- list()
	varNamesNoSOW.all <- gsub('\\[\\d+','',varNames)
	varNamesNoSOW <- unique(varNamesNoSOW.all)
	for(v.i in 1:length(varNamesNoSOW)){
		varName <- cleanNames(varNamesNoSOW[v.i])
		varsIdc.lst[[varName]] <- which(varNamesNoSOW.all==varNamesNoSOW[v.i])
		numSOW <- length(varsIdc.lst[[cleanNames(varNamesNoSOW[v.i])]])
		if(numSOW>1){
			perVarData[[varName]] <- data.frame(matrix(NA,ncol=2+length(outputDataYears),nrow=workUnitLength*numSOW))
			colnames(perVarData[[varName]]) <- c('polID','sowID',outputDataYears)
		} else {
			perVarData[[varName]] <- data.frame(matrix(NA,ncol=1+length(outputDataYears),nrow=workUnitLength))
			colnames(perVarData[[varName]]) <- c('id',outputDataYears)
		}
	}
	varNames <- names(perVarData)
	# runStatus is keyed the same way as the data it describes: per (polID,sowID)
	# when a run carries several states of the world, per run otherwise. Its first
	# column has to be the one that ascends across chunk files, so that
	# verifyChunkOrder can take the streaming merge path, which is polID.
	# taken from the status the run itself reported where there is one, so that the
	# frame allocated here and the rows filled into it cannot disagree
	numSOW.status <- if(!is.null(parOutput[[1]]$runStatus)){
		nrow(parOutput[[1]]$runStatus)
	} else {
		length(varsIdc.lst[[1]])
	}
	if(numSOW.status>1){
		runStatus <- data.frame(polID=rep(NA,workUnitLength*numSOW.status),
														sowID=rep(NA,workUnitLength*numSOW.status),
														completed=rep(NA,workUnitLength*numSOW.status),
														failYear=rep(NA,workUnitLength*numSOW.status),
														likelihoodOK=rep(NA,workUnitLength*numSOW.status))
	} else {
		runStatus <- data.frame(id=rep(NA,workUnitLength),
														completed=rep(NA,workUnitLength),
														failYear=rep(NA,workUnitLength),
														likelihoodOK=rep(NA,workUnitLength))
	}
	if(verbosity>0){
		cat('Reading parOutput:\n')
	}
	for(run.i in 1:length(parOutput)){
		if(verbosity>0){
			cat(sprintf('\rReading run %i of %i.',run.i, workUnitLength))
		}
		runDat <- parOutput[[run.i]]$runDat
		for(varName in varNames){
			if(length(varsIdc.lst[[varName]])>1){
				perVarDataIndices <- (run.i+(run.i-1)*(numSOW-1)):((run.i+(run.i-1)*(numSOW-1))+numSOW-1)
				perVarData[[varName]][perVarDataIndices,'polID'] <- parOutput[[run.i]]$parmsIndex
				perVarData[[varName]][perVarDataIndices,'sowID'] <- 1:numSOW
				varDataFromRun <- unname(t(parOutput[[run.i]]$runDat[,varsIdc.lst[[varName]]]))
				# catches varDataFromRun being shorter than the output data frame
				perVarData[[varName]][perVarDataIndices,3:(2+ncol(varDataFromRun))] <- varDataFromRun
			} else{
				varDataFromRun <- unname(unlist(c(parOutput[[run.i]]$parmsIndex,runDat[varName])))
				# writing just into the entries for which we have data from run
				# to catch varDataFromRun being shorter than the output data frame
				# we start from one here instead 3 as above, because our ourputstructure
				# does not have the columns for polID and sowID
				perVarData[[varName]][run.i,1:length(varDataFromRun)] <- varDataFromRun
			}
		}
		logLike[run.i,] <- c(parOutput[[run.i]]$parmsIndex,parOutput[[run.i]]$logLike)
		# parOutput of a run that predates the runStatus output carries no status,
		# fall back to decoding it out of the log likelihood marker
		runStatus.run <- parOutput[[run.i]]$runStatus
		if(is.null(runStatus.run)){
			runStatus.run <- funDecodeLogLikeRunStatus(parOutput[[run.i]]$logLike,
																								 ids=1,policyMode=policyMode)
			runStatus.run$sowID <- 1
			runStatus.run <- runStatus.run[rep(1,numSOW.status),]
			runStatus.run$sowID <- 1:numSOW.status
		}
		if(numSOW.status>1){
			runStatusIndices <- ((run.i-1)*numSOW.status+1):(run.i*numSOW.status)
			runStatus[runStatusIndices,'polID'] <- parOutput[[run.i]]$parmsIndex
			runStatus[runStatusIndices,'sowID'] <- runStatus.run$sowID
			runStatus[runStatusIndices,c('completed','failYear','likelihoodOK')] <-
				runStatus.run[,c('completed','failYear','likelihoodOK')]
		} else {
			runStatus[run.i,] <- c(parOutput[[run.i]]$parmsIndex,
														 runStatus.run$completed[1],
														 runStatus.run$failYear[1],
														 runStatus.run$likelihoodOK[1])
		}
	}
	if(verbosity>0){
		cat('\nWriting to files\n')
	}
	# The chunks are always written as a single plain uncompressed csv, no matter
	# which final formats perVarOutputTypes asks for. mergePerVarFiles derives all
	# of those from these intermediates by streaming them together, which needs
	# them uncompressed and in one format only. See writePerVarChunkFile.
	chunkFolder <- file.path(baseWD,location.output,'detectedParmSpace','PerVarChunks')
	fullPrecision <- if(exists('perVarFullPrecision')){perVarFullPrecision}else{TRUE}
	for(varName in varNames){
		if(verbosity>0){
			cat(sprintf('\rWriting var %i of %i: %s %s', run.i, workUnitLength, varName, rep(' ',100)))
		}
		dir.create(file.path(chunkFolder,varName),showWarnings = F,recursive = T)
		writePerVarChunkFile(perVarData[[varName]],
												 file.path(chunkFolder,varName,
												 					paste0(varName,'-',workUnit.i,'-',workerID,'.csv')),
												 fullPrecision=fullPrecision)
	}
	if(verbosity>0){
		cat(sprintf('\rWriting logLike %s', rep(' ',100)))
	}
	dir.create(file.path(chunkFolder,'logLike'),showWarnings = F,recursive = T)
	writePerVarChunkFile(logLike,
											 file.path(chunkFolder,'logLike',
											 					paste0('logLike','-',workUnit.i,'-',workerID,'.csv')),
											 fullPrecision=fullPrecision)
	if(verbosity>0){
		cat(sprintf('\rWriting runStatus %s', rep(' ',100)))
	}
	dir.create(file.path(chunkFolder,'runStatus'),showWarnings = F,recursive = T)
	writePerVarChunkFile(runStatus,
											 file.path(chunkFolder,'runStatus',
											 					paste0('runStatus','-',workUnit.i,'-',workerID,'.csv')),
											 fullPrecision=fullPrecision)
	if(verbosity>0){
		cat('\n')
	}
	return(list(logLike=logLike,runStatus=runStatus))
}

# per variable chunk and merged file handling ####
# The workers write one plain uncompressed csv per (variable, work unit, worker)
# into
#   <location.output>/detectedParmSpace/PerVarChunks/<varName>/
# mergePerVarFiles then assembles one file per variable and per requested final
# format into
#   <location.output>/detectedParmSpace/PerVarFiles-<outputType>/<varName>.<ext>
# The merge streams the chunks into the final file rather than holding them in
# memory, so its memory use does not grow with numSample. See
# workerMergePerVarFiles.

# fwritePerVarCsv ####
# The one place per variable csv is written, so that chunks and merged files
# agree on their format down to the byte.
# fwrite writes doubles with 15 significant digits, which is not enough to
# reproduce a double exactly. fullPrecision spells them out in full instead,
# which costs about three times the writing time and 10% more bytes. It matters
# here because the merged files, the RDS included, are built from these csvs
# rather than from the doubles themselves.
# The marker for a failed run is chosen to survive either setting, see
# logLike.failedRun.
fwritePerVarCsv <- function(varData,file,fullPrecision=TRUE,compress='none'){
	if(fullPrecision){
		varData <- data.table::as.data.table(varData)
		for(j in seq_along(varData)){
			if(is.double(varData[[j]])){
				data.table::set(varData,j=j,value=sprintf('%.17g',varData[[j]]))
			}
		}
	}
	data.table::fwrite(varData,file,quote=FALSE,na='NA',nThread=1,
										 showProgress=FALSE,compress=compress)
}

# writePerVarChunkFile ####
# Writes a single chunk of one variable. Always plain uncompressed csv with an
# unquoted header, so that all chunks of a variable share a byte identical
# header line and can be concatenated without being parsed.
writePerVarChunkFile <- function(varData,file,fullPrecision=TRUE){
	fwritePerVarCsv(varData,file,fullPrecision=fullPrecision,compress='none')
}

# peekCsvBounds ####
# Cheap probe of a chunk csv. Returns its header line, the byte offset at which
# the data starts, and the id (first column) of its first and its last data row,
# without parsing the file. Returns NULL when the file is empty, has no data
# rows, or does not end on a complete line, which is what a worker that died
# mid write leaves behind.
peekCsvBounds <- function(file,tailBytes=65536L){
	size <- file.info(file)$size
	if(is.na(size) || size<=0){return(NULL)}
	con <- file(file,'rb')
	on.exit(close(con),add=TRUE)
	# header line, and with it the offset at which the data starts
	headBuf <- readBin(con,'raw',n=min(size,tailBytes))
	nlPos <- which(headBuf==as.raw(10L))
	if(length(nlPos)<2){return(NULL)} # need a header and at least one data line
	header <- sub('\r$','',rawToChar(headBuf[1:(nlPos[1]-1L)]))
	firstLine <- rawToChar(headBuf[(nlPos[1]+1L):(nlPos[2]-1L)])
	# last data line, read from the tail of the file. Widen the window in the
	# unlikely case that it does not hold two line endings.
	window <- min(size,as.numeric(tailBytes))
	repeat{
		seek(con,size-window)
		tailBuf <- readBin(con,'raw',n=window)
		nlTail <- which(tailBuf==as.raw(10L))
		if(length(nlTail)>=2 || window>=size){break}
		window <- min(size,window*8)
	}
	if(length(nlTail)<2 || nlTail[length(nlTail)]!=length(tailBuf)){
		return(NULL) # truncated, the file does not end on a complete line
	}
	lastLine <- rawToChar(tailBuf[(nlTail[length(nlTail)-1L]+1L):(nlTail[length(nlTail)]-1L)])
	firstId <- suppressWarnings(as.numeric(sub(',.*','',firstLine)))
	lastId <- suppressWarnings(as.numeric(sub(',.*','',lastLine)))
	if(is.na(firstId) || is.na(lastId)){return(NULL)}
	return(list(file=file,header=header,dataOffset=nlPos[1],
							firstId=firstId,lastId=lastId,size=size))
}

# verifyChunkOrder ####
# Decides whether the chunks of one variable can simply be concatenated in id
# order. Returns a list with
#   ok     TRUE when concatenating is safe
#   reason why it is not, when it is not
#   bounds the chunk bounds ordered by their first id
#   warn   things that are odd but do not make concatenating wrong
# Concatenating is safe when all chunks agree on their header, each chunk is
# itself ascending in its first column, and the chunks do not overlap. Gaps in
# the id coverage do not make the result wrong, only incomplete, so they are
# reported as a warning rather than rejected. That is also all the fallback
# path could do about them.
verifyChunkOrder <- function(bounds,expectedFirstId=1){
	warn <- character(0)
	if(length(bounds)==0){
		return(list(ok=FALSE,reason='no chunk files',bounds=bounds,warn=warn))
	}
	bad <- vapply(bounds,is.null,logical(1))
	if(any(bad)){
		return(list(ok=FALSE,
								reason=sprintf('%i of %i chunk files are empty or truncated',
															 sum(bad),length(bad)),
								bounds=bounds,warn=warn))
	}
	headers <- vapply(bounds,function(b){b$header},character(1))
	if(any(headers!=headers[1])){
		return(list(ok=FALSE,
								reason=sprintf('chunk files disagree on their header (%i distinct ones)',
															 length(unique(headers))),
								bounds=bounds,warn=warn))
	}
	firstIds <- vapply(bounds,function(b){b$firstId},double(1))
	lastIds <- vapply(bounds,function(b){b$lastId},double(1))
	if(any(lastIds<firstIds)){
		return(list(ok=FALSE,reason='a chunk file is not sorted by its first column',
								bounds=bounds,warn=warn))
	}
	ord <- order(firstIds)
	bounds <- bounds[ord]
	firstIds <- firstIds[ord]
	lastIds <- lastIds[ord]
	if(length(bounds)>1){
		overlap <- firstIds[-1]<=lastIds[-length(lastIds)]
		if(any(overlap)){
			return(list(ok=FALSE,
									reason=sprintf('the id ranges of %i chunk files overlap',sum(overlap)+1),
									bounds=bounds,warn=warn))
		}
		gaps <- sum(firstIds[-1]!=(lastIds[-length(lastIds)]+1))
		if(gaps>0){
			warn <- c(warn,sprintf('%i gap(s) in the id coverage of the chunk files',gaps))
		}
	}
	if(firstIds[1]!=expectedFirstId){
		warn <- c(warn,sprintf('the chunk files start at id %g rather than %g',
													 firstIds[1],expectedFirstId))
	}
	return(list(ok=TRUE,reason=NA_character_,bounds=bounds,warn=warn))
}

# coercePerVarTypes ####
# fread infers each column type from the values it happens to see, so an id
# column, or a variable that holds only whole numbers in this particular run,
# would come back as integer. The per variable data is a double matrix when the
# workers build it, so put every column back to double.
# integer64 needs its own branch. A column of whole numbers too large for an
# int32, which is what money and volume quantities look like, comes back from
# fread as bit64::integer64, and that stores a 64 bit integer in the bits of a
# double. So is.double() is TRUE for it and the plain test below would leave it
# alone, and unclassing it reinterprets those bits: 2608405563540 reads back as
# 1.288724e-311 and an NA reads back as 0, which is.na() then does not
# recognise. Converting has to go through bit64 for that reason.
# The fread calls all pass integer64='double' so this should not come up, but a
# call added without it would otherwise corrupt the data silently.
coercePerVarTypes <- function(d){
	for(cn in names(d)){
		col <- d[[cn]]
		if(inherits(col,'integer64')){
			if(!requireNamespace('bit64',quietly=TRUE)){
				stop('an integer64 column needs bit64 to be converted without corrupting it\n')
			}
			data.table::set(d,j=cn,value=bit64::as.double.integer64(col))
		} else if(!is.double(col)){
			data.table::set(d,j=cn,value=as.double(col))
		}
	}
	invisible(d)
}

# rbindChunkList ####
# Row binds the chunks of one variable, used by the paths that cannot simply
# concatenate. When all chunks agree on their number of columns the names of the
# first are imposed on all of them, which is what makes headers that differ only
# in how they were written (X1980 versus 1980) mergeable. Chunks of differing
# width are matched by name and padded instead.
rbindChunkList <- function(chunks){
	nCol <- vapply(chunks,ncol,integer(1))
	if(all(nCol==nCol[1])){
		for(k in seq_along(chunks)){
			data.table::setnames(chunks[[k]],names(chunks[[1]]))
		}
	}
	return(data.table::rbindlist(chunks,use.names=TRUE,fill=TRUE))
}

# concatChunkCsvs ####
# Streams the data rows (header line skipped) of the given chunk files into an
# open binary connection. Nothing is parsed and never more than blockSize bytes
# are held at once, so the memory used is independent of how large the variable
# is.
concatChunkCsvs <- function(bounds,outCon,blockSize=32L*1024L*1024L){
	for(b in bounds){
		con <- file(b$file,'rb')
		seek(con,b$dataOffset)
		repeat{
			buf <- readBin(con,'raw',n=blockSize)
			if(length(buf)==0){break}
			writeBin(buf,outCon)
		}
		close(con)
	}
	invisible(NULL)
}

# openCsvGzSink ####
# A binary connection that gzips into file. Prefers piping through the gzip
# binary, so that compressing happens in its own process alongside the reading
# R process, and falls back to Rs own gzfile.
openCsvGzSink <- function(file){
	gzBin <- ''
	if(nzchar(Sys.which('pigz'))){
		gzBin <- 'pigz -p 1 -c'
	} else if(nzchar(Sys.which('gzip'))){
		gzBin <- 'gzip -c'
	}
	if(nzchar(gzBin)){
		con <- try(pipe(paste0(gzBin,' > ',shQuote(file)),'wb'),silent=TRUE)
		if(!inherits(con,'try-error')){return(con)}
	}
	return(gzfile(file,'wb'))
}

# workerMergePerVarFiles ####
# Merges all chunk files of one variable into one final file per requested
# output type.
# The fast path concatenates the chunk csvs byte wise in id order, which needs
# neither parsing nor sorting nor holding the variable in memory. It is taken
# whenever verifyChunkOrder says the chunks line up, which is what the work unit
# and per worker splitting in clusterRunFridaForSamplePoints produces. Otherwise
# everything is read in, sorted and written out, which is correct but needs the
# whole variable in memory.
workerMergePerVarFiles <- function(v.i,varNames,chunkFolder,outputFolder,
																	 outputTypes=perVarOutputTypes,
																	 verbosity=0,compressCsv=T,rdsCompress=TRUE,
																	 fullPrecision=TRUE){
	varName <- varNames[v.i]
	perVarSubfolder <- file.path(chunkFolder,varName)
	fileList <- file.path(perVarSubfolder,naturalsort(list.files(perVarSubfolder)))
	if(length(fileList)==0){
		if(verbosity>0){cat(sprintf('No chunk files for %s, skipping\n',varName))}
		return(invisible(NULL))
	}
	if(verbosity>0){cat(sprintf('Processing %i files of %s...',length(fileList),varName))}
	if(length(outputTypes)==0 || !all(outputTypes %in% c('RDS','csv'))){
		stop(sprintf('outputTypes must be a non empty subset of RDS and csv, is c(%s)\n',
								 paste0("'",outputTypes,"'",collapse=',')))
	}
	outFile <- list()
	for(outputType in outputTypes){
		dir.create(file.path(outputFolder,paste0('PerVarFiles-',outputType)),
							 showWarnings = F,recursive = T)
		outFile[[outputType]] <- file.path(outputFolder,paste0('PerVarFiles-',outputType),varName)
	}
	wantCsv <- 'csv' %in% outputTypes
	wantRDS <- 'RDS' %in% outputTypes
	check <- verifyChunkOrder(lapply(fileList,peekCsvBounds))
	for(w in check$warn){
		warning(sprintf('%s: %s',varName,w),call.=FALSE,immediate.=TRUE)
	}
	if(check$ok){
		if(verbosity>0){cat('streaming...')}
		# an uncompressed csv is needed as an intermediate whenever an RDS is
		# wanted. It is gzipped afterwards, or removed when no csv was asked for.
		if(wantCsv){
			plainCsv <- paste0(outFile[['csv']],'.csv')
		} else {
			plainCsv <- file.path(outputFolder,paste0('PerVarFiles-',outputTypes[1]),
														paste0(varName,'-merging.csv'))
		}
		if(wantRDS || !compressCsv){
			outCon <- file(plainCsv,'wb')
		} else {
			outCon <- openCsvGzSink(paste0(outFile[['csv']],'.csv.gz'))
		}
		writeBin(charToRaw(paste0(check$bounds[[1]]$header,'\n')),outCon)
		concatChunkCsvs(check$bounds,outCon)
		# for a pipe this is the exit status of the compressor, a truncated output
		# would otherwise go unnoticed
		status <- close(outCon)
		if(is.numeric(status) && length(status)==1 && !is.na(status) && status!=0){
			stop(sprintf('compressing the merged csv of %s failed with status %i\n',varName,status))
		}
		if(wantRDS){
			if(verbosity>0){cat('RDS...')}
			varData <- data.table::fread(plainCsv,nThread=1,showProgress=FALSE,
																	 integer64='double')
			coercePerVarTypes(varData)
			data.table::setDF(varData)
			saveRDS(varData,paste0(outFile[['RDS']],'.RDS'),compress=rdsCompress)
			rm(varData)
		}
		if(wantCsv && compressCsv && file.exists(plainCsv)){
			if(verbosity>0){cat('gzip...')}
			status <- system2('gzip',c('-f',shQuote(plainCsv)))
			if(status!=0){
				stop(sprintf('gzip of the merged csv of %s failed with status %i\n',varName,status))
			}
		} else if(!wantCsv){
			unlink(plainCsv,force=TRUE)
		}
	} else {
		warning(sprintf('%s: cannot concatenate the chunk files (%s), reading them all instead',
										varName,check$reason),call.=FALSE,immediate.=TRUE)
		if(verbosity>0){cat(sprintf('fallback (%s)...',check$reason))}
		varData <- rbindChunkList(
			lapply(fileList,data.table::fread,nThread=1,showProgress=FALSE,
						 integer64='double'))
		data.table::setorderv(varData,names(varData)[1])
		coercePerVarTypes(varData)
		data.table::setDF(varData)
		for(outputType in outputTypes){
			writePerVarFile(varData,outFile[[outputType]],outputType=outputType,
											compressCsv=compressCsv,rdsCompress=rdsCompress,
											fullPrecision=fullPrecision)
		}
		rm(varData)
	}
	if(verbosity>0){cat('removing split files...')}
	unlink(perVarSubfolder,recursive = T,force = T)
	if(verbosity>0){cat('done\n')}
	invisible(NULL)
}

# workerMergePerVarFilesIndepProc ####
# Same as workerMergePerVarFiles but in a process of its own, so that the memory
# of the fallback path is handed back to the OS after every variable.
workerMergePerVarFilesIndepProc <- function(v.i,varNamesFileName,chunkFolder,outputFolder,
																						outputTypes,verbosity=0,compressCsv=T,
																						rdsCompress=TRUE,fullPrecision=TRUE,
																						baseWD=getwd()){
	command <- paste('Rscript',shQuote(file.path(baseWD,'workerFileMergeScript.R')),
									 v.i,shQuote(varNamesFileName),shQuote(chunkFolder),
									 shQuote(outputFolder),shQuote(paste(outputTypes,collapse=',')),
									 verbosity,compressCsv,rdsCompress,fullPrecision,shQuote(baseWD))
	if(verbosity>1){
		cat(command)
		cat('\n')
	}
	status <- system(command)
	if(status!=0){
		warning(sprintf('merge subprocess for variable %i exited with status %i',v.i,status),
						call.=FALSE)
	}
	invisible(status)
}

# workerMergePerVarFilesLegacy ####
# Merges chunk files that sit per output type in PerVarFiles-<outputType>/<varName>/
# rather than in the csv tree in PerVarChunks/. Reads them all at once, so it
# needs the whole variable in memory.
workerMergePerVarFilesLegacy <- function(v.i,varNames,outputTypeFolder,outputType,
																				 verbosity=0,compressCsv=T,rdsCompress=TRUE,
																				 fullPrecision=TRUE){
	varName <- varNames[v.i]
	perVarSubfolder <- file.path(outputTypeFolder,varName)
	fileList <- file.path(perVarSubfolder,naturalsort(list.files(perVarSubfolder)))
	if(length(fileList)==0){
		if(verbosity>0){cat(sprintf('No legacy chunk files for %s, skipping\n',varName))}
		return(invisible(NULL))
	}
	if(verbosity>0){cat(sprintf('Processing %i legacy files of %s...',length(fileList),varName))}
	varData <- rbindChunkList(
		lapply(fileList,function(f){
			d <- data.table::as.data.table(readPerVarFile(f,outputType))
			data.table::setnames(d,gsub('(^X)([0-9]{4})','\\2',names(d),perl = T))
			d
		}))
	data.table::setorderv(varData,names(varData)[1])
	coercePerVarTypes(varData)
	data.table::setDF(varData)
	writePerVarFile(varData,file.path(outputTypeFolder,varName),outputType=outputType,
									compressCsv=compressCsv,rdsCompress=rdsCompress,
									fullPrecision=fullPrecision)
	unlink(perVarSubfolder,recursive = T,force = T)
	if(verbosity>0){cat('done\n')}
	invisible(NULL)
}

# startFileMergeCluster ####
# A PSOCK cluster set up to run the merge workers. Sourcing is done from baseWD
# rather than from the inherited working directory, which is not necessarily the
# project root.
startFileMergeCluster <- function(numWorkersFileMerge,baseWD){
	clFileMerge <- makePSOCKcluster(numWorkersFileMerge)
	clusterExport(clFileMerge,'baseWD',envir=environment())
	gobble <- clusterEvalQ(clFileMerge,{
		suppressPackageStartupMessages(
			library(data.table,quietly=TRUE,warn.conflicts=FALSE))
		# we parallelise over variables, no additional OpenMP threads please
		data.table::setDTthreads(1)
		source(file.path(baseWD,'funRunFRIDA.R'))
		source(file.path(baseWD,'naturalsort.R'))
	})
	return(clFileMerge)
}

# mergePerVarFiles ####
# Assembles the per chunk files the workers wrote into one file per variable,
# for every format in outputTypes.
# parStrat 1 processes the variables one after another in this process, 2 (the
# default) spreads them over a PSOCK cluster of numWorkersFileMerge workers and
# 3 gives every variable its own R process.
mergePerVarFiles <- function(verbosity=1,parStrat=2,compressCsv=T,
														 outputTypeFoldersOverride = NULL,
														 varNamesOverride = NULL,
														 outputTypes=perVarOutputTypes,
														 rdsCompress=if(exists('perVarRdsCompress')){perVarRdsCompress}else{TRUE},
														 fullPrecision=if(exists('perVarFullPrecision')){perVarFullPrecision}else{TRUE}){
	if(verbosity>0){
		cat('Merging per Var files\n')
	}
	# ensure path is interepreted correctly in case location.output is provided absolute
	if(grepl('^/',location.output) && exists('baseWD')){
		location.output <- system(paste0('realpath --relative-to="',baseWD,'" "',location.output,'"'),intern = T)
	}
	baseWD <- if(exists('baseWD')){baseWD}else{getwd()}
	outputFolder <- file.path(baseWD,location.output,'detectedParmSpace')
	# ensure detectedParmSpace is not duplicated
	while(grepl('/detectedParmSpace/detectedParmSpace',outputFolder)){
		outputFolder <- gsub('/detectedParmSpace/detectedParmSpace','/detectedParmSpace',outputFolder)
	}
	chunkFolder <- file.path(outputFolder,'PerVarChunks')
	if(!dir.exists(chunkFolder)){
		if(verbosity>0){
			cat(sprintf('No %s, looking for chunk folders of the previous layout\n',chunkFolder))
		}
		return(invisible(mergePerVarFilesLegacy(outputFolder=outputFolder,baseWD=baseWD,
																						verbosity=verbosity,parStrat=parStrat,
																						compressCsv=compressCsv,
																						rdsCompress=rdsCompress,
																						fullPrecision=fullPrecision,
																						outputTypeFoldersOverride=outputTypeFoldersOverride,
																						varNamesOverride=varNamesOverride)))
	}
	if(!is.null(varNamesOverride)){
		varNames <- varNamesOverride
	} else {
		varNames <- basename(list.dirs(chunkFolder,recursive = F))
	}
	if(verbosity>0){
		cat(sprintf('Entering %s\n',chunkFolder))
		cat(sprintf('Found %i variable sub folder(s), writing %s\n',
								length(varNames),paste(outputTypes,collapse=' and ')))
	}
	if(length(varNames)==0){
		return(invisible(NULL))
	}
	if(parStrat==1){
		for(v.i in 1:length(varNames)){
			if(verbosity>0){cat(sprintf('(%i of %i) ',v.i,length(varNames)))}
			workerMergePerVarFiles(v.i,varNames=varNames,chunkFolder=chunkFolder,
														 outputFolder=outputFolder,outputTypes=outputTypes,
														 verbosity=verbosity,compressCsv=compressCsv,
														 rdsCompress=rdsCompress,
														 fullPrecision=fullPrecision)
		}
	} else if(parStrat==2){
		if(verbosity>0){cat(sprintf('Parallel proccessing all vars with %i workers\n',numWorkersFileMerge))}
		clFileMerge <- startFileMergeCluster(numWorkersFileMerge,baseWD)
		gobble <- parLapplyLB(clFileMerge,1:length(varNames),workerMergePerVarFiles,
													varNames=varNames,
													chunkFolder=chunkFolder,
													outputFolder=outputFolder,
													outputTypes=outputTypes,
													compressCsv=compressCsv,
													rdsCompress=rdsCompress,
													fullPrecision=fullPrecision,
													chunk.size = 1)
		stopCluster(clFileMerge)
	} else if(parStrat==3){
		if(verbosity>0){cat(sprintf('Parallel proccessing all vars with %i workers in independent processes\n',numWorkersFileMerge))}
		varNamesFileName <- file.path(baseWD,paste0('tempVarNamesListForFileMerge',
																								gsub('\\.','-',format(Sys.time(), "%Y-%m-%d-%H-%M-%OS3")),
																								'.RDS'))
		saveRDS(varNames,varNamesFileName)
		clFileMerge <- makePSOCKcluster(numWorkersFileMerge)
		gobble <- parLapplyLB(clFileMerge,1:length(varNames),workerMergePerVarFilesIndepProc,
													varNamesFileName=varNamesFileName,
													chunkFolder=chunkFolder,
													outputFolder=outputFolder,
													outputTypes=outputTypes,
													verbosity=verbosity,
													compressCsv=compressCsv,
													rdsCompress=rdsCompress,
													fullPrecision=fullPrecision,
													baseWD=baseWD,
													chunk.size = 1)
		stopCluster(clFileMerge)
		unlink(varNamesFileName,force = T)
	} else {
		stop('unkown parStrat\n')
	}
	# every worker removes its own variable sub folder, so this should be empty
	if(length(list.files(chunkFolder,all.files = T,no.. = T))==0){
		unlink(chunkFolder,recursive = T,force = T)
	} else {
		warning(sprintf('%s is not empty after merging, leaving it in place',chunkFolder),
						call.=FALSE,immediate.=TRUE)
	}
	if(verbosity>0){cat('done\n')}
	invisible(NULL)
}

# mergePerVarFilesLegacy ####
# The merge for the chunk layout where every output type has its own tree of
# chunk files, rather than the single csv tree in PerVarChunks/.
mergePerVarFilesLegacy <- function(outputFolder,baseWD,verbosity=1,parStrat=2,
																	 compressCsv=T,rdsCompress=TRUE,fullPrecision=TRUE,
																	 outputTypeFoldersOverride = NULL,
																	 varNamesOverride = NULL){
	if(!is.null(outputTypeFoldersOverride)){
		outputTypeFolders <- outputTypeFoldersOverride
	} else {
		outputTypeFolders <- basename(list.dirs(outputFolder,recursive = F))
	}
	for(outputTypeFolder in outputTypeFolders){
		outputType <- strsplit(outputTypeFolder,'-')[[1]][2]
		if(is.na(outputType) || !outputType %in% perVarOutputTypes){
			next
		}
		outputTypeFolder <- file.path(outputFolder,outputTypeFolder)
		if(verbosity>0){cat(paste('Entering',outputTypeFolder,'\n'))}
		if(!is.null(varNamesOverride)){
			varNames <- varNamesOverride
		} else {
			varNames <- basename(list.dirs(outputTypeFolder,recursive = F))
		}
		if(verbosity>0){cat(sprintf('Found %i variable sub folder(s)\n',length(varNames)))}
		if(length(varNames)==0){
			next
		}
		if(parStrat==1){
			for(v.i in 1:length(varNames)){
				if(verbosity>0){cat(sprintf('(%i of %i) ',v.i,length(varNames)))}
				workerMergePerVarFilesLegacy(v.i,varNames=varNames,
																		 outputTypeFolder=outputTypeFolder,
																		 outputType=outputType,verbosity=verbosity,
																		 compressCsv=compressCsv,rdsCompress=rdsCompress,
																		 fullPrecision=fullPrecision)
			}
		} else {
			if(verbosity>0){cat(sprintf('Parallel proccessing all vars with %i workers\n',numWorkersFileMerge))}
			clFileMerge <- startFileMergeCluster(numWorkersFileMerge,baseWD)
			gobble <- parLapplyLB(clFileMerge,1:length(varNames),workerMergePerVarFilesLegacy,
														varNames=varNames,
														outputTypeFolder=outputTypeFolder,
														outputType=outputType,
														compressCsv=compressCsv,
														rdsCompress=rdsCompress,
														fullPrecision=fullPrecision,
														chunk.size = 1)
			stopCluster(clFileMerge)
		}
		if(verbosity>0){cat('done\n')}
	}
	invisible(NULL)
}

# readPerVarFile ####
readPerVarFile <- function(file,outputType=NULL){
	if(is.null(outputType)){
		outputType <- tools::file_ext(file)
		if(outputType==''){
			stop('No outputType and no file ext\n')
		}
		if(outputType=='gz'){
			outputType <- 'csv'
		}
	}
	# remove file extensions
	fileNoExt <- gsub('\\.gz','',file)
	fileNoExt <- gsub('\\.csv','',fileNoExt)
	fileNoExt <- gsub('\\.RDS','',fileNoExt)
	if(outputType=='csv'){
		csvFile <- paste0(fileNoExt,'.csv')
		if(!file.exists(csvFile) && file.exists(paste0(fileNoExt,'.csv.gz'))){
			csvFile <- paste0(fileNoExt,'.csv.gz')
		}
		varData <- data.table::fread(csvFile,nThread=1,showProgress=FALSE,
																 integer64='double')
		coercePerVarTypes(varData)
		# the callers index the result with data.frame semantics (varData[,-1]),
		# which means something else on a data.table
		data.table::setDF(varData)
		# some files carry the X that read.csv puts in front of the year column names
		colnames(varData) <- gsub('^X(?=\\d+$)','',colnames(varData),perl = T)
		return(varData)
	} else if(outputType=='RDS'){
		return(readRDS(paste0(fileNoExt,'.RDS')))
	}
}

# writePerVarFile ####
writePerVarFile <- function(varData,file,outputType=NULL,compressCsv=T,rdsCompress=TRUE,
													fullPrecision=TRUE){
	if(is.null(outputType)){
		outputType <- tools::file_ext(file)
		if(outputType==''){
			stop('No outputType and no file ext\n')
		}
		if(outputType=='gz'){
			outputType <- 'csv'
			compressCsv <- T
		}
	}
	# remove file extensions
	fileNoExt <- gsub('\\.gz','',file)
	fileNoExt <- gsub('\\.csv','',fileNoExt)
	fileNoExt <- gsub('\\.RDS','',fileNoExt)
	if(outputType=='csv'){
		if(compressCsv){
			fwritePerVarCsv(varData,paste0(fileNoExt,'.csv.gz'),
										fullPrecision=fullPrecision,compress='gzip')
			# fwrite compresses straight into the .gz, so an uncompressed leftover of
			# an earlier run would otherwise shadow it in readPerVarFile
			unlink(paste0(fileNoExt,'.csv'),force=TRUE)
		} else {
			fwritePerVarCsv(varData,paste0(fileNoExt,'.csv'),
										fullPrecision=fullPrecision,compress='none')
		}
	} else if(outputType=='RDS'){
		saveRDS(varData,paste0(fileNoExt,'.RDS'),compress=rdsCompress)
	} else {
		stop('No outputType\n')
	}
}



















