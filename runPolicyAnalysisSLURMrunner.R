# expID and unitID will be overwritten by the main runner
# load the config
source('initialise.R')
source('configPolicyAnalysis.R')

if(!exists('restartFailed')){
	restartFailed <- F
}

# Check the setup
# setup writes out sample points for the workers
# slurm runner, this file, checks the chunk folders
# check consists of
# 	- verifying the existence of sample points
# 	- checking if slurm run is already active and/or results already exist
# 		- sending of a run for the chunk if it is not active
# 		- skipping if the run is active or results exist

# check that files exist
if(!file.exists(file.path(location.output,'pdpMeta.RDS')) ||
	!file.exists(file.path(location.output,'pdp.lst.RDS')) ||
	!file.exists(file.path(location.output,'jointPolicies.RDS'))){
	stop('Missing policy files runPolicyAnalysisSLURM has not completed\n')
}

# check list of work units
workUnitDirs <- list.files(file.path(location.output,'workUnits'),
												pattern = 'workUnit-')
statuses <- c()
numJobsQueued <- 0
for(workUnitDir in workUnitDirs){
	unitID <- as.numeric(gsub('workUnit-','',workUnitDir))
	if(!file.exists(file.path(location.output,'workUnits',workUnitDir,'samplePoints.RDS')) ||
		 !file.exists(file.path(location.output,'workUnits',workUnitDir,'status.txt'))){
		stop('Missing files for ', workUnitDir,'\nrunPolicyAnalysisSLURM has not completed\n')
	}
	# check for already running process
	status <- readLines(file.path(location.output,'workUnits',workUnitDir,'status.txt'))
	if(status=='submitted'||status=='started'||status=='running'){
		numJobsQueued <- numJobsQueued+1
	}	else {
		if(status=='failed' & restartFailed){
			status=='prepared'
			write('prepared',file.path(location.output,'workUnits',paste0('workUnit-',i),'status.txt'))
		}
		if(status=='prepared'){
			# TODO: submit slurm job for this workUnit	
			if(useSLURM){
				if(numJobsQueued < maxJobsQueue){
					numJobsQueued <- numJobsQueued + 1
					system(paste0('./submit_PolicyAnalysisSLURMJob.sh',
												' -h ',SLURMhours,
												' -m ',SLURMminutes,
												' -M ',SLURMmemorySize,
												' -a ',SLURMaccount,
												' -p ',SLURMpartition,
												' -e ',SLURMemail,
												' -o ',location.output,
												' -s ',name.output,
												' -u ',unitID,
												' -w ',numWorkers),
						intern = T, wait = F,ignore.stdout = T)
					status <- 'submitted (new)'
				}
			} else {
				system(paste0('Rscript runPolicyAnalysisWorkUnit.R ',unitID))
			}
		}
	}
	# other status e.g. failed is simply reported in the table
	statuses[workUnitDir] <- status
}
cat("summary of work units' statuses\n")
print(table(statuses))

# run the merger
if('completed' %in% names(table(statuses))){
	numComplete <- table(statuses)['completed']
} else {
	numComplete <- 0
}
if(numComplete < length(statuses)){
	mergeStatus <- 'waiting for runs to complete'
} else if(file.exists(file.path(location.output,'mergeStatus.txt'))){
	mergeStatus <- readLines(file.path(location.output,'mergeStatus.txt'))
} else {
	mergeStatus <- 'ready for merge'
} 
write(mergeStatus,file.path(location.output,'mergeStatus.txt'))
if(mergeStatus=='ready for merge'){
	cat('submitting file merge \n')
	system(paste0('./submit_PolicyAnalysisMergeSLURMjob.sh',
								' -h ',SLURMhours,
								' -m ',SLURMminutes,
								' -M ',SLURMmemorySize,
								' -a ',SLURMaccount,
								' -p ',SLURMpartition,
								' -e ',SLURMemail,
								' -o ',location.output,
								' -s ',name.output,
								' -w ',numWorkersFileMerge),
				 intern = T, wait = F,ignore.stdout = T)
} else if(mergeStatus=='submitted'){
	cat('Merge already submitted and waiting for allocation, check squeue\n')
} else if(mergeStatus=='running'){
	cat('Merge already running, check squeue\n')
} else if(mergeStatus=='sbatch failed'){
	cat('Failed to submit merge job\n')
} else if(mergeStatus=='failed'){
	cat('Merge failed debug by running runPolicyAnalysisMerger.R interactively\n')
} else if(mergeStatus=='completed'){
	cat('Merge completed\n')
}

# run the plots

