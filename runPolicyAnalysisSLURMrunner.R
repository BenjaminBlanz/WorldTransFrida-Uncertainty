# expID and unitID will be overwritten by the main runner
# load the config
source('initialise.R')
source('configPolicyAnalysis.R')


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
	} else if(status=='prepared'){
		# TODO: submit slurm job for this workUnit	
		if(useSLURM){
			if(numJobsQueued < maxJobsQueue){
				numJobsQueued <- numJobsQueued + 1
				system(paste0('./submit_PolicyAnalysisSLURMJob.sh',
					' -o ',location.output,
					' -u ',unitID,
					' -w ',numWorkers),
					intern= F, wait = F)
				status <- 'submitted (new)'
			}
		} else {
			system(paste0('Rscript runPolicyAnalysisWorkUnit.R ',unitID))
		}
	}
	# other status e.g. failed is simply reported in the table
	statuses[workUnitDir] <- status
}
cat("summary of work units' statuses\n")
print(table(statuses))
