# expID and unitID will be overwritten by the main runner
# load the config
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
if(!file.exists(file.path(baseWD,location.output,paste0('workunit-',i),'pdpMeta.RDS'))||
	 !file.exists(file.path(baseWD,location.output,paste0('workunit-',i),'pdp.lst.RDS'))||
	 !file.exists(file.path(baseWD,location.output,paste0('workunit-',i),'jointPolicies.RDS'))){
	stop('Missing files, runPolicyAnalysisSLURM has not completed\n')	
}

# check list of work units
workUnitDirs <- list.files(file.path(location.output),
												pattern = 'workunit-')
statuses <- c()
for(workUnitDir in workUnitDirs){
	if(!file.exists(file.path(location.output,paste0('workunit-',i),'samplePoints.RDS'))){
		stop('Missing files for', workUnitDir,'\n')
	}
	# check for already running process
	#TODO read status file
	if(status=='started' || status=='running'){
		# skip this unit it is being processed 
		break
	}
	if(status=='failed'){
		write('cleanup',file.path(location.output,paste0('workunit-',i),'log.txt',append=T))
		cleanup(workUnitDir)
		status=='prepared'
	}
	if(status=='prepared'){
		# TODO: submit slurm job for this workunit	
	}
	statuses[workUnitDir] <- status
}
