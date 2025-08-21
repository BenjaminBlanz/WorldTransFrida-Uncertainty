# settings for runPolicyAnalysis.R

source('funParmSpace.R')

expID <- 'testing2'
addAutoNameToExpID <- T

# number of joint policy scenarios
numInitialJointPol <- 1e5

# probability of status quo policy per policy domain
# in sampling joint policies
# This ensures that not every joint policy run includes 
# variation in all policies.
# Higher values cause more sparse joint policies
nullPolProb <- 0.5

# parallel ####
useSLURM <- TRUE
maxJobsQueue <- 10
numWorkers <- parallel::detectCores()
numWorkersFileMerge <- floor(numWorkers/3)
# How large the chunks of work are, smaller means more frequent pauses to write out
# itermediate results (and update the diagnostic output).
# smaller also means increaseing the total number of files, which can cause file system 
# slow down.
chunkSizePerWorker <- 50
# tyoe of cluster. PSOCK allows connections across a network
# FORK forks the currently running process, but with copy on write memory
# sharing
clusterType <- 'psock'
# name of the directory containing the policy worker dirs
origName.workDir <- name.workDir <- 'policy-WorkDirs'
# basename of the worker dirs within the above dir.
name.workerDirBasename <- 'policy-WorkerDir_'
# tmpfs location for the worker directories to not churn the hard drive
# and be faster
# typical options on linux are /dev/shm or /run/user/####/ where #### is the uid
# if both of these are unavailable use notTMPFS or some other arbitrary location on disk
# tmpfsBaseDir <- paste0('/run/user/',system('id -u',intern = T),'/rwork')
tmpfsBaseDir <- paste0('/dev/shm/',system('id -u',intern = T),'/rwork')
# tmpfsBaseDir <- 'notTMPFS'

# output ####
perVarOutputTypes <- c('csv','RDS')
# gzip the merged csv files
compressCsv <- TRUE

# locations and names ####
# location of frida/stella for running
baselocation.frida <-location.frida <- './FRIDAforPolicyAnalysis'
# location for setting parameters for FRIDA
# e.g. turnig climate feedbacks on or off
# or policy
location.frida.configs <- './FRIDA-configs'
baselocation.stella <- location.stella <- './Stella_Simulator_Linux'
# location frida/stella is stored while the above is located in tmpfs
location.frida.storage <- './FRIDAforPolicyAnalysis-store'
location.stella.storage <- './Stella_Simulator_Linux-store'
locaion.baseOutput <- file.path('policy-workOutput')
# names of files written to the FRIDA Data directory for the running and export
name.fridaExportVarsFile <- 'varsForExport.txt'
name.fridaInputFile <- 'uncertainty_analysis_paramter_values.csv'
name.fridaOutputFile <- 'uncertainty_analysis_exported_variables.csv'
location.frida.info <- './FRIDA-info'
name.frida_info <- 'frida_info.csv'
name.frida_extra_variables_to_export_list <- 'frida_extra_variables_to_export_list.csv'


location.singleDomainPolicyFiles <- file.path('policy-singleDomainPolicyMatrices')

policyFiles <- list.files(location.singleDomainPolicyFiles)

if(addAutoNameToExpID){
	name.output <- gsub('\\.','_',paste0(expID,'-N-',numInitialJointPol,'-nPr-',nullPolProb,'-polFiles-',
																		 paste(cleanNames(tools::file_path_sans_ext(policyFiles)),collapse='-')))
} else {
	name.output <- expID
}
name.workDir <- paste0(name.workDir,'-',name.output)
expID <- name.output
location.output <- file.path(locaion.baseOutput,name.output)
origTmpfsDir <- tmpfsDir <- file.path(tmpfsBaseDir,name.output)
dir.create(location.output,recursive = T,showWarnings = F)

# hardcoding this here, as the below code relies on it.
# this should become standard also in uncertainty in the future
# Do not change!
writePerWorkerFiles <- TRUE
doNotReturnRunDataSavePerWorkerOnly <- TRUE
