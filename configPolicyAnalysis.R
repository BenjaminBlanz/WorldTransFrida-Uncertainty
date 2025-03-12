# settings for runPolicyAnalysis.R

source('funParmSpace.R')

# number of joint policy scenarios
numInitialJointPol <- 1e5

# probability of status quo policy per policy domain
# in sampling joint policies
# This ensures that not every joint policy run includes 
# variation in all policies.
# Higher values cause more sparse joint policies
nullPolProb <- 0.5

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
name.workerDirBasename <- 'policy-WorkerDir_'
name.workDir <- 'policy-WorkDirs'
# tmpfs location for the worker directories to not churn the hard drive
# and be faster
# typical options on linux are /dev/shm or /run/user/####/ where #### is the uid
# if both of these are unavailable use notTMPFS or some other arbitrary location on disk
tmpfsBaseDir <- paste0('/run/user/',system('id -u',intern = T))
# tmpfsBaseDir <- paste0('/dev/shm/',system('id -u',intern = T))
# tmpfsBaseDir <- 'notTMPFS'

# locations and names ####
# location of frida/stella for running
location.frida <- './FRIDAforPolicyAnalysis'
location.stella<- './Stella_Simulator_Linux'
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


location.singleDomainPoliciyFiles <- file.path('policy-singleDomainPolicyMatrices')


name.output <- gsub('\\.','_',paste0('N-',numInitialJointPol,'-nPr-',nullPolProb))
location.output <- file.path(locaion.baseOutput,name.output)
tmpfsDir <- file.path(tmpfsBaseDir,'rwork',name.output)
dir.create(location.output,recursive = T,showWarnings = F)
