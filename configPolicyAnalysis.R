# settings for runPolicyAnalysis.R

source('funParmSpace.R')

# Set the expID to a something descriptive of the run.
# The expID is used to name the output folder, among other things.
# addAutoNameToExpID can be used to automatically generate something descriptive based 
# on the run configuration but this will cause too long filenames and the run will crash.
expID <- 'testing2'
# This adds a description of the configured policies to sample to the expID
# Set this to false if you intend to sample a wide range of policies as
# the filenames will be too long to handle and the run will crash.
addAutoNameToExpID <- F

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
# How long will it take approximately? Job will be killed after this time!
# But max 8:00 hours and the shorter the run is, the earlier the job gets run
SLURMhours <- '8'
SLURMminutes <- '00'
# Does the job need larger memory?
# (the larger the memory you request, the longer the job might sit in the queue, if the machine is full)
memorySize='256G' # can be ['256G' '512G', '1024G'], '256G' sometimes fails with 100,000 samples
SLURMmemorySize <- '256G'
# Use a different group account for ressources? Which partition?
SLURMaccount <- 'uc1275'
SLURMpartition <- 'compute'
# Enter Email here in case you want to receive a mail, when the job failed
SLURMemail <- ''
# bumber of workers for each of the slurm jobs
numWorkers <- 256
numWorkersFileMerge <- min(numWorkers,ceiling(1e6/numInitialJointPol))
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
# valid types are 'RDS' and 'csv'
perVarOutputTypes <- c('RDS') #'csv'
# gzip the csv files
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

# Policies to sample ####
# By default all policy files present in location.singleDomainPolicyFiles will be sampled
# To sample only a subset either specify a location containing only the desired sample above
# or provide a vector of strings with the file names of those you want to sample here.
policyFiles <- list.files(location.singleDomainPolicyFiles)


# plotting ####

filterSpec <- list()
filterSpec$inflation_inflation_rate <- c('gtabs',0.2)
filterSpec$gdpgr <- c('gtabs',0.2)
filterSpec$stagr <- c('gtabs',0.2)
filterSpec$gdp_real_gdp_in_2021c <- c('ltval',100)

numPlotThreads <- 10
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
plotType <- 2
plot.numColLevels <- 26
plot.palletteName <- "grDevices::RdPu"
plot.palletteNameSOW <- "grDevices::Greens"
plotWidth <- 20
plotHeight <- 20
plotUnit <- 'cm'
plotRes <- 150
# plot lines are
# 1. Area of the SOW outcomes
# 2. Area of the median outcomes
plot.lty <- c('dotted','solid')#,'dotdash','dotted')
plot.lwd <- c(1,3)
plot.lcol <- c(1,1) # colors of the outlines
plot.col <- c(gray(0.7,0.5),gray(0.2,0.5)) # colors of the area fills
plot.relyrange <- c(0.0001,0.9999)

# automatic stuff ####
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
