# settings for runPolicyAnalysis.R

source('funParmSpace.R')

# Set the expID to a something descriptive of the run.
# The expID is used to name the output folder, among other things.
# addAutoNameToExpID can be used to automatically generate something descriptive based 
# on the run configuration but this will cause too long filenames and the run will crash.
expID <- 'AllPolicies1e6-moreExports'
# This adds a description of the configured policies to sample to the expID
# Set this to false if you intend to sample a wide range of policies as
# the filenames will be too long to handle and the run will crash.
addAutoNameToExpID <- F

# number of joint policy scenarios
numInitialJointPol <- 1e6

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

medianSOWid <- 6

# TODO: explain the filterSpec
# these are the filters that drop policies that are outside the functional scope
# of the model (while the math may work in some of these cases the conceptual model 
# and hence the interpretability does not)
# Types:
#  gtabs: absolute value greater than the filter level
#  ltabs: absolute value less than the filter level
#  gtval: value greater than the filter level
#  ltval: value less than the filter level
filterSpec <- list()
filterSpec$inflation_inflation_rate <- list()
filterSpec$inflation_inflation_rate$type <- 'gtabs'
filterSpec$inflation_inflation_rate$level <- 0.5
filterSpec$inflation_inflation_rate$allowedTransgressions <- 2
filterSpec$real_gdp_growth_rate <- list()
filterSpec$real_gdp_growth_rate$type <- 'ltval'
filterSpec$real_gdp_growth_rate$level <- -0.1
filterSpec$real_gdp_growth_rate$sowID <- medianSOWid # the median SOW
filterSpec$surface_temperature_anomaly_growth_rate <- list()
filterSpec$surface_temperature_anomaly_growth_rate$type <- 'gtabs'
filterSpec$surface_temperature_anomaly_growth_rate$level <- 0.5
filterSpec$surface_temperature_anomaly_growth_rate$allowedTransgressions <- 0
filterSpec$gdp_real_gdp_in_2021c <- list()
filterSpec$gdp_real_gdp_in_2021c$type <- 'ltval'
filterSpec$gdp_real_gdp_in_2021c$level <- 5e4 # roughly the 1980 level of GDP
filterSpec$gdp_real_gdp_in_2021c$allowedTransgressions <- 0
filterSpec$gdp_real_gdp_in_2021c$years <- 2023:2150 # the years in which this filter applies

#these are the filters that drop policies with "undersirable" outcomes
desiredFilterSpec <- list()
desiredFilterSpec$gdp_real_gdp_in_2021c <- list()
desiredFilterSpec$gdp_real_gdp_in_2021c$type <- 'ltval'
desiredFilterSpec$gdp_real_gdp_in_2021c$level <- 1.5e5 # a bit less than the 2023 median level
desiredFilterSpec$gdp_real_gdp_in_2021c$allowedTransgressions <- 4 # ~37% of cases with default 11 runs
desiredFilterSpec$gdp_real_gdp_in_2021c$years <- 2023:2150 # the years in which this filter applies
desiredFilterSpec$energy_balance_model_surface_temperature_anomaly <- list()
desiredFilterSpec$energy_balance_model_surface_temperature_anomaly$type <- 'gtval'
desiredFilterSpec$energy_balance_model_surface_temperature_anomaly$level <- 2
desiredFilterSpec$energy_balance_model_surface_temperature_anomaly$allowedTransgressions <- 5 # 5/11 allowed to transgress
desiredFilterSpec$energy_balance_model_surface_temperature_anomaly$years <- 2150 # allow overshoot
desiredFilterSpec$gdp_future_year_in_recession <- list()
desiredFilterSpec$gdp_future_year_in_recession$type <- 'gtval'
desiredFilterSpec$gdp_future_year_in_recession$level <- 10
desiredFilterSpec$gdp_future_year_in_recession$allowedTransgressions <- 4 # ~37% of cases with default 11 runs
# after which consecutively applied desired filter should we run plots?
# default is to only plot with all the filters applied, to save time
filtersToPlot <- c(length(desiredFilterSpec))

#select a specific run to highlight and output the policies
selectedRunSpec <- list()
selectedRunSpec$var <- 'gdp_real_gdp_in_2021c'
selectedRunSpec$year <- 2150
selectedRunSpec$optimize <- 'max'
selectedRunSpec$sow <- 6

selectedRun.lty <- 1
selectedRun.lwd <- 2
selectedRun.col <- 'black'
selectedRunEnsemble.lty <- 2
selectedRunEnsemble.lwd <- 1
selectedRunEnsemble.col <- 'black'

numPlotThreads <- 15
location.plots <- 'figures'
# plot 
plotTypes <- c(2,3)
plot.numColLevels <- 26
plot.palletteName <- "grDevices::RdPu"
plot.palletteNameSOW <- "grDevices::Greens"
plot.palletteOmmitEntries <- 1:5
plot.palletteReplaceEntries.idc <- c(1) # the first entry is for values >=0 and <1
plot.palletteReplaceEntries.cols <- c('white')
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
# pre define plot ylimits for selected variables (uses the cleanName(varName) as the key)
plot.ylimOverrides <- list()
plot.ylimOverrides$gdp_real_gdp_in_2021c <- c(0,2.5e6)

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
