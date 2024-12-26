require(parallel)

if(!exists('location.output')){
	stop('run this only after config.R\n')
}

# config ####
numWorkers <- detectCores()
clusterType <- 'psock'

# tmpfs location for the worker directories to not churn the hard drive
# and be faster
# alternative path to /dev/shm would be /run/user/####/ where #### is the uid
tmpfsDir <- paste0('/run/user/',system('id -u',intern = T),'/rwork/')


## cluster cleanup ####
cat('cluster cleanup...')
# stop cluster
tryCatch(stopCluster(cl),error=function(e){})
# clean up working directories
unlink('workerDirs',recursive=T,force=T)
unlink(tmpfsDir,recursive=T,force=T)
cat('done\n')



# cluster setup ####
cat('cluster setup...')

dir.create(tmpfsDir,recursive = F)
system(paste('ln -s',tmpfsDir,'workerDirs'))

baseWD <- getwd()
workDirBasename <- 'workDir_'
# start cluster
if(clusterType=='fork'){
	cl <- makeForkCluster(numWorkers)
} else if (clusterType=='psock'){
	cl <- makePSOCKcluster(numWorkers)
}
gobble <- clusterEvalQ(cl,tools::psnice(value=15))
gobble <- clusterExport(cl,list('location.output','baseWD',
											'chunkSizePerWorker','workDirBasename'))
gobble <- clusterEvalQ(cl,source(file.path(baseWD,'initialise.R')))
# make working directories
workers <- 1:length(cl)
gobble <- clusterApply(cl,workers,function(i){
	workerID <- i
	dir.create(file.path('workerDirs',paste0(workDirBasename,i)),showWarnings = F,recursive = T)
	setwd(file.path('workerDirs',paste0(workDirBasename,i)))
})
# clusterEvalQ(cl,getwd())
# copy over the model and simulator
gobble <- clusterApply(cl,workers,function(i){
	file.copy(file.path(baseWD,location.frida),getwd(),recursive=T)
	file.copy(file.path(baseWD,location.stella),getwd(),recursive=T)
	file.copy(file.path(baseWD,'frida_info.csv'),getwd())
})

# copy extra vars if present (for restarting the cluster during testing)
if(exists('sampleParms')){clusterExport(cl,list('sampleParms'))}
if(exists('calDat')){clusterExport(cl,list('calDat'))}
if(exists('resSigma')){clusterExport(cl,list('resSigma'))}




cat('done\n')
