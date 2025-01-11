
library(parallel,quietly=T,warn.conflicts = F)


baseWD <- getwd()
workDirBasename <- 'workDir_'

if(!exists('location.output')){
	stop('run this only after config.R\n')
}

# do not stop and start the cluster, instead just reinitialise
Sys.sleep(1) # prevent race condition if the cluster was busy or sth.
if(exists('cl')){
	result <- try(clusterEvalQ(cl,1+1)[[1]],silent=T)
}
if(exists('cl')&&is.numeric(result)&&result==2){
	gobble <- clusterApply(cl,workers,function(i){
		workerID <<- i
	})
	gobble <- clusterExport(cl,list('baseWD','workDirBasename'))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'config.R')))
	if(exists('sampleParms')){clusterExport(cl,list('sampleParms'))}
	if(exists('calDat')){clusterExport(cl,list('calDat'))}
	if(exists('resSigma')){clusterExport(cl,list('resSigma'))}
	gobble <- clusterEvalQ(cl,tools::psnice(value=15))
	gobble <- clusterExport(cl,list('location.output','baseWD',
																	'chunkSizePerWorker','workDirBasename'))
	gobble <- clusterEvalQ(cl,suppressPackageStartupMessages({
		library(Rmpfr,quietly=T,warn.conflicts = F) # use to calculate the likelihood from loglikelihood
		library(optimx,quietly=T,warn.conflicts = F)
		library(tictoc)
		library(SobolSequence,quietly = T)
	}))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'funRunFRIDA.R')))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'funPlot.R')))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'funParmSpace.R')))
	# copy over the model and simulator
	gobble <- clusterApply(cl,workers,function(i){
		Sys.sleep(workerID*0.1) # stagger the file copying to not DOS the server
		file.copy(file.path(baseWD,location.frida),getwd(),recursive=T,overwrite = T)
		file.copy(file.path(baseWD,location.stella),getwd(),recursive=T,overwrite = T)
		file.copy(file.path(baseWD,location.frida.info,name.frida_info),getwd(),overwrite = T)
	})
} else {
	# cluster setup ####
	cat('cluster setup...')
	
	if(!file.exists('workerDirs')){
		source('setupTMPFS.R')
	}
	
	# start cluster
	if(clusterType=='fork'){
		cl <- makeForkCluster(numWorkers)
	} else if (clusterType=='psock'){
		cl <- makePSOCKcluster(numWorkers,setup_strategy='sequential') # setup sequential to not DOS the server
	}
	gobble <- clusterExport(cl,list('baseWD','workDirBasename'))
	gobble <- clusterEvalQ(cl,tools::psnice(value=15))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'initialise.R')))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'config.R')))
	# make working directories
	workers <- 1:length(cl)
	gobble <- clusterApply(cl,workers,function(i){
		workerID <<- i
		Sys.sleep(workerID*0.1) # stagger the file copying to not DOS the server
		system(paste('rm -r',file.path('workerDirs',paste0(workDirBasename,i))),ignore.stdout = T,ignore.stderr = T)
		dir.create(file.path('workerDirs',paste0(workDirBasename,i)),showWarnings = F,recursive = T)
		setwd(file.path('workerDirs',paste0(workDirBasename,i)))
	})
	# clusterEvalQ(cl,getwd())
	# copy over the model and simulator
	gobble <- clusterApply(cl,workers,function(i){
		Sys.sleep(workerID*0.1) # stagger the file copying to not DOS the server
		file.copy(file.path(baseWD,location.frida),getwd(),recursive=T,overwrite = T)
		file.copy(file.path(baseWD,location.stella),getwd(),recursive=T,overwrite = T)
		file.copy(file.path(baseWD,location.frida.info,name.frida_info),getwd(),overwrite = T)
	})
	
	# copy extra vars if present (for restarting the cluster during testing)
	if(exists('sampleParms')){clusterExport(cl,list('sampleParms'))}
	if(exists('calDat')){clusterExport(cl,list('calDat'))}
	if(exists('resSigma')){clusterExport(cl,list('resSigma'))}
	cat('done\n')
}
