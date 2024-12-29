
library(parallel,quietly=T,warn.conflicts = F)

if(!exists('location.output')){
	stop('run this only after config.R\n')
}

# do not stop and start the cluster, instead just reinitialise
if(exists('cl')&&try(clusterEvalQ(cl,1+1)[[1]],silent=T)==2){
	gobble <- clusterApply(cl,workers,function(i){
		workerID <<- i
	})
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
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'config.R')))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'funRunFRIDA.R')))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'funPlot.R')))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'funParmSpace.R')))
} else {
	# cluster setup ####
	cat('cluster setup...')
	
	if(!file.exists('workerDirs')){
		source('setupTMPFS.R')
	}
	
	baseWD <- getwd()
	workDirBasename <- 'workDir_'
	# start cluster
	if(clusterType=='fork'){
		cl <- makeForkCluster(numWorkers)
	} else if (clusterType=='psock'){
		cl <- makePSOCKcluster(numWorkers,setup_strategy='sequential') # setup sequential to not DOS the server
	}
	gobble <- clusterEvalQ(cl,tools::psnice(value=15))
	gobble <- clusterExport(cl,list('location.output','baseWD',
												'chunkSizePerWorker','workDirBasename'))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'initialise.R')))
	# make working directories
	workers <- 1:length(cl)
	gobble <- clusterApply(cl,workers,function(i){
		workerID <<- i
		Sys.sleep(workerID*0.5) # stagger the file copying to not DOS the server
		system(paste('rm -r',file.path('workerDirs',paste0(workDirBasename,i))),ignore.stdout = T,ignore.stderr = T)
		dir.create(file.path('workerDirs',paste0(workDirBasename,i)),showWarnings = F,recursive = T)
		setwd(file.path('workerDirs',paste0(workDirBasename,i)))
	})
	# clusterEvalQ(cl,getwd())
	# copy over the model and simulator
	gobble <- clusterApply(cl,workers,function(i){
		Sys.sleep(workerID*0.5) # stagger the file copying to not DOS the server
		file.copy(file.path(baseWD,location.frida),getwd(),recursive=T)
		file.copy(file.path(baseWD,location.stella),getwd(),recursive=T)
		file.copy(file.path(baseWD,'frida_info.csv'),getwd())
	})
	
	# copy extra vars if present (for restarting the cluster during testing)
	if(exists('sampleParms')){clusterExport(cl,list('sampleParms'))}
	if(exists('calDat')){clusterExport(cl,list('calDat'))}
	if(exists('resSigma')){clusterExport(cl,list('resSigma'))}
	cat('done\n')
}
