
library(parallel,quietly=T,warn.conflicts = F)


baseWD <- getwd()
if(!exists('name.workerDirBasename')){
	name.workerDirBasename <- 'workDir_'
}
if(!exists('name.workDir')){
	name.workDir <- 'workerDirs'
}

if(!exists('location.output')){
	stop('run this only after config.R\n')
}

# do not stop and start the cluster, instead just reinitialise
if(exists('cl')){
	result <- try(clusterEvalQ(cl,1+1)[[1]],silent=T)
}
if(exists('cl')&&is.numeric(result)&&result==2){
	gobble <- clusterApply(cl,workers,function(i){
		workerID <<- i
	})
	if(exists('sampleParms')){clusterExport(cl,list('sampleParms'))}
	if(exists('calDat')){clusterExport(cl,list('calDat'))}
	if(exists('resSigma')){clusterExport(cl,list('resSigma'))}
	gobble <- clusterEvalQ(cl,tools::psnice(value=15))
	gobble <- clusterExport(cl,list('baseWD',
																	'location.output',
																	'name.workDir','name.workerDirBasename',
																	'location.frida','location.stella',
																	'location.frida.configs',
																	'name.frida_info',
																	'name.fridaExportVarsFile',
																	'name.fridaInputFile',
																	'name.fridaOutputFile',
																	'name.frida_extra_variables_to_export_list',
																	'chunkSizePerWorker'))
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
		# Sys.sleep(workerID*0.1) # stagger the file copying to not DOS the server
		file.copy(file.path(baseWD,location.frida),getwd(),recursive=T,overwrite = T)
		file.copy(file.path(baseWD,location.stella),getwd(),recursive=T,overwrite = T)		
		if(exists('location.frida.info')){
			file.copy(file.path(baseWD,location.frida.info,name.frida_info),getwd(),overwrite = T)
		}
	})
} else {
	# cluster setup ####
	cat('cluster setup...')
	
	if(!file.exists(name.workDir)){
		source('setupTMPFS.R')
	}
	
	# start cluster
	cat('start cluster...')
	if(clusterType=='fork'){
		cl <- makeForkCluster(numWorkers)
	} else if (clusterType=='psock'){
		cl <- makePSOCKcluster(numWorkers)#,setup_strategy='sequential') # setup sequential to not DOS the server
		# add ,outfile='' to enable cluster output on the console
	}
	cat('exports...')
	gobble <- clusterExport(cl,list('baseWD',
																	'location.output',
																	'name.workDir','name.workerDirBasename',
																	'location.frida','location.stella',
																	'name.frida_info',
																	'name.fridaExportVarsFile',
																	'name.fridaInputFile',
																	'name.fridaOutputFile',
																	'name.frida_extra_variables_to_export_list',
																	'chunkSizePerWorker'))
	if(exists('location.frida.info')){
		gobble <- clusterExport(cl,'location.frida.info')	
	}
	gobble <- clusterEvalQ(cl,tools::psnice(value=15))
	gobble <- clusterEvalQ(cl,source(file.path(baseWD,'initialise.R')))
	cat('work dirs...')
	# make working directories
	workers <- 1:length(cl)
	gobble <- clusterApply(cl,workers,function(i){
		workerID <<- i
		# Sys.sleep(workerID*0.1) # stagger the file copying to not DOS the server
		system(paste('rm -r',file.path(name.workDir,paste0(name.workerDirBasename,i))),ignore.stdout = T,ignore.stderr = T)
		dir.create(file.path(name.workDir,paste0(name.workerDirBasename,i)),showWarnings = F,recursive = T)
		setwd(file.path(name.workDir,paste0(name.workerDirBasename,i)))
	})
	# clusterEvalQ(cl,getwd())
	# copy over the model and simulator
	gobble <- clusterApply(cl,workers,function(i){
		# Sys.sleep(workerID*0.1) # stagger the file copying to not DOS the server
		file.copy(file.path(baseWD,location.frida),getwd(),recursive=T,overwrite = T)
		file.copy(file.path(baseWD,location.stella),getwd(),recursive=T,overwrite = T)
		if(exists('location.frida.info')){
			file.copy(file.path(baseWD,location.frida.info,name.frida_info),getwd(),overwrite = T)
		}
	})
	if(!exists('skipExtraVars')){
		cat('extra vars...')
		# copy extra vars if present (for restarting the cluster during testing)
		if(exists('sampleParms')){clusterExport(cl,list('sampleParms'))}
		if(exists('calDat')){clusterExport(cl,list('calDat'))}
		if(exists('resSigma')){clusterExport(cl,list('resSigma'))}
	}
	cat('done\n')
}
