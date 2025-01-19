# clean up first
if(!exists('cl')<0&&try(clusterEvalQ(cl,1+1)[[1]],silent=T)!=2){
	if(file.exists('workerDirs')){
		source('cleanup.R')
	}
	# create the tmpfsDir and link to workerDirs
	if(!file.exists(tmpfsDir)){
		dir.create(tmpfsDir,recursive = T,showWarnings = F)
	}
	if(!file.exists('workerDirs')){
		system(paste('ln -s',tmpfsDir,'workerDirs'))
	}
	
	# set up the tmpfs for the single threaded runs
	if(!file.exists(location.frida.storage)){
		system(paste('mv',location.frida,location.frida.storage))
		system(paste('cp -r',location.frida.storage,file.path(tmpfsDir,location.frida)))
		system(paste('ln -s',file.path(tmpfsDir,location.frida),location.frida))
	}
	if(!file.exists(location.stella.storage)){
		system(paste('mv',location.stella,location.stella.storage))
		system(paste('cp -r',location.stella.storage,file.path(tmpfsDir,location.stella)))
		system(paste('ln -s',file.path(tmpfsDir,location.stella),location.stella))
	}
}
