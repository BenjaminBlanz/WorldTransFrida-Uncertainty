
# if there is a cluster running, we must already have setup the TMPFS,
# so do nothing further.
if(!exists('cl') &&
	 try(clusterEvalQ(cl,1+1)[[1]],silent=T)!=2){
	# clean up first
	source('cleanup.R')
	cat('Setting up directories...')	
	# create the tmpfsDir and link to workerDirs
	# include workunit.i in the filenames so multiple instances can work side by side
	if(exists('workUnit.i')){
		tmpfsDir <- paste0(origTmpfsDir,'-workUnit-',workUnit.i)
		workDirLocation.frida <- paste0(baselocation.frida,'-',workUnit.i)
		workDirLocation.stella <- paste0(baselocation.stella,'-',workUnit.i)
		name.workDir <- paste0(origName.workDir,'-',workUnit.i)
	} else {
		workDirLocation.frida <- baselocation.frida
		workDirLocation.stella <- baselocation.stella
	}
	dir.create(tmpfsDir,recursive = T,showWarnings = F)
	system(paste('ln -s',tmpfsDir,name.workDir))
	system(paste('cp -r',baselocation.frida,file.path(tmpfsDir,workDirLocation.frida)))
	system(paste('cp -r',baselocation.stella,file.path(tmpfsDir,workDirLocation.stella)))
	location.frida <- paste0(baselocation.frida,'-',name.output)
	location.stella <- paste0(baselocation.stella,'-',name.output)
	if(exists('workUnit.i')){
		location.frida <- paste0(location.frida,'-',workUnit.i)
		location.stella <- paste0(location.stella,'-',workUnit.i)
	}
	system(paste('ln -s',file.path(tmpfsDir,workDirLocation.frida),location.frida))
	system(paste('ln -s',file.path(tmpfsDir,workDirLocation.stella),location.stella))
	cat('done\n')
} else {
	cat('Using existing directories\n')
}

if(disk.free(location.frida)< 2e4){
	stop('less than 20mib in frida location\n')
}
