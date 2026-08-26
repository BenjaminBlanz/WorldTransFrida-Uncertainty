
# the frida/stella locations of this run (and of this work unit). Determined here
# rather than only where the directories are created, because they have to be set on
# both paths below. config.R resets them to the base locations every time it is
# sourced, so when we skip the setup nothing else points them back at the directories
# belonging to this run, and every job started from this working directory would read
# and write the same shared frida directory.
runLocation.frida <- paste0(baselocation.frida,'-',name.output)
runLocation.stella <- paste0(baselocation.stella,'-',name.output)
if(exists('workUnit.i')){
	runLocation.frida <- paste0(runLocation.frida,'-',workUnit.i)
	runLocation.stella <- paste0(runLocation.stella,'-',workUnit.i)
}

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
	system(paste('ln -s',file.path(tmpfsDir,workDirLocation.frida),runLocation.frida))
	system(paste('ln -s',file.path(tmpfsDir,workDirLocation.stella),runLocation.stella))
	cat('done\n')
} else {
	cat('Using existing directories\n')
}
# after the branch, and after cleanup.R which resets these to the base locations
location.frida <- runLocation.frida
location.stella <- runLocation.stella

if(disk.free(location.frida)< 2e4){
	stop('less than 20mib in frida location\n')
}
