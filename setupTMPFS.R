
# if there is a cluster running, we must already have setup the TMPFS,
# so do nothing further.
if(!exists('cl') &&
	 try(clusterEvalQ(cl,1+1)[[1]],silent=T)!=2){
	# clean up first
	source('cleanup.R')
	cat('Setting up directories...')	
	# create the tmpfsDir and link to workerDirs
	dir.create(tmpfsDir,recursive = T,showWarnings = F)
	system(paste('ln -s',tmpfsDir,name.workDir))
	# set up the tmpfs for the single threaded runs
	system(paste('cp -r',baselocation.frida,file.path(tmpfsDir,baselocation.frida)))
	location.frida <- paste0(baselocation.frida,'-',name.output)
	system(paste('ln -s',file.path(tmpfsDir,baselocation.frida),location.frida))
	system(paste('cp -r',baselocation.stella,file.path(tmpfsDir,baselocation.stella)))
	location.stella <- paste0(baselocation.stella,'-',name.output)
	system(paste('ln -s',file.path(tmpfsDir,baselocation.stella),location.stella))
	cat('done\n')
} else {
	cat('Using existing directories\n')
}
