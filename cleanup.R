# cluster cleanup ####
cat('cluster cleanup...')
# stop cluster
tryCatch(stopCluster(cl),error=function(e){})
# clean up working directories
unlink('workerDirs',recursive=T,force=T)
unlink(tmpfsDir,recursive=T,force=T)
cat('done\n')

# single treaded tmpfs cleanup ####
if(file.exists(location.frida.storage)){
	system(paste('rm', location.frida))
	system(paste('mv',location.frida.storage,location.frida))
}
if(file.exists(location.stella.storage)){
	system(paste('rm', location.stella))
	system(paste('mv',location.stella.storage,location.stella))
}
