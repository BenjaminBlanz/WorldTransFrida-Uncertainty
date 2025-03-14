# single threaded tmpfs cleanup ####
if(file.exists(location.frida.storage)){
	system(paste('rm', location.frida))
	system(paste('mv',location.frida.storage,location.frida))
}
if(file.exists(location.stella.storage)){
	system(paste('rm', location.stella))
	system(paste('mv',location.stella.storage,location.stella))
}


# cluster cleanup ####
cat('cluster cleanup...')
# stop cluster
if(exists('cl')){
	tryCatch(stopCluster(cl),error=function(e){})
	rm(cl)
}
# clean up working directories
unlink(name.workDir,recursive=T,force=T)
unlink(tmpfsBaseDir,recursive=T,force=T)
cat('done\n')
