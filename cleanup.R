# cluster cleanup ####
# stop cluster
cat('ensuring cluster is not running...')
if(exists('cl')){
	tryCatch(stopCluster(cl),error=function(e){})
	rm(cl)
}
cat('done\n')
# directory cleanup ####
cat('cleaning up directories...')
# single threaded tmpfs cleanup ####
location.frida <- paste0(baselocation.frida,'-',name.output)
location.stella <- paste0(baselocation.stella,'-',name.output)
if(exists('workUnit.i')){
	location.frida <- paste0(location.frida,'-',workUnit.i)
	location.stella <- paste0(location.stella,'-',workUnit.i)
}
# use %in% list.files instead of file.exists as file.exists will return FALSE for dangling links
if(basename(location.frida)%in%list.files(dirname(baselocation.frida))){
	system(paste('rm -rf', location.frida))
	location.frida <- baselocation.frida
}
if(basename(location.stella)%in%list.files(dirname(location.stella))){
	system(paste('rm -rf', location.stella))
	location.stella <- baselocation.stella
}
# clean up working directories
if(exists('workUnit.i')){
	name.workDir <- paste0(origName.workDir,'-',workUnit.i)
}
if(basename(name.workDir)%in%list.files(dirname(name.workDir))){
	system(paste('rm -rf',name.workDir))
}
if(file.exists(tmpfsDir)){
	system(paste('rm -rf',tmpfsDir))
}
# only remove the tmpfsBaseDir if it is empty (other parallel runs may still have content)
if(file.exists(tmpfsBaseDir)&&length(list.files(tmpfsBaseDir))==0){
	system(paste('rm -rf',tmpfsBaseDir))
}
cat('done\n')
