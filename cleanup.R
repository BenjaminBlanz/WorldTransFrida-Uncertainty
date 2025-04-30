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
# use %in% list.files instead of file.exists as file.exists will return FALSE for dangling links
if(basename(location.frida)%in%list.files(dirname(baselocation.frida))){
	system(paste('rm', location.frida))
	location.frida <- baselocation.frida
}
location.stella <- paste0(baselocation.stella,'-',name.output)
if(basename(location.stella)%in%list.files(dirname(location.stella))){
	system(paste('rm', location.stella))
	location.stella <- baselocation.stella
}
# clean up working directories
if(basename(name.workDir)%in%list.files(dirname(name.workDir))){
	system(paste('rm -r',name.workDir))
}
if(file.exists(tmpfsDir)){
	system(paste('rm -r',tmpfsDir))
}
# only remove the tmpfsBaseDir if it is empty (other parallel runs may still have content)
if(file.exists(tmpfsBaseDir)&&length(list.files(tmpfsBaseDir))==0){
	system(paste('rm -r',tmpfsBaseDir))
}
cat('done\n')
