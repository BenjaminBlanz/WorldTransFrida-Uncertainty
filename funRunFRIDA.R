#
# Functions to run frida
#

# fun write frida input ####
# uses location.frida and name.fridaInputFile from the global env.
writeFRIDAInput <- function(variables,values){
	parmValues <- data.frame(Variable=variables,Value=values)
	write.table(parmValues,file = file.path(location.frida,'Data',name.fridaInputFile),
							row.names = F,col.names = F,sep=',')
}

# fun runFridaParmsByIndex ####
# Uses sampleParms,samplePoints,location.frida, and name.fridaInputFile
# from the global environment
runFridaParmsByIndex <- function(i){
	writeFRIDAInput(sampleParms$Variable,samplePoints[,i])
	system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',
							 file.path(location.frida,'FRIDA.stmx')),
				 ignore.stdout = T,ignore.stderr = T,wait = T)
	return(read.csv(file.path(location.frida,'Data',name.fridaOutputFile)))
}

# fun runFridaDefaultParms ####
# Uses sampleParms,samplePoints,location.frida, and name.fridaInputFile
# from the global environment
runFridaDefaultParms <- function(){
	system(paste('rm',file.path(location.frida,'data',name.fridaInputFile)),
				 ignore.stdout = T, ignore.stderr = T)
	system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',
							 file.path(location.frida,'FRIDA.stmx')),
				 ignore.stdout = T,ignore.stderr = T,wait = T)
	return(read.csv(file.path(location.frida,'Data',name.fridaOutputFile)))
}

