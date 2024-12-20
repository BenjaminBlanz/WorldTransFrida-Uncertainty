#
# Functions to run frida
#

# fun write frida input ####
writeFRIDAInput <- function(variables,values,location.frida,name.fridaInputFile){
	parmValues <- data.frame(Variable=variables,Value=values)
	write.table(parmValues,file = file.path(location.frida,'Data',name.fridaInputFile),
							row.names = F,col.names = F,sep=',')
}

# fun runFridaOnceByIndex ####
# Uses sampleParms,samplePoints,location.frida, and name.fridaInputFile
# from the global environment!
runFridaParmsByIndex <- function(i){
	writeFRIDAInput(sampleParms$Variable,samplePoints[,i],location.frida,name.fridaInputFile)
	system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',
							 file.path(location.frida,'FRIDA.stmx')),
				 ignore.stdout = T,ignore.stderr = T,wait = T)
	return(read.csv(file.path(location.frida,'Data',name.fridaOutputFile)))
}

