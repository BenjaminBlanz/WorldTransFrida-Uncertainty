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
	retDat <- read.csv(file.path(location.frida,'Data',name.fridaOutputFile))
	colnames(retDat) <- cleanNames(colnames(retDat))
	rownames(retDat) <- retDat$year
	retDat <- retDat[,-1]
	return(retDat)
}

# fun runFridaDefaultParms ####
# Uses sampleParms,samplePoints,location.frida, and name.fridaInputFile
# from the global environment
runFridaDefaultParms <- function(){
	system(paste('rm',file.path(location.frida,'Data',name.fridaInputFile)),
				 ignore.stdout = T, ignore.stderr = T)
	system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',
							 file.path(location.frida,'FRIDA.stmx')),
				 ignore.stdout = T,ignore.stderr = T,wait = T)
	retDat <- read.csv(file.path(location.frida,'Data',name.fridaOutputFile))
	colnames(retDat) <- cleanNames(colnames(retDat))
	rownames(retDat) <- retDat$year
	retDat <- retDat[,-1]
	return(retDat)
}

# fun cleanNames ####
# takes a vector of e.g. column names and brings them into 
# a comparable standard format
# also drops the trailing 1 of the run id which we do not use
# as we only have single run setups of frida
cleanNames <- function(colNames){
	gsub('time','year',
			 gsub('_$','',
			 		 gsub('_1','',
			 		 		 gsub('\\[1]','_',
			 		 		 		 gsub('_+','_',
			 		 		 		 		 gsub('[. ]','_',
			 		 		 		 		 		 tolower(colNames)))))))
}

# fun idxOfVarName ####
idxOfVarName <- function(varNames,vecOfVarNames){
	varNames <- cleanNames(varNames)
	vecOfVarNames <- cleanNames(vecOfVarNames)
	idx <- c()
	for(i in 1:length(varNames)){
		idx[i] <- which(vecOfVarNames==varNames[i])
	}
	return(idx)
}
