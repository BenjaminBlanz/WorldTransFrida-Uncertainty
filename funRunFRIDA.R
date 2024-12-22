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

# taken from funchir package
# convert plot region sizes to inches
xdev2in <- function (x = 1) 
{
	x * par("pin")[1L]/diff(par("usr")[1L:2L])
}
ydev2in <- function (y = 1) 
{
	y * par("pin")[2L]/diff(par("usr")[3L:4L])
}
xydev2in <-	function (xy = 1) 
{
	u = par("usr")
	xy * par("pin")/c(u[2L] - u[1L], u[4L] - u[3L])
}
dist.f <- function(oi,yi,ys,offsets,keepOutSize){
	return(min(sqrt((ydev2in(ys-yi))^2+(xdev2in(offsets-oi))^2))-keepOutSize)
}

# funValidRange ####
# returns the first and last index in the variable that has a data point
funValidRange <- function(x){
	validRange <- c(1,length(x))
	if(is.na(x[1])){
		validRange[1] <- which(diff(is.na(x))==-1)[1]+1
	}
	if(is.na(x[length(x)])){
		validRange[2] <- max(which(diff(is.na(x))==1))
	}
	return(validRange)
}

# funLikelihood ####
funLikelihood <- function(){
	
}
