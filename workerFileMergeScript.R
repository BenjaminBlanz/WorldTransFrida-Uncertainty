source('funRunFRIDA.R')

args <- commandArgs(T)
v.i <- as.numeric(args[1])
outputType <- as.character(args[2])
outputTypeFolder <- as.character(args[3])
varNames <- readRDS(as.character(args[4]))
verbosity<-as.numeric(args[5])
compressCsv<-as.logical(args[6])

varName <- varNames[v.i]
perVarSubfolder <- file.path(outputTypeFolder,varName)
fileList <- list.files(perVarSubfolder)
if(verbosity>0){cat(sprintf('Processing %i files of %s...',length(fileList),varName))}
if(verbosity>0){cat('reading and merging...')}

filesContents.lst <- list()
filesContents.lst[[1]] <- readPerVarFile(file.path(perVarSubfolder,fileList[1]),outputType)
for(f.i in 2:length(fileList)){
	filesContents.lst[[f.i]] <- readPerVarFile(file.path(perVarSubfolder,fileList[f.i]),outputType)
	# hack to deal with incomplete runs messing up column headers
	# proper fix is in data generation, but this will make old results work
	colnames(filesContents.lst[[f.i]]) <- colnames(filesContents.lst[[1]])
}
rbindStr <- paste0('varData <- base::rbind(filesContents.lst[[',
									 paste(1:length(filesContents.lst),
									 			collapse = ']],filesContents.lst[['),
									 ']])')
eval(parse(text=rbindStr))
varData <- sort_by(varData,varData[,1])
colnames(varData) <- gsub('(^X)([0-9]{4})','\\2',colnames(varData),perl = T)
if(verbosity>0){cat('writing...')}
writePerVarFile(varData,file.path(outputTypeFolder,varName),
								outputType=outputType,compressCsv=compressCsv)
if(verbosity>0){cat('removing split files...')}
unlink(perVarSubfolder,recursive = T,force = T)
if(verbosity>0){cat('done\n')}
rm(list=c('varData','filesContents.lst'))
