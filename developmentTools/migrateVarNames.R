# migrateVarNames.R ####
#
# Renames the per var files and plot files of existing runs, and the variable
# names inside their result files, to the names the current cleanNames gives.
# Before its rewrite cleanNames dropped every _1 (Aged 1 to 20 Years became
# aged_to_20_years) and turned every time into year (future time in recession
# became future_year_in_recession).
#
# A file's old and new name are both worked out from the FRIDA name of its
# variable, found in the export lists, the calibration data and the plot data of
# the runs. An old name that stands for more than one new name, or that is the
# correct name of another variable, is left alone and reported, as are per var
# files whose name matches no known FRIDA name.
#
# Run from the repository root, on the machine the results are on, while no job
# writes to them:
#   Rscript developmentTools/migrateVarNames.R <workOutput dir> [...]           reports
#   Rscript developmentTools/migrateVarNames.R --apply <workOutput dir> [...]   renames
# --apply logs every change to varNameMigration-<date>.tsv in each workOutput
# dir, which is what it takes to undo it.

source('funRunFRIDA.R')

args <- commandArgs(trailingOnly=TRUE)
applyChanges <- '--apply'%in%args
roots <- args[args!='--apply']
if(length(roots)==0){
	stop('usage: Rscript developmentTools/migrateVarNames.R [--apply] <workOutput dir> [...]\n')
}
roots <- normalizePath(roots,mustWork=TRUE)

# cleanNames before the rewrite, which named the existing files
cleanNamesBefore <- function(colNames){
	gsub('time','year',
			 gsub('_+$','',
			 		 gsub('_+','_',
			 		 		 gsub(',','_',
			 		 		 		 gsub('\\$','',
			 		 		 		 		 gsub('_1','',
			 		 		 		 		 		 gsub('\\]','_',
			 		 		 		 		 		 		 gsub('\\[\\*','_',
			 		 		 		 		 		 		 		 gsub('\\[\\d+','_',
			 		 		 		 		 		 		 		 		 gsub('[. ]','_',
			 		 		 		 		 		 		 		 		 		 tolower(colNames)))))))))))
}

# runs still writing would race the renames ####
if(applyChanges){
	jobs <- tryCatch(system2('squeue',c('-h','-u',Sys.getenv('USER'),'-o','%j'),
													 stdout=TRUE,stderr=FALSE),
									 error=function(e){character(0)})
	busy <- intersect(jobs,unlist(lapply(roots,function(r){basename(list.dirs(r,recursive=FALSE))})))
	if(length(busy)>0){
		stop(sprintf('jobs are running for %s, apply once they are done\n',paste(busy,collapse=', ')))
	}
}

# files named after a variable ####
perVarExt <- '\\.(csv\\.gz|csv|RDS)$'
plotDataSuffix <- '-[^-]+-[^-]+-weighted\\.(RDS|csv)$'
pngSuffix <- '-[^-]+-weighted\\.png$'
# the variable part of each path, NA for a path not named after a variable
varPartOf <- function(paths){
	base <- basename(paths)
	dir <- basename(dirname(paths))
	res <- rep(NA_character_,length(paths))
	perVar <- grepl('^PerVarFiles-',dir) & grepl(perVarExt,base)
	res[perVar] <- sub(perVarExt,'',base[perVar])
	chunks <- dir=='PerVarChunks'
	res[chunks] <- base[chunks]
	plotData <- dir=='plotData' & grepl(plotDataSuffix,base)
	res[plotData] <- sub(plotDataSuffix,'',base[plotData])
	png <- grepl('/CI-plots/',paths) & dir!='plotData' & grepl(pngSuffix,base)
	res[png] <- sub(pngSuffix,'',base[png])
	target <- dir=='subSample.TargetVars' & grepl('\\.png$',base)
	res[target] <- sub('\\.png$','',base[target])
	return(res)
}
kindOf <- function(paths){
	dir <- basename(dirname(paths))
	ifelse(grepl('^PerVarFiles-',dir),'per var file',
				 ifelse(dir=='PerVarChunks','chunk folder',
				 			 ifelse(dir=='plotData','plot data',
				 			 			 ifelse(dir=='subSample.TargetVars','target figure','figure'))))
}

cat('listing the runs...')
runs <- unlist(lapply(roots,function(r){list.dirs(r,recursive=FALSE)}))
files <- unlist(lapply(runs,function(run){
	dps <- file.path(run,'detectedParmSpace')
	perVarDirs <- list.dirs(dps,recursive=FALSE)
	perVarDirs <- perVarDirs[grepl('^PerVarFiles-',basename(perVarDirs))]
	c(unlist(lapply(perVarDirs,list.files,full.names=TRUE)),
		list.dirs(file.path(dps,'PerVarChunks'),recursive=FALSE),
		list.files(file.path(run,'figures'),recursive=TRUE,full.names=TRUE))
}))
files <- files[!is.na(varPartOf(files))]
varParts <- varPartOf(files)
cat(sprintf('done, %i runs, %i files named after a variable\n',length(runs),length(files)))

# FRIDA names ####
cat('collecting FRIDA names...')
readFQN <- function(f){
	fqn <- tryCatch(read.csv(f,check.names=FALSE)$FRIDA.FQN,error=function(e){NULL})
	fqn[!is.na(fqn)&nchar(fqn)>4]
}
exportLists <- c(file.path(dirname(roots),'FRIDA-info','frida_extra_variables_to_export_list.csv'),
								 unlist(lapply(runs,function(run){
								 	list.files(file.path(run,'runScriptsAndConfiguration','input','FRIDA-info'),
								 						 pattern='\\.csv$',full.names=TRUE)
								 })))
rawNames <- unlist(lapply(exportLists[file.exists(exportLists)],readFQN))
calFiles <- file.path(runs,'Calibratio_Data_Cleaned_and_Transposed.csv')
rawNames <- c(rawNames,unlist(lapply(calFiles[file.exists(calFiles)],function(f){
	tryCatch(names(read.csv(f,check.names=FALSE,nrows=1))[-1],error=function(e){NULL})
})))
stemsOf <- function(raw,fun){fun(gsub('\\[\\d+','',gsub(' \\[(\\d)\\]','[\\1]',gsub(' +$','',raw))))}
# a plotted variable carries its FRIDA name, one file each is enough for the
# variables no other source knows
known <- c(stemsOf(rawNames,cleanNamesBefore),stemsOf(rawNames,cleanNames))
unknownPlots <- files[kindOf(files)=='plot data' & grepl('\\.RDS$',files) & !varParts%in%known]
unknownPlots <- unknownPlots[!duplicated(varPartOf(unknownPlots))]
rawNames <- c(rawNames,unlist(lapply(unknownPlots,function(f){
	tryCatch(as.character(readRDS(f)$varName.orig)[1],error=function(e){NULL})
})))
rawNames <- unique(rawNames[!is.na(rawNames)])
cat(sprintf('done, %i\n',length(rawNames)))

# old name -> new name ####
pairs <- unique(data.frame(old=stemsOf(rawNames,cleanNamesBefore),new=stemsOf(rawNames,cleanNames),
													 raw=rawNames))
# an array entry names no file of its own, its elements do
pairs <- pairs[!grepl('*',pairs$new,fixed=TRUE),]
correct <- unique(pairs$new)
changes <- pairs[pairs$old!=pairs$new,]
newCount <- vapply(split(changes$new,changes$old),function(x){length(unique(x))},integer(1))
ambiguous <- names(newCount)[newCount>1]
clashing <- intersect(changes$old,correct)
skipped <- changes[changes$old%in%c(ambiguous,clashing),]
changes <- changes[!changes$old%in%c(ambiguous,clashing),]
newOf <- setNames(changes$new,changes$old)[!duplicated(changes$old)]

cat('\nnames to change\n')
for(o in names(newOf)){
	cat(sprintf('  %s -> %s\n',o,newOf[[o]]))
}
if(nrow(skipped)>0){
	cat('\nleft alone, the old name does not say which variable a file belongs to\n')
	for(o in unique(skipped$old)){
		cat(sprintf('  %s: %s\n',o,paste(unique(skipped$raw[skipped$old==o]),collapse=' | ')))
	}
}

# files ####
toRename <- files[varParts%in%names(newOf)]
parts <- varPartOf(toRename)
targets <- file.path(dirname(toRename),
										 paste0(newOf[parts],substring(basename(toRename),nchar(parts)+1)))
taken <- file.exists(targets)
rootOf <- function(path){roots[startsWith(path,paste0(roots,'/'))][1]}
cat('\nfiles to rename\n')
if(length(toRename)>0){
	print(table(checkout=basename(dirname(vapply(toRename,rootOf,character(1)))),
							kind=kindOf(toRename)))
} else {
	cat('  none\n')
}
if(any(taken)){
	cat(sprintf('\n%i of them left alone, a file of the new name exists already:\n',sum(taken)))
	cat(paste0('  ',toRename[taken],collapse='\n'),'\n')
}
unmapped <- unique(varParts[kindOf(files)=='per var file' &
													 	!varParts%in%c(known,correct,names(newOf),'logLike','runStatus')])
if(length(unmapped)>0){
	cat(sprintf('\n%i per var file names match no known FRIDA name, left as they are:\n',length(unmapped)))
	cat(paste0('  ',sort(unmapped),collapse='\n'),'\n')
}

# names inside the result files ####
# data frames, matrices and lists of them, with names from newOf
renameInside <- function(x){
	if(is.list(x)&&!is.data.frame(x)){
		return(lapply(x,renameInside))
	}
	dn <- dimnames(x)
	if(!is.null(dn)){
		dimnames(x) <- lapply(dn,function(n){
			if(is.null(n)){return(n)}
			hit <- n%in%names(newOf)
			n[hit] <- newOf[n[hit]]
			n
		})
	}
	return(x)
}
runFiles <- unlist(lapply(runs,function(run){file.path(run,c('calDat.RDS','sigma.RDS','sigma-indepParms.RDS'))}))
runFiles <- runFiles[file.exists(runFiles)]
toRewrite <- runFiles[vapply(runFiles,function(f){
	obj <- tryCatch(readRDS(f),error=function(e){NULL})
	!is.null(obj)&&!identical(renameInside(obj),obj)
},logical(1))]
cat(sprintf('\n%i run files hold old names and are rewritten, %i renamed plot data files get their $variable set\n',
						length(toRewrite),sum(kindOf(toRename)=='plot data'&grepl('\\.RDS$',toRename)&!taken)))

if(!applyChanges){
	cat('\nnothing changed, run with --apply to make these changes\n')
	quit(status=0)
}

# apply ####
logFiles <- setNames(file.path(roots,sprintf('varNameMigration-%s.tsv',format(Sys.Date()))),roots)
logChange <- function(action,from,to){
	logFile <- logFiles[[rootOf(from)]]
	if(!file.exists(logFile)){
		cat('action\tfrom\tto\n',file=logFile)
	}
	cat(sprintf('%s\t%s\t%s\n',action,from,to),file=logFile,append=TRUE)
}
cat('\nrenaming...')
for(i in which(!taken)){
	if(!file.rename(toRename[i],targets[i])){
		stop(sprintf('could not rename %s\n',toRename[i]))
	}
	logChange('rename',toRename[i],targets[i])
	if(kindOf(targets[i])=='plot data'&&grepl('\\.RDS$',targets[i])){
		p <- readRDS(targets[i])
		if(identical(p$variable,parts[i])){
			p$variable <- newOf[[parts[i]]]
			saveRDS(p,targets[i])
			logChange('set $variable',targets[i],newOf[[parts[i]]])
		}
	}
}
cat('done\nrewriting...')
for(f in toRewrite){
	saveRDS(renameInside(readRDS(f)),f)
	logChange('rename inside',f,f)
}
cat('done\n')
cat(sprintf('logged to %s\n',paste(logFiles[file.exists(logFiles)],collapse=', ')))
