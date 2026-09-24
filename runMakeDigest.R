# runMakeDigest.R ####
#
# Writes a digest of result folders, for Zenodo uploads and for the code behind
# paper figures: the metadata, parameter space, sample points, representative
# sample and plot data of a run, without the per var files and the figures.
# A 10k sample run of 4.6 GB gives a digest of 85 MB, 54 MB of it sample points.
#
# The digest holds
#   every top level file of the run up to 50 MB, samplePoints only as
#   samplePoints.csv.gz
#   repSample/
#   figures/ without the png and pdf files, i.e. the plot data and the config of
#   the plotted runs; --figures copies all of figures/
#   digest.txt, naming the source, what was left out, and the md5 of each file
#
# Run from the repository root:
#   Rscript runMakeDigest.R [--figures] [--overwrite] [--out <dir>] <run or folder of runs> [...]
# The digest of a run goes next to it, or into --out. An existing digest is only
# replaced with --overwrite.

args <- commandArgs(trailingOnly=TRUE)
withFigures <- '--figures'%in%args
overwrite <- '--overwrite'%in%args
outDir <- NA
if('--out'%in%args){
	out.i <- match('--out',args)
	outDir <- args[out.i+1]
	args <- args[-c(out.i,out.i+1)]
}
roots <- args[!args%in%c('--figures','--overwrite')]
if(length(roots)==0||any(is.na(roots))){
	stop('usage: Rscript runMakeDigest.R [--figures] [--overwrite] [--out <dir>] <run or folder of runs> [...]\n')
}
roots <- normalizePath(roots,mustWork=TRUE)
if(!is.na(outDir)){
	dir.create(outDir,F,T)
	outDir <- normalizePath(outDir)
}
maxFileSize <- 50*1024^2

# runs ####
# a run folder has run metadata, figures or detectedParmSpace in it. Runs may be
# grouped in folders below the root.
isRun <- function(d){
	file.exists(file.path(d,'runMetadata.txt'))||
		dir.exists(file.path(d,'figures'))||dir.exists(file.path(d,'detectedParmSpace'))
}
findRuns <- function(d,depth){
	if(isRun(d)){
		return(d)
	}
	sub <- list.dirs(d,recursive=FALSE)
	sub <- sub[basename(sub)!='overlayed'&!grepl('-digest(\\.partial)?$',sub)]
	unlist(lapply(sub,function(s){
		if(isRun(s)){
			s
		} else if(depth>1){
			findRuns(s,depth-1)
		}
	}))
}

# digest ####
fmtSize <- function(bytes){
	units <- c('B','KB','MB','GB','TB')
	p <- pmax(0,pmin(length(units)-1,floor(log(pmax(bytes,1),1024))))
	sprintf('%.1f %s',bytes/1024^p,units[p+1])
}
dirSize <- function(d){
	sum(file.size(list.files(d,recursive=TRUE,full.names=TRUE,all.files=TRUE)),na.rm=TRUE)
}

makeDigest <- function(run){
	target <- paste0(if(is.na(outDir)) run else file.path(outDir,basename(run)),'-digest')
	if(dir.exists(target)){
		if(!overwrite){
			cat(sprintf('%s exists, skipping it (--overwrite replaces it)\n',target))
			return(invisible(NULL))
		}
		unlink(target,recursive=TRUE)
	}
	# built under a temporary name, so an interrupted digest is not taken for a complete one
	partial <- paste0(target,'.partial')
	unlink(partial,recursive=TRUE)
	dir.create(partial,F,T)
	leftOut <- data.frame(path=character(0),size=numeric(0))
	copyRel <- function(rel){
		dir.create(dirname(file.path(partial,rel)),F,T)
		file.copy(file.path(run,rel),file.path(partial,rel),copy.date=TRUE)
	}

	entries <- list.files(run,all.files=TRUE,no..=TRUE)
	isDir <- dir.exists(file.path(run,entries))

	# top level files
	files <- entries[!isDir&!entries%in%c('samplePoints.csv','samplePoints.RDS')]
	sizes <- file.size(file.path(run,files))
	for(f in files[sizes<=maxFileSize]){
		copyRel(f)
	}
	leftOut <- rbind(leftOut,data.frame(path=files[sizes>maxFileSize],size=sizes[sizes>maxFileSize]))

	# sample points
	if(file.exists(file.path(run,'samplePoints.csv'))){
		R.utils::gzip(file.path(run,'samplePoints.csv'),file.path(partial,'samplePoints.csv.gz'),
									remove=FALSE)
	} else if(file.exists(file.path(run,'samplePoints.RDS'))){
		write.csv(readRDS(file.path(run,'samplePoints.RDS')),gzfile(file.path(partial,'samplePoints.csv.gz')))
	}

	# repSample and figures
	if(dir.exists(file.path(run,'repSample'))){
		for(f in list.files(file.path(run,'repSample'),recursive=TRUE,all.files=TRUE)){
			copyRel(file.path('repSample',f))
		}
	}
	if(dir.exists(file.path(run,'figures'))){
		figFiles <- file.path('figures',list.files(file.path(run,'figures'),recursive=TRUE,all.files=TRUE))
		isImage <- grepl('\\.(png|pdf)$',figFiles,ignore.case=TRUE)
		if(withFigures){
			isImage[] <- FALSE
		}
		for(f in figFiles[!isImage]){
			copyRel(f)
		}
		if(any(isImage)){
			leftOut <- rbind(leftOut,data.frame(path=sprintf('figures/ png and pdf files (%i)',sum(isImage)),
																					size=sum(file.size(file.path(run,figFiles[isImage])))))
		}
	}

	# other folders
	for(d in entries[isDir&!entries%in%c('repSample','figures')]){
		leftOut <- rbind(leftOut,data.frame(path=paste0(d,'/'),size=dirSize(file.path(run,d))))
	}

	# digest.txt
	digestFiles <- sort(list.files(partial,recursive=TRUE,all.files=TRUE))
	md5 <- tools::md5sum(file.path(partial,digestFiles))
	writeLines(c('Digest of a FRIDA uncertainty analysis run',
							 '==========================================',
							 '',
							 sprintf('source       %s',run),
							 sprintf('host         %s',Sys.info()[['nodename']]),
							 sprintf('written      %s',format(Sys.time(),'%Y-%m-%d %H:%M:%S')),
							 sprintf('figures      %s',if(withFigures) 'all' else 'plot data only, no png or pdf'),
							 '',
							 'samplePoints.csv.gz is samplePoints.csv of the run, gzipped. runMetadata.txt',
							 'describes the model and scripts the run used.',
							 '',
							 'left out',
							 '--------',
							 if(nrow(leftOut)>0) sprintf('%-12s %s',fmtSize(leftOut$size),leftOut$path) else 'nothing',
							 '',
							 'md5',
							 '---',
							 sprintf('%s  %s',unname(md5),digestFiles)),
						 file.path(partial,'digest.txt'))
	file.rename(partial,target)
	cat(sprintf('%s: %s -> %s, %i item(s) left out\n  %s\n',basename(run),fmtSize(dirSize(run)),
							fmtSize(dirSize(target)),nrow(leftOut),target))
}

runs <- unique(unlist(lapply(roots,findRuns,depth=3)))
if(length(runs)==0){
	stop('no runs found\n')
}
cat(sprintf('writing the digest of %i run(s)\n',length(runs)))
for(run in runs){
	makeDigest(run)
}
