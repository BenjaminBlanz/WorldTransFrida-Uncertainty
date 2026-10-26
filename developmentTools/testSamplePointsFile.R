# testSamplePointsFile.R ####
#
# A csv of sample points, e.g. a Saltelli design, is run as given. Its rows are
# the runs, numbered in file order, its optional first column Case labels them,
# and the frida_info file of the run is built from its columns. Anything that
# would run other points than the file holds has to stop the run instead.
#
# Run from the repository root:
#   Rscript developmentTools/testSamplePointsFile.R

library(data.table)
source('funRunFRIDA.R')
source('funParmSpace.R')

tmp <- tempfile('samplePointsFileTest')
dir.create(tmp)
on.exit(unlink(tmp,recursive=TRUE))

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),''))); fail <<- fail+1}
}
# the error message of expr, NA if it ran through
errorOf <- function(expr){
	tryCatch({force(expr); NA_character_},error=function(e){conditionMessage(e)})
}
checkError <- function(label,expr,pattern){
	msg <- errorOf(expr)
	check(label,!is.na(msg)&&grepl(pattern,msg),ifelse(is.na(msg),'no error',msg))
}
csvFile <- function(name,lines){
	f <- file.path(tmp,name)
	writeLines(lines,f)
	return(f)
}

cat('reading the csv\n')
f <- csvFile('withCase.csv',c('Case,p.a[1],p.b[1]','A:1,1,2','B:1,3,4','A_B[p.a[1]]:1,3,2'))
r <- funReadSamplePointsFile(f)
check('one row per run, one column per parameter',identical(dim(r$samplePoints),c(3L,2L)))
check('the parameter names are kept as written',
			identical(colnames(r$samplePoints),c('p.a[1]','p.b[1]')),
			paste(colnames(r$samplePoints),collapse=', '))
check('the runs are numbered in file order',identical(rownames(r$samplePoints),c('1','2','3')))
check('the Case column labels them',identical(r$cases,c('A:1','B:1','A_B[p.a[1]]:1')))
check('the values are numbers',is.numeric(r$samplePoints)&&r$samplePoints[3,'p.a[1]']==3)

fGz <- file.path(tmp,'withCase.csv.gz')
R.utils::gzip(f,fGz,remove=FALSE)
check('a gzipped csv reads the same',identical(funReadSamplePointsFile(fGz),r))

rNoCase <- funReadSamplePointsFile(csvFile('noCase.csv',c('p.a[1],p.b[1]','1,2','3,4')))
check('without a Case column there are no labels',is.null(rNoCase$cases))
check('and every column is a parameter',identical(colnames(rNoCase$samplePoints),c('p.a[1]','p.b[1]')))

checkError('a Case used twice stops, naming the row',
					 funReadSamplePointsFile(csvFile('dupCase.csv',c('Case,p.a[1]','A:1,1','A:1,2'))),
					 'Case.*row 2')
checkError('an empty Case stops, naming the row',
					 funReadSamplePointsFile(csvFile('emptyCase.csv',c('Case,p.a[1]','A:1,1',',2'))),
					 'Case.*row 2')
checkError('a value that is not a number stops, naming the column',
					 funReadSamplePointsFile(csvFile('text.csv',c('Case,p.a[1],p.b[1]','A:1,1,x','B:1,2,3'))),
					 'p\\.b\\[1\\]')
checkError('a missing value stops, naming the column',
					 funReadSamplePointsFile(csvFile('na.csv',c('Case,p.a[1],p.b[1]','A:1,1,','B:1,2,3'))),
					 'p\\.b\\[1\\]')
big <- funReadSamplePointsFile(csvFile('big.csv',c('p.a[1],p.b[1]','3000000000,1','3000000001,2')))
check('whole numbers beyond the integer range keep their value',
			big$samplePoints[1,'p.a[1]']==3e9,big$samplePoints[1,'p.a[1]'])
numCase <- funReadSamplePointsFile(csvFile('numCase.csv',c('Case,p.a[1]','1.10,1','1.1,2')))
check('a Case that looks like a number is kept as written',identical(numCase$cases,c('1.10','1.1')))
checkError('an infinite value stops, naming the column',
					 funReadSamplePointsFile(csvFile('inf.csv',c('Case,p.a[1],p.b[1]','A:1,1,Inf','B:1,2,3'))),
					 'p\\.b\\[1\\]')
checkError('a parameter given twice stops, naming it',
					 funReadSamplePointsFile(csvFile('twice.csv',c('Case,p.a[1],p.a[1]','A:1,1,2','B:1,2,3'))),
					 'p\\.a\\[1\\]')

cat('the frida_info of the sample points\n')
info <- data.frame(Variable=c('p.a[1]','p.b[1]','p.c[1]'),Value=c(2,5,1))
printed <- capture.output(p <- funSamplePointsParms(r$samplePoints,info))
check('one row per column, in column order',identical(p$Variable,c('p.a[1]','p.b[1]')))
check('Value is the model default',identical(p$Value,c(2,5)))
check('Min and Max span the sampled values',identical(p$Min,c(1,2))&&p$Max[1]==3)
check('a default outside them widens the range to include it',p$Max[2]==5,p$Max[2])
check('and that is printed',any(grepl('p.b[1]',printed,fixed=TRUE))&&
				!any(grepl('p.a[1]',printed,fixed=TRUE)),paste(printed,collapse='\n'))
checkError('a parameter the model does not have stops, naming it',
					 funSamplePointsParms(r$samplePoints,info[-1,]),'p\\.a\\[1\\]')
constant <- r$samplePoints
constant[,'p.b[1]'] <- 2
checkError('a parameter with a single value stops, naming it',
					 funSamplePointsParms(constant,info),'p\\.b\\[1\\]')

cat('the run loads the supplied points\n')
out <- file.path(tmp,'output')
dir.create(out)
sampleParms <- data.frame(Variable=c('p.b[1]','p.a[1]'),Value=c(2,2),Min=c(0,0),Max=c(5,5))
checkError('without suppliedSamplePoints.RDS the run stops, saying how to make it',
					 funLoadSuppliedSamplePoints(out,sampleParms,3),'runPrepareSamplePoints')
saveRDS(r$samplePoints,file.path(out,'suppliedSamplePoints.RDS'))
sp <- funLoadSuppliedSamplePoints(out,sampleParms,3)
check('the columns follow the sampled parameters',identical(colnames(sp),c('p.b[1]','p.a[1]')))
check('with their values',sp[1,'p.b[1]']==2&&sp[3,'p.a[1]']==3)
check('and the runs keep their numbers',identical(rownames(sp),c('1','2','3')))
checkError('a different number of runs than numSample stops',
					 funLoadSuppliedSamplePoints(out,sampleParms,4),'3.*4|4.*3')
otherParms <- sampleParms
otherParms$Variable[2] <- 'p.c[1]'
checkError('other parameters than the sampled ones stop, naming both sides',
					 funLoadSuppliedSamplePoints(out,otherParms,3),'p\\.c\\[1\\].*p\\.a\\[1\\]|p\\.a\\[1\\].*p\\.c\\[1\\]')

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('the sample points file is not read as it should be\n')
}
