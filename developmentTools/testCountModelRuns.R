# testCountModelRuns.R ####
#
# The run counter is the instrument every other measurement in this work depends
# on, so it needs checking like anything else: that it is silent until switched
# on, that it counts every call once, and that it attributes calls to the section
# that was open when they happened.
#
# Run from the repository root:
#   Rscript developmentTools/testCountModelRuns.R

source('developmentTools/countModelRuns.R')

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),'')))
		fail <<- fail+1}
}

# the hook exactly as it appears in runFRIDASpecParms
hookedRun <- function(){
	if(exists('devTools.countModelRuns')&&isTRUE(devTools.countModelRuns)){
		devTools.modelRunCount <<- if(exists('devTools.modelRunCount')){
			devTools.modelRunCount+1
		} else {
			1
		}
	}
	invisible(NULL)
}

cat('the hook is inert until it is switched on\n')
if(exists('devTools.countModelRuns')){rm(devTools.countModelRuns)}
if(exists('devTools.modelRunCount')){rm(devTools.modelRunCount)}
for(i in 1:5){hookedRun()}
check('nothing is counted and no variable is created',!exists('devTools.modelRunCount'))

cat('once switched on it counts every call\n')
devToolsCountModelRuns()
for(i in 1:7){hookedRun()}
check('seven calls give seven',devToolsModelRunCount()==7,
			sprintf('got %d',devToolsModelRunCount()))

cat('sections get the calls made while they were open\n')
devTools.sectionCounts <- list()
devTools.currentSection <- NULL
devToolsMarkSection('parscale')
for(i in 1:11){hookedRun()}
devToolsMarkSection('range finding')
for(i in 1:4){hookedRun()}
counts <- devToolsReportModelRuns(file=file.path(tempdir(),'counts.csv'))
check('the first section got its eleven',counts[['parscale']]==11,
			sprintf('got %d',counts[['parscale']]))
check('the second got its four',counts[['range finding']]==4,
			sprintf('got %d',counts[['range finding']]))
check('and the calls before any section were not attributed to one',
			sum(counts)==15,sprintf('total %d',sum(counts)))

cat('a section opened twice accumulates\n')
devTools.sectionCounts <- list()
devTools.currentSection <- NULL
devToolsMarkSection('parscale')
for(i in 1:3){hookedRun()}
devToolsMarkSection('other')
hookedRun()
devToolsMarkSection('parscale')
for(i in 1:2){hookedRun()}
counts <- devToolsReportModelRuns(file=file.path(tempdir(),'counts2.csv'))
check('both visits are added up',counts[['parscale']]==5,
			sprintf('got %d',counts[['parscale']]))

cat('two runs can be compared\n')
before <- file.path(tempdir(),'before.csv')
after <- file.path(tempdir(),'after.csv')
write.csv(data.frame(section=c('parscale','range finding'),stellaRuns=c(4000,1200)),
					before,row.names=FALSE)
write.csv(data.frame(section=c('parscale','range finding'),stellaRuns=c(1000,600)),
					after,row.names=FALSE)
cmp <- devToolsCompareModelRuns(before,after)
check('the saving per section is reported',
			isTRUE(all.equal(sort(cmp$saving),c(2,4))))

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('the model run counter does not count correctly\n')
}
