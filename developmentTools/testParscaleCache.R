# testParscaleCache.R ####
#
# The parscale cache has to distinguish a parameter that was never tried from one
# that was tried and could not be determined, has to keep working with a bare
# numeric parscale.RDS, and has to match cached entries by name, not by position.
#
# Exercises the load logic on synthetic files, so it needs neither stella nor
# a determination. Run from the repository root:
#   Rscript developmentTools/testParscaleCache.R

tmp <- tempfile('parscaleCacheTest')
dir.create(tmp)
on.exit(unlink(tmp,recursive=TRUE))

source('funParmSpace.R')

# the real load, driven the way runMLEandParmSpace.R drives it
loadParscaleCache <- function(jParVectNames,location.output,
															redoAllCalc=FALSE,redoFailedParscales=FALSE){
	if(redoAllCalc){
		parscale <- rep(NA_real_,length(jParVectNames))
		names(parscale) <- jParVectNames
		notDetermined <- rep(FALSE,length(jParVectNames))
		names(notDetermined) <- jParVectNames
		return(list(parscale=parscale,notDetermined=notDetermined))
	}
	funReadCachedParscale(file.path(location.output,'parscale.RDS'),
												jParVectNames,redoFailedParscales=redoFailedParscales)
}

# what the determination loop would pick up as work to do
workToDo <- function(res,parscaleSkip){
	which((is.na(res$parscale)|is.infinite(res$parscale))&
					!parscaleSkip&!res$notDetermined)
}

nms <- c('a','b','c','d')
skip <- setNames(c(FALSE,FALSE,FALSE,TRUE),nms)
ok <- 0
fail <- 0
check <- function(label,got,want){
	if(isTRUE(all.equal(got,want))){
		cat(sprintf('  ok   %s\n',label)); ok <<- ok+1
	} else {
		cat(sprintf('  FAIL %s\n       got  %s\n       want %s\n',label,
								paste(format(got),collapse=','),paste(format(want),collapse=',')))
		fail <<- fail+1
	}
}

cat('new format, one determined, one failed, one never tried\n')
saveRDS(list(parscale=setNames(c(1e-3,NA,NA,NA),nms),
						 status=setNames(c('determined','notDetermined',NA,'skippedExternalRange'),nms),
						 savedAt=Sys.time()),
				file.path(tmp,'parscale.RDS'))
r <- loadParscaleCache(nms,tmp)
check('only the never-tried parameter is work',unname(workToDo(r,skip)),3L)
check('the failed one is remembered as failed',unname(r$notDetermined),c(F,T,F,F))
check('the determined value is loaded',unname(r$parscale[1]),1e-3)

cat('redoFailedParscales retries the failure\n')
r <- loadParscaleCache(nms,tmp,redoFailedParscales=TRUE)
check('failed and never-tried are both work',unname(workToDo(r,skip)),c(2L,3L))

cat('redoAllCalc ignores the cache entirely\n')
r <- loadParscaleCache(nms,tmp,redoAllCalc=TRUE)
check('everything unskipped is work',unname(workToDo(r,skip)),c(1L,2L,3L))

cat('old bare-numeric format\n')
saveRDS(setNames(c(1e-3,NA,NA,NA),nms),file.path(tmp,'parscale.RDS'))
r <- loadParscaleCache(nms,tmp)
check('an NA without a status stays work',unname(workToDo(r,skip)),c(2L,3L))
check('nothing is claimed as a remembered failure',unname(r$notDetermined),c(F,F,F,F))

cat('cache holding the same parameters in a different order\n')
saveRDS(list(parscale=setNames(c(4,3,2,1),c('d','c','b','a')),
						 status=setNames(rep('determined',4),c('d','c','b','a')),
						 savedAt=Sys.time()),
				file.path(tmp,'parscale.RDS'))
r <- loadParscaleCache(nms,tmp)
check('values follow names, not positions',unname(r$parscale),c(1,2,3,4))

cat('cache missing a parameter that exists now\n')
saveRDS(list(parscale=setNames(c(1,2),c('a','b')),
						 status=setNames(c('determined','notDetermined'),c('a','b')),
						 savedAt=Sys.time()),
				file.path(tmp,'parscale.RDS'))
r <- loadParscaleCache(nms,tmp)
check('the new parameter is the only work',unname(workToDo(r,skip)),3L)

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('parscale cache load does not behave as intended\n')
}
