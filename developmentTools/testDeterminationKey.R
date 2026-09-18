# testDeterminationKey.R ####
#
# A cached determination is rejected unless it was computed from the same model,
# the same calibration data and the same likelihood settings, so that the
# parscales and ranges of one model are never handed to another in silence.
#
# Every field of the key, one at a time, has to reject the cache and say which
# one moved. The files the analysis writes into the model directory, and a new
# copy of the same model, must not.
#
# Run from the repository root:
#   Rscript developmentTools/testDeterminationKey.R

source('funRunFRIDA.R')
source('funParmSpace.R')

tmp <- tempfile('determinationKeyTest')
dir.create(tmp)
dir.create(file.path(tmp,'frida','FRIDA_Modules'),recursive=TRUE)
dir.create(file.path(tmp,'frida','Data'))
dir.create(file.path(tmp,'info'))
on.exit(unlink(tmp,recursive=TRUE))

writeLines('a model',con=file.path(tmp,'frida','FRIDA.stmx'))
writeLines('a module',con=file.path(tmp,'frida','FRIDA_Modules','Climate.itmx'))
writeLines('year,x',con=file.path(tmp,'frida','Data','frida_input_data.csv'))
writeLines('a policy',con=file.path(tmp,'frida','Data','policyParameters.csv'))
writeLines('a variable',con=file.path(tmp,'frida','Data','varsForExport.txt'))
exclude <- c('policyParameters.csv','varsForExport.txt')
writeLines('Variable,Value',con=file.path(tmp,'info','frida_info.csv'))

calDat <- data.frame(a=1:5,b=6:10)
resSigma <- diag(2)
parNames <- c('p1','p2','p3')
settings <- list(treatVarsAsIndep=TRUE,likeCutoffRatio=1000,rangeTol=1e-15,
								 ignoreParBounds=FALSE,forceParBounds=FALSE,
								 rangeRootTol=1e-4,rangeRootMaxIter=60)

buildKey <- function(calDat=get('calDat',envir=parent.frame()),
										 resSigma=get('resSigma',envir=parent.frame()),
										 parNames=get('parNames',envir=parent.frame()),
										 settings=get('settings',envir=parent.frame()),
										 baseNegLL=100,frida=file.path(tmp,'frida')){
	funDeterminationKey(frida,file.path(tmp,'info'),'frida_info.csv',
											calDat,resSigma,parNames,settings,baseNegLL=baseNegLL,
											exclude=exclude)
}

base <- buildKey()

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),'')))
		fail <<- fail+1}
}
# the mismatch must name what moved, not just that something did
namesIt <- function(mismatch,what){
	length(mismatch)>0 && any(grepl(what,mismatch,fixed=TRUE))
}

cat('an unchanged key is accepted\n')
check('no mismatch against itself',length(funDeterminationKeyMismatch(base,buildKey()))==0)
check('and the parameter list agrees',!funDeterminationParNameMismatch(base,buildKey()))

cat('every field is checked\n')
# puts the file back afterwards, so the checks that follow are not comparing
# against a file this test just changed
withChangedFile <- function(path,newContent){
	original <- readLines(path,warn=FALSE)
	writeLines(newContent,con=path)
	m <- funDeterminationKeyMismatch(base,buildKey())
	writeLines(original,con=path)
	m
}
m <- withChangedFile(file.path(tmp,'frida','FRIDA.stmx'),'a different model')
check('a changed model file is caught and named',namesIt(m,'FRIDA model'),
			paste(m,collapse='; '))
check('and the file is restored, so nothing after this sees it as changed',
			length(funDeterminationKeyMismatch(base,buildKey()))==0)

m <- withChangedFile(file.path(tmp,'frida','FRIDA_Modules','Climate.itmx'),'another module')
check('a changed module file is caught and named',namesIt(m,'FRIDA model'),
			paste(m,collapse='; '))

m <- withChangedFile(file.path(tmp,'frida','Data','frida_input_data.csv'),'year,y')
check('changed model input data is caught and named',namesIt(m,'FRIDA model'),
			paste(m,collapse='; '))

m <- withChangedFile(file.path(tmp,'info','frida_info.csv'),'Variable,Value,Extra')
check('changed frida_info is caught and named',namesIt(m,'frida_info'),
			paste(m,collapse='; '))

m <- funDeterminationKeyMismatch(base,buildKey(calDat=data.frame(a=1:5,b=11:15)))
check('changed calibration data is caught and named',namesIt(m,'calibration data'),
			paste(m,collapse='; '))

m <- funDeterminationKeyMismatch(base,buildKey(resSigma=diag(2)*2))
check('a changed residual covariance is caught and named',
			namesIt(m,'residual covariance'),paste(m,collapse='; '))

for(setting in names(settings)){
	changed <- settings
	changed[[setting]] <- if(is.logical(settings[[setting]])){
		!settings[[setting]]
	} else {
		settings[[setting]]*2
	}
	m <- funDeterminationKeyMismatch(base,buildKey(settings=changed))
	check(sprintf('a changed %s is caught and named',setting),namesIt(m,setting),
				paste(m,collapse='; '))
}

m <- funDeterminationKeyMismatch(base,buildKey(baseNegLL=100.5))
check('a moved starting likelihood is caught and named',
			namesIt(m,'likelihood at the starting parameters'),paste(m,collapse='; '))

cat('what the analysis writes and where the model sits do not count\n')
m <- withChangedFile(file.path(tmp,'frida','Data','policyParameters.csv'),'another policy')
check('a different scenario file is not a different model',length(m)==0,
			paste(m,collapse='; '))
m <- withChangedFile(file.path(tmp,'frida','Data','varsForExport.txt'),'another variable')
check('a different export list is not a different model',length(m)==0,
			paste(m,collapse='; '))
# every job copies the model afresh, which gives every file a new mtime
Sys.sleep(1.1)
dir.create(file.path(tmp,'copy'))
invisible(file.copy(file.path(tmp,'frida'),file.path(tmp,'copy'),recursive=TRUE,copy.date=FALSE))
m <- funDeterminationKeyMismatch(base,buildKey(frida=file.path(tmp,'copy','frida')))
check('a fresh copy of the same model is the same model',length(m)==0,
			paste(m,collapse='; '))

cat('the parameter list is reported apart from the rest\n')
newPars <- buildKey(parNames=c('p1','p2','p3','p4'))
check('a new parameter does not invalidate everything else',
			length(funDeterminationKeyMismatch(base,newPars))==0)
check('but is reported on its own',funDeterminationParNameMismatch(base,newPars))

cat('a cache with no key at all\n')
check('is rejected, saying so',
			identical(funDeterminationKeyMismatch(NULL,base),'no key was recorded with it'))
check('and counts as a parameter mismatch too',funDeterminationParNameMismatch(NULL,base))

cat('the readers act on it\n')
saveRDS(list(parscale=setNames(c(1,2,3),parNames),
						 status=setNames(rep('determined',3),parNames),
						 key=base),file.path(tmp,'parscale.RDS'))
kept <- funReadCachedParscale(file.path(tmp,'parscale.RDS'),parNames,currentKey=base)
check('a matching key lets the parscales through',all(kept$parscale==c(1,2,3)))
dropped <- funReadCachedParscale(file.path(tmp,'parscale.RDS'),parNames,
																 currentKey=buildKey(baseNegLL=999))
check('a mismatched key drops them',all(is.na(dropped$parscale)))

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('the determination key does not reject what it should\n')
}
