# testParameterFileColumns.R ####
#
# Parameter files, i.e. frida_info, the external ranges and the integer parms,
# name their parameters in a Variable or a Parameter column. Stella's own
# Parameter Info.csv has both, Variable for the names and Parameter as a flag,
# and keeps them as they are.
#
# Run from the repository root:
#   Rscript developmentTools/testParameterFileColumns.R

source('funRunFRIDA.R')

tmp <- tempfile('parameterFileColumnsTest')
dir.create(tmp)
on.exit(unlink(tmp,recursive=TRUE))

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),''))); fail <<- fail+1}
}
csvFile <- function(name,lines){
	f <- file.path(tmp,name)
	writeLines(lines,f)
	return(f)
}

cat('reading a parameter file\n')
p <- funReadParameterFile(csvFile('parameter.csv',
	c('Parameter,Value,Min,Max,Min Multi,Max Multi','Finance.loan maturation time[1],10,9,13,0.9,1.3')))
check('a Parameter column becomes Variable',identical(p$Variable,'Finance.loan maturation time[1]'),
			paste(colnames(p),collapse=', '))
check('the other columns are read as read.csv reads them',all(c('Value','Min','Max','Min.Multi')%in%colnames(p)))
v <- funReadParameterFile(csvFile('stella.csv',c('Variable,Value,Min,Max,Parameter,Policy','x[1],1,0,2,1,')))
check('with a Variable column, a Parameter column stays a column of its own',
			identical(v$Variable,'x[1]')&&identical(v$Parameter,1L))

cat('sampling the parameters of a user frida_info with a Parameter column\n')
location.frida.info <- tmp
location.output <- tmp
name.frida_info <- 'parameter.csv'
frida_info_type <- 'user'
redoAllCalc <- FALSE
name.frida_parameter_exclusion_list <- 'none.csv'
invisible(capture.output(sampleParms <- prepareSampleParms()))
check('the parameter is sampled',identical(sampleParms$Variable,'Finance.loan maturation time[1]')&&
				sampleParms$Min==9&&sampleParms$Max==13)

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('parameter files are not read as they should be\n')
}
