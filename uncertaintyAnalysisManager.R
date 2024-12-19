#
# Manager script for the uncertainty analysis of FRIDA
# 

require(SobolSequence)
require(tictoc)
require(parallel)
require(scales)

# config ####
cat('Config...')
numWorkers <- 10
numSample <- 5e5

# How large the chunks of work are, smaller means more frequent pauses to write out
# itermediate results (and update the diagnostic output).
chunkSizePerWorker <- 100

# by default sobol sequence covers the entire range between min and max with 
# equal density.
# However we might want to ensure that there are similar number of points above and 
# below the Value in our baseline calibration, our prior.
restretchSamplePoints <- T

plotWhileRunning <- T

location.frida <- './FRIDAforUncertaintyAnalysis'
location.stella<- './Stella_Simulator_Linux'
name.fridaInputFile <- 'uncertainty_analysis_paramter_values.csv'
name.fridaOutputFile <- 'uncertainty_analysis_exported_variables.csv'

location.output <- file.path('workOutput',paste0('NumSample-',numSample,'-chunkSizePerWorker-',chunkSizePerWorker))
dir.create(location.output,showWarnings=F,recursive=T)
cat('done\n')


# specify sampling parameters ####
cat('Specify sampling parameter...')
# read in the parameters in frida that have ranges defined
frida_info <- read.csv("frida_info.csv")
columnsThatAreFlags <- c(2,3,4,5,6,7,8,9,10)
# select the parameters to be sampled
sampleParms <- frida_info[rowSums(frida_info[,columnsThatAreFlags])>0 & frida_info$No.Sensi==0 & frida_info$Policy==0,-columnsThatAreFlags]
cat('done\n')

# generate sobol sequence ####
cat('Generate sobol sequence...')
# sobolSequence.points generates points on the unit interval for each var
# transformed, so vars are in rows samples in cols, makes the next steps easier
samplePoints <- sobolSequence.points(nrow(sampleParms),31,numSample)
if(sum(duplicated(samplePoints))<numSample){
	stop('Not enough unique sample points. Check the sobol generation\n')
}
samplePoints <- t(samplePoints)

if(!restretchSamplePoints){
	# Substract the min and multiply by max-min to strecth the unit interval to the 
	# actual sampling range.
	samplePointsStretched <- samplePoints*(sampleParms$Max-sampleParms$Min) + sampleParms$Min
	# plot(samplePointsStretched[1,],samplePointsStretched[2,])
	# abline(v=sampleParms$Value[1],h=sampleParms$Value[2],col='red')
	samplePoints <- samplePointsStretched
} else {
	# stretch the sample points to be left and right of the mean centre value of the 
	# description file
	lowIdc <- samplePoints<0.5
	highIdc <- samplePoints>=0.5
	samplePointsLow <- samplePointsHigh <- samplePoints
	samplePointsLow[highIdc] <- NA
	samplePointsLow <- samplePointsLow*2*(sampleParms$Value-sampleParms$Min) + sampleParms$Min
	samplePointsHigh[lowIdc] <- NA
	samplePointsHigh <- (samplePoints-0.5)*2*(sampleParms$Max-sampleParms$Value) + sampleParms$Value
	samplePointsReStretched <- samplePoints
	samplePointsReStretched[lowIdc] <- samplePointsLow[lowIdc]
	samplePointsReStretched[highIdc] <- samplePointsHigh[highIdc]
	# plot(samplePointsReStretched[1,],samplePointsReStretched[2,])
	# abline(v=sampleParms$Value[1],h=sampleParms$Value[2],col='red')
	samplePoints <- samplePointsReStretched
}
cat('done\n')

# run FRIDA with the samples ####
## fun write frida input ####
writeFRIDAInput <- function(variables,values){
	parmValues <- data.frame(Variable=variables,Value=values)
	write.table(parmValues,file = file.path(location.frida,'Data',name.fridaInputFile),
							row.names = F,col.names = F,sep=',')
}
## cluster ####
### fun runFridaPar ####
runFridaPar <- function(i){
	writeFRIDAInput(sampleParms$Variable,samplePoints[,i])
	system(paste(file.path(location.stella,'stella_simulator'),'-i','-x','-q',file.path(location.frida,'FRIDA.stmx')),
				 ignore.stdout = T,ignore.stderr = T,wait = T)
	return(read.csv(file.path(location.frida,'Data',name.fridaOutputFile)))
}
### cluster setup ####
cat('cluster setup...')
baseWD <- getwd()
workDirBasename <- 'workDir_'
# start cluster
cl <- makeForkCluster(numWorkers,renice=15)
workers <- 1:length(cl)
# make working directories
gobble <- clusterApply(cl,workers,function(i){dir.create(file.path('workerDirs',paste0(workDirBasename,i)),showWarnings = F,recursive = T)}) 
gobble <- clusterApply(cl,workers,function(i){setwd(file.path('workerDirs',paste0(workDirBasename,i)))})
#clusterEvalQ(cl,getwd())
# copy over the model and simulator
gobble <- clusterApply(cl,workers,function(i){system(paste('cp -r',file.path(baseWD,location.frida),getwd()))})
gobble <- clusterApply(cl,workers,function(i){system(paste('cp -r',file.path(baseWD,location.stella),getwd()))})
cat('done\n')

### cluster run ####
cat('cluster run...')
if(plotWhileRunning){
	plot(0,0,type='n',ylim=c(0,1e6),xlim=c(1980,2130),
			 xlab='year',ylab='Real GDP in 2021 bn intl$',
			 xaxt='n')
	axis(1,at=seq(1980,2130,10))
}
workUnitBoundaries <- seq(1,numSample,chunkSizePerWorker*length(cl))
# add one to the last work unit boundary, as during running we always deduct one from the next boundary
workUnitBoundaries[length(workUnitBoundaries)] <- numSample+1
cat(sprintf('  Run of %i runs split up into %i work units.\n',numSample,length(workUnitBoundaries)-1))
chunkTimes <- c()
for(i in 1:(length(workUnitBoundaries)-1)){
	if(file.exists(file.path(location.output,paste0('workUnit-',i,'.RDS')))){
		cat(sprintf('\r   Skipping unit %i already exists',i))
	} else {
		cat(sprintf('\r   Running unit %i: samples %i to %i', i, workUnitBoundaries[i],workUnitBoundaries[i+1]-1))
		if(length(chunkTimes>1)){
			cat(sprintf(', average duration per unit so far %f, expect completion in %f',
									mean(chunkTimes,na.rm=T),mean(chunkTimes,na.rm=T)*(length(workUnitBoundaries)-1-i)))
		}
		tic()
		workUnit <- workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)
		runDat <- parLapply(cl=cl,X=workUnit,fun = runFridaPar)
		timing <- toc(quiet=T)
		chunkTimes[i] <- timing$toc-timing$tic
		saveRDS(runDat,file.path(location.output,paste0('workUnit-',i,'.RDS')))
		if(plotWhileRunning){
			# readline(prompt="Press [enter] to continue")
			for(l in 1:length(runDat)){
				lines(runDat[[l]]$Year,runDat[[l]]$GDP.Real.GDP.in.2021c..1.,col=alpha(i,0.25))
			}
		}
		rm(runDat)
	}
}
cat(sprintf('\r    runs completed average chunk time %f, over all run time %f\n',mean(chunkTimes,na.rm=T),sum(chunkTimes,na.rm=T)))
cat('done')

### cluster cleanup ####
cat('cluster cleanup...')
# stop cluster
stopCluster(cl)
# clean up working directories
system('rm -r workerDirs')
cat('done\n')


