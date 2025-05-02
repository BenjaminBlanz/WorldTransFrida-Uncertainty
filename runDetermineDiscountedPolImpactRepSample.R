source('config-RunDetermineRepresentiveSamples.R')

# baseline run ####
# send off the baseline run
expIDpreString <- 'discountedPolicyImpact'
if(is.na(preexsistingBaselineFolder)){
	baselineExpID <- paste0(expIDpreString,'_Baseline-S',numSample,'-',specFileForScenaro,'-ClimateFeedback_On-ClimateSTAOverride_Off')
}
statusFile <- file.path('workOutput',baselineExpID,'status')
if(!file.exists(file.path('workOutput',baselineExpID))){
	system(paste('./submit_UncertaintyAnalysisLevante.sh',
							 '-n',numSample,
							 '--pol',specFileForScenaro,
							 '-s',paste0(expIDpreString,'_Baseline')))
	setwd(aggDamWD)
	stop('Baseline run has been submitted to SLURM, please restart this script once the baseline run has completed\n')
}
if(!file.exists(statusFile)||
	 !readChar(statusFile,file.info(statusFile)$size-1)=='completed'){
	setwd(aggDamWD)
	stop('Baseline run has not completed yet. Run this script again once it has completed.\n')
} else {
	cat('Baseline completed continuing\n')
}

# policy runs ####
# send off the policy runs
poldRunsExpIDpreSring <- ''
polRuns <- data.frame(policyName=specFilesForPols,expID=NA,status=NA)
for(pol.i in 1:length(specFilesForPols)){
	expID <- paste0(expIDpreString,'-S',numSample,'-',polRuns$policyName[pol.i],'-ClimateFeedback_On-ClimateSTAOverride_Off')
	polRuns$expID[pol.i] <- expID
	statusFile <- file.path('workOutput',expID,'status')
	if(!file.exists(statusFile)){
		polRuns$status[pol.i] <- 'not present'
	} else if (readChar(statusFile,file.info(statusFile)$size-1)=='completed'){
		polRuns$status[pol.i] <- 'completed'
	} else {
		polRuns$status[pol.i] <- 'presumed running'
	}
}
if(sum(polRuns$status=='not present')>0){
	for(pol.i in which(polRuns$status=='not present')){
		expID <- polRuns$expID[pol.i]
		system(paste('./submit_UncertaintyAnalysisLevante.sh',
								 '-n',numSample,
								 '--pol',polRuns$policyName[pol.i],
								 '--cfb','ClimateFeedback_On.csv',
								 '-s',expIDpreString,
								 '--cid',baselineExpID,
								 '--cpps','true',
								 '--cpsp','true'))
	}
	setwd(aggDamWD)
	stop('Missing pol runs have been submitted to SLURM. Please restart this script once the pol runs have completed\n')
} else if(sum(polRuns$status=='completed')<nrow(polRuns)){
	setwd(aggDamWD)
	stop('pol runs have not completed yet. Please restart this script when all of them have completed\n')
} else {
	cat('pol runs completed continuing\n')
}

# collect time series ####
cat('reading vars from run data...')
baselineValues <- list()
for(v.i in 1:length(subSample.TargetVars)){
	baselineValues[[subSample.TargetVars[v.i]]] <- 
		readRDS(file.path(baselineExpID,'detectedParmSpace','PerVarFiles-RDS',subSample.TargetVars[v.i]))
}
policyValues <- list()
for(pol.i in 1:length(specFilesForPols)){
	policyValues[[specFilesForPols[pol.i]]] <- list()
	for(v.i in 1:length(subSample.TargetVars)){
		policyValues[[specFilesForPols[pol.i]]][[subSample.TargetVars[v.i]]] <- 
			readRDS(file.path(polRuns$expID[pol.i],'detectedParmSpace','PerVarFiles-RDS',subSample.TargetVars[v.i]))
	}	
}
cat('done\n')

# calculate policy impacts ####
policyDiffs <- list()
policyDiscountedDiffs <- list()
repSampleIdx <-list()
deflator <- 1/(1+discountFactor)^(c(rep(0,2025-1980),1:(ncol(baselineValues[[1]])-(2025-1980))))
for(pol.i in 1:length(specFilesForPols)){
	policyDiffs[[specFilesForPols[pol.i]]] <- list()
	policyDiscountedDiffs[[specFilesForPols[pol.i]]] <- list()
	repSampleIdx[[specFilesForPols[pol.i]]] <- list()
	for(v.i in 1:length(subSample.TargetVars)){
		diff <- policyValues[[specFilesForPols[pol.i]]][[subSample.TargetVars[v.i]]] - 
			baselineValues[[subSample.TargetVars[v.i]]]
		policyDiffs[[specFilesForPols[pol.i]]][[subSample.TargetVars[v.i]]] <- diff
		policyDiscountedDiffs[[specFilesForPols[pol.i]]][[subSample.TargetVars[v.i]]] <-
			diff %*% deflator
		repSampleIdx[[specFilesForPols[pol.i]]][[subSample.TargetVars[v.i]]] <- c()
		quantiles <- quantile(policyDiscountedDiffs[[specFilesForPols[pol.i]]][[subSample.TargetVars[v.i]]],
													subSample.Ps,na.rm=T)
		for(p.i in 1:length(subSample.Ps)){
			repSampleIdx[[specFilesForPols[pol.i]]][[subSample.TargetVars[v.i]]][p.i] <-
				which.min((policyDiscountedDiffs[[specFilesForPols[pol.i]]][[subSample.TargetVars[v.i]]]-quantiles[p.i])^2)
		}
	}	
}

samplePoints <- as.data.frame(readRDS(file.path(baselineExpID,'samplePoints.RDS')))
repSample <- samplePoints[unique(unlist(repSampleIdx)),]
colnames(repSample) <- gsub('\\[1\\]','',colnames(repSample))
colnames(repSample) <- gsub('\\[1,','[*,',colnames(repSample))

cat('writing out to subSampleParameterValues.csv ...')
write.table(repSample,file.path(location.output,'subSampleParameterValuesDPIS.csv'),
						append = F,sep = ',',row.names = F)
cat('done\n')
