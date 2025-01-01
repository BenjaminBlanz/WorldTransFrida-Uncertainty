library('DescTools',quietly = T)
source('initialise.R')
source('config.R')
source('runInitialiseData.R')
source('setupTMPFS.R')
cat(sprintf('processing the results in\n %s\n',
						file.path(location.output,'detectedParmSpace')))
calDat <- readRDS(file.path(location.output,'calDat.RDS'))$calDat
resSigma <- readRDS(file.path(location.output,'sigma-indepParms.RDS'))

# run files ####
# for input
location.runFiles <- file.path(location.output,'detectedParmSpace')
runFilesList <- list.files(location.runFiles,pattern = 'workUnit-[0-9]+\\.RDS')

if(length(runFilesList)==0){
	stop('no run files to process\nHave you run runMLEandParmSpace?\n')
}

# read ####
cat('reading sampleParms...')
sampleParms <- readRDS(file.path(location.output,'sampleParmsParscaleRanged.RDS'))
cat('done\nreading samplePoints...')
samplePoints <- as.data.frame(readRDS(file.path(location.output,'samplePoints.RDS')))
if(nrow(samplePoints)!=numSample){
	stop('Number of sample points and numSample do not match\n')
}
cat('done\n')

# plot vars ####
defRun <- runFridaDefaultParms()
yearsToPlot <- rownames(defRun)
varsToPlot.lst <- list()
workUnitBoundaries <- seq(1,nrow(calDat)+1,5)
for(i in 1:(length(workUnitBoundaries)-1)){
	varsToPlot.lst[[i]] <- colnames(calDat)[workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)]
}

for(plotWeightType in plotWeightTypes){
	# weighting ####
	cat('')
	if(plotWeightType=='likelihood'){
		samplePoints$plotWeight <- exp(logLike)
	} else if(plotWeightType == 'logCutoff'){
		# somewhat wrong likelihood weighting
		# log like ####
		cat('reading log likelihoods...\n')
		logLike <- rep(NA,numSample)
		completeRunsSoFar <- 0
		for(f.i in 1:length(runFilesList)){
			cat(sprintf('\r chunk %i of %i',f.i,length(runFilesList)))
			parOutput <- readRDS(file.path(location.output,'detectedParmSpace',paste0('workUnit-',f.i,'.RDS')))
			for(l in 1:length(parOutput)){
				logLike[parOutput[[l]]$parmsIndex] <- parOutput[[l]]$logLike
				if(!is.na(parOutput[[l]]$runDat[[1]][length(parOutput[[l]]$runDat[[1]])])){
					completeRunsSoFar <- completeRunsSoFar + 1
				}
			}
			rm(parOutput)
		}
		cat(sprintf('\r read %i files, collected %i sample log likes, %i runs in data where complete\n',
								length(runFilesList),numSample,completeRunsSoFar))
		samplePoints$logLike <- logLike
		logLike.ecdf <- ecdf(logLike)
		samplePoints$plotWeight <- logLike.ecdf(logLike)
	} else if(plotWeightType == 'equaly'){
		# equal weighting
		samplePoints$plotWeight <- rep(1,length(logLike))
	} else if(plotWeightType == 'linearly'){
		samplePoints$plotWeight <- order(logLike)/nrow(samplePoints)
	} else {
		stop('unknown plotWeightType\n'	)
	}
	for(vars.i in 1:length(varsToPlot.lst)){
		varsToPlot <- varsToPlot.lst[[vars.i]]
		# collect time series ####
		cat(sprintf('reading %i vars...\n',length(varsToPlot)))
		cat(paste0(varsToPlot,collapse='\n'))
		cat('\n')
		# dimensions time in rows, run IDs in columns,  variables to read
		runsData <- array(NA,dim=c(length(yearsToPlot),nrow(samplePoints),length(varsToPlot)))
		dimnames(runsData) <- list(yearsToPlot,1:nrow(samplePoints),varsToPlot)
		for(f.i in 1:length(runFilesList)){
			cat(sprintf('\r chunk %i of %i',f.i,length(runFilesList)))
			parOutput <- readRDS(file.path(location.output,'detectedParmSpace',paste0('workUnit-',f.i,'.RDS')))
			for(l in 1:length(parOutput)){
				runsData[,parOutput[[l]]$parmsIndex,] <- unlist(parOutput[[l]]$runDat[,varsToPlot])
			}
			rm(parOutput)
		}
		cat('done\n')
		
		# CI plots ####
		location.plots.ci <- file.path(location.output,location.plots,plotWeightType,'CI-plots')
		dir.create(location.plots.ci,F,T)
		for(var.i in 1:length(varsToPlot)){
			cat(sprintf('Plotting %s\n',varsToPlot[var.i]))
			## calculate the  CIs ####
			cat('  calculating CI bourders...')
			CIsToPlot <- sort(CIsToPlot)
			ciBoundQs <- unique(c(rev((1-CIsToPlot)/2),1-(1-CIsToPlot)/2))
			ciBoundQs.lty <- c(rev(CIsToPlot.lty),CIsToPlot.lty[-1])
			ciBoundQs.lwd <- c(rev(CIsToPlot.lwd),CIsToPlot.lwd[-1])
			ciBoundQs.lcol <- c(rev(CIsToPlot.lcol),CIsToPlot.lcol[-1])
			medianQIdx <- length(CIsToPlot)+1
			ciBounds <- array(NA,dim=c(length(yearsToPlot),length(ciBoundQs)))
			colnames(ciBounds) <- ciBoundQs
			rownames(ciBounds) <- yearsToPlot																		
			for(year.i in 1:length(yearsToPlot)){
				ciBounds[year.i,] <- Quantile(runsData[year.i,,var.i],
																			weights = samplePoints$plotWeight,
																			probs = ciBoundQs,
																			na.rm = T)
			}
			cat('done\n')
			## draw ####
			cat('  drawing...')
			png(file.path(location.plots.ci,paste(varsToPlot[var.i],plotWeightType,'weighted',sep='-')),
					width = plotWidth,height = plotHeight,units = plotUnit,res = plotRes)
			varsForExport.cleanNames.orig <- cleanNames(varsForExport.cleanNames.orig)
			varName <- varsForExport.fridaNames.orig[which(varsForExport.cleanNames.orig==varsToPlot[var.i])]
			plot(yearsToPlot,ciBounds[,medianQIdx],
					 ylim=range(ciBounds[,c(2,length(ciBoundQs)-1)]),
					 type='n',
					 xlab='year',
					 ylab=varName,
					 xaxt='n',
					 main=varName)
			mtext(paste('Samples',plotWeightType,'weighted'),3,0.1)
			axis(1,at=seq(as.numeric(yearsToPlot[1]),as.numeric(yearsToPlot[length(yearsToPlot)]),10))
			legend('topleft',legend=c('median',paste0(CIsToPlot[-1]*100,'% CI'),'Data'),
						 lty=c(CIsToPlot.lty,NA),lwd=c(CIsToPlot.lwd,NA),
						 pch=c(rep(NA,length(CIsToPlot)),20), col = c(CIsToPlot.lcol,calDat.col))
			for(ci.i in length(CIsToPlot.col):1){
				if(CIsToPlot[ci.i]==0){
					#skip
				} else {
					if((length(ciBoundQs)%%2)==0){
						idxOfLowCiBounds1 <- length(ciBoundQs)/2
						ciBound.low <- ciBounds[,idxOfLowCiBounds1-ci.i+1]
						ciBound.high <- ciBounds[,idxOfLowCiBounds1+ci.i]
					} else {
						idxOfLowCiBounds1 <- length(ciBoundQs)/2
						ciBound.low <- ciBounds[,idxOfLowCiBounds1-ci.i+1.5]
						ciBound.high <- ciBounds[,idxOfLowCiBounds1+ci.i-.5]
					}
					polygon(c(yearsToPlot,rev(yearsToPlot)),c(ciBound.low,rev(ciBound.high)),
									col=CIsToPlot.col[ci.i],lty = 0)
				}
			}
			for(q.i in 1:length(ciBoundQs)){
				lines(yearsToPlot,ciBounds[,q.i],
							lty=ciBoundQs.lty[q.i],
							lwd=ciBoundQs.lwd[q.i],
							col=ciBoundQs.lcol[q.i])
			}
			points(rownames(calDat),calDat[[varsToPlot[var.i]]],
						 col='red',pch=20)
			dev.off()
			cat('done\n')
		}
		
		rm(runsData)
		gc()
		# for(var.i in 1:length(varsToPlot)){
			# # fan plots ####
			# location.plots.lines <- file.path(location.output,location.plots,'CI-plots')
			# dir.create(location.plots.lines,F,T)
			# png(file.path(location.plots.lines,paste(varsToPlot[var.i],plotWeightType,'weighted',sep='-')),
			# 		width = plotWidth,height = plotHeight,units = plotUnit,res = plotRes)
			# varName <- varsForExport.fridaNames.orig[which(varsForExport.cleanNames.orig==varsToPlot[var.i])]
			# plot(yearsToPlot,ciBounds[,medianQIdx],
			# 		 ylim=range(ciBounds[,c(2,length(ciBoundQs)-1)]),
			# 		 type='n',
			# 		 xlab='year',
			# 		 ylab=varName,
			# 		 xaxt='n',
			# 		 main=varName)
			# mtext(paste('Samples',plotWeightType,'weighted'),3,0.1)
			# axis(1,at=seq(as.numeric(yearsToPlot[1]),as.numeric(yearsToPlot[length(yearsToPlot)]),10))
			# legend('topleft',legend=c('median','runs','Data'),
			# 			 lty=c(CIsToPlot.lty[1],'solid',NA),lwd=c(CIsToPlot.lwd[1],1,NA),
			# 			 pch=c(NA,NA,20), col = c(CIsToPlot.col[1],CIsToPlot.col[2],calDat.col))
			# for(l.i in 1:100){#1:nrow(samplePoints)){
			# 	lines(yearsToPlot,runsData[,l.i,var.i],col=adjustcolor(1,alpha=samplePoints$plotWeight))
			# }
			# points(rownames(calDat),calDat[[varsToPlot[var.i]]],
			# 			 col='red',pch=20)
			# dev.off()
		# }
	}
}
	


	
	
	
