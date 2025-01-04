source('initialise.R')
source('config.R')
source('runInitialiseData.R')
source('setupTMPFS.R')
cat(sprintf('processing the results in\n %s\n',
						file.path(location.output,'detectedParmSpace')))
calDat <- readRDS(file.path(location.output,'calDat.RDS'))$calDat
resSigma <- readRDS(file.path(location.output,'sigma-indepParms.RDS'))
colnames(resSigma) <- rownames(resSigma) <- colnames(calDat)

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
yearsToPlot.lst <- list()
for (y.i in 1:length(yearsToPlot.names)){
	if(yearsToPlot.names[y.i] == 'allYears'){
		yearsToPlot.lst[[y.i]] <- rownames(defRun)
	} else {
		startEndYears <- as.numeric(StrSplit(yearsToPlot.names[y.i],split = '-'))
		yearsToPlot.lst[[y.i]] <- as.character(startEndYears[1]:startEndYears[2])
	}
}
varsToPlot.lst <- list()
workUnitBoundaries <- seq(1,ncol(calDat)+1,5)
workUnitBoundaries <- c(workUnitBoundaries,ncol(calDat)+1)
for(i in 1:(length(workUnitBoundaries)-1)){
	varsToPlot.lst[[i]] <- colnames(calDat)[workUnitBoundaries[i]:(workUnitBoundaries[i+1]-1)]
}
# varsToPlot.lst <- rev(varsToPlot.lst)


for(plotWeightType in plotWeightTypes){
	# weighting ####
	cat(sprintf('Plotting %s weighted\n',plotWeightType))
	if(!exists('logLike')&&plotWeightType %in% c('likelihood','logCutoff','linearly')){
		# log like ####
		cat(' reading log likelihoods...\n')
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
	}
	
	if(plotWeightType=='likelihood'){
		samplePoints$plotWeight <- exp(logLike)
	} else if(plotWeightType == 'logCutoff'){
		# somewhat wrong likelihood weighting
		samplePoints$plotWeight <- logLike.ecdf(logLike)
	} else if(plotWeightType == 'equaly'){
		# equal weighting
		samplePoints$plotWeight <- rep(1,nrow(samplePoints))
	} else if(plotWeightType == 'linearly'){
		samplePoints$plotWeight <- order(logLike)/nrow(samplePoints)
	} else {
		stop('unknown plotWeightType\n'	)
	}
	
	location.progress <- file.path(location.output,location.plots,'CI-plots',
																 paste0(plotWeightType,'Weighted'),'plotProgress')
	dir.create(location.progress,F,T)
	for(vars.i in 1:length(varsToPlot.lst)){
		varsToPlot <- varsToPlot.lst[[vars.i]]
		if(sum(file.exists(file.path(location.progress,varsToPlot)))==length(varsToPlot)){
			cat(sprintf(' vars %i to %i have been completed already (according to plotProgress files)',
									workUnitBoundaries[vars.i],workUnitBoundaries[vars.i+1]-1))
		} else {
			# collect time series ####
			cat(sprintf('Reading %i vars: %i to %i of %i total vars...\n   ',
									length(varsToPlot),workUnitBoundaries[vars.i],workUnitBoundaries[vars.i+1]-1,ncol(calDat)))
			cat(paste0(varsToPlot,collapse='\n   '))
			cat('\n')
			# dimensions time in rows, run IDs in columns,  variables to read
			# read all the years, selectively plot later
			runsData <- array(NA,dim=c(nrow(defRun),nrow(samplePoints),length(varsToPlot)))
			dimnames(runsData) <- list(rownames(defRun),1:nrow(samplePoints),varsToPlot)
			for(f.i in 1:length(runFilesList)){
				cat(sprintf('\r  reading chunk %i of %i',f.i,length(runFilesList)))
				parOutput <- readRDS(file.path(location.output,'detectedParmSpace',paste0('workUnit-',f.i,'.RDS')))
				for(l in 1:length(parOutput)){
					runsData[,parOutput[[l]]$parmsIndex,] <- unlist(parOutput[[l]]$runDat[,varsToPlot])
				}
				rm(parOutput)
			}
			cat('\r  reading done                                                            \n')
			
			# CI plots ####
			for(var.i in 1:length(varsToPlot)){
				cat(sprintf('  %i-%i: %s...',vars.i,var.i,varsToPlot[var.i]))
				## calculate the  CIs ####
				cat('CI borders...')
				CIsToPlot <- sort(CIsToPlot)
				ciBoundQs <- unique(c(rev((1-CIsToPlot)/2),1-(1-CIsToPlot)/2))
				ciBoundQs.lty <- c(rev(CIsToPlot.lty),CIsToPlot.lty[-1])
				ciBoundQs.lwd <- c(rev(CIsToPlot.lwd),CIsToPlot.lwd[-1])
				ciBoundQs.lcol <- c(rev(CIsToPlot.lcol),CIsToPlot.lcol[-1])
				medianQIdx <- length(CIsToPlot)+1
				ciBounds <- array(NA,dim=c(nrow(defRun),length(ciBoundQs)))
				means <- rep(NA,nrow(defRun))
				names(means) <- rownames(ciBounds) <- rownames(defRun)
				colnames(ciBounds) <- ciBoundQs
				for(year.i in 1:nrow(defRun)){
					means[year.i] <- weighted.mean(runsData[year.i,,varsToPlot[var.i]],
																				 w = samplePoints$plotWeight,
																				 na.rm=T)
					ciBounds[year.i,] <- Quantile(runsData[year.i,,var.i],
																				weights = samplePoints$plotWeight,
																				probs = ciBoundQs,
																				na.rm = T)
				}
				cat('\n    ')
				means.store <- means
				ciBounds.store <- ciBounds
				for(uncertaintyType in uncertaintiesToPlot){
					cat(sprintf('%s...',uncertaintyType))
					for(year.i in 1:nrow(defRun)){
						if(uncertaintyType=='all uncertainty'||uncertaintyType=='fit uncertainty'){
							means[year.i] <- means.store[year.i]
							ciBounds[year.i,] <- ciBounds.store[year.i,]
						} else if (uncertaintyType=='noise uncertainty'){
							means[year.i] <- defRun[year.i,varsToPlot[var.i]]
							ciBounds[year.i,] <- defRun[year.i,varsToPlot[var.i]]
						} else {
							stop('unkown uncertaintyType\n')
						}
						if(uncertaintyType=='noise uncertainty'||uncertaintyType=='all uncertainty'){
							if(treatVarsAsIndep){
								ciBounds[year.i,] <- ciBounds[year.i,] + qnorm(ciBoundQs,mean = 0, sd=sqrt(resSigma[varsToPlot[var.i],varsToPlot[var.i]]))
							} else {
								stop('Noise uncertainty not implemented yet for dependent vars\n')
							}
						}
					}
					## save data ####
					plotData <- list()
					plotData$variable <- varsToPlot[var.i]
					plotData$years <- rownames(defRun)
					plotData$ciBoundQs <- ciBoundQs
					plotData$ciBounds <- ciBounds
					plotData$uncertaintyType <- uncertaintyType
					plotData$means <- means
					plotData$defaultRun <- defRun
					dir.create(file.path(location.output,location.plots,'CI-plots',
															 paste0(plotWeightType,'Weighted'),'plotData'),F,T)
					saveRDS(plotData,file.path(location.output,location.plots,'CI-plots',
																		 paste0(plotWeightType,'Weighted'),'plotData',
																		 paste0(paste(varsToPlot[var.i],plotWeightType,'weighted',sep='-'),'.RDS')))
					## draw ####
					cat('drawing...')
					for(years.i in 1:length(yearsToPlot.lst)){
						for(alsoPlotMean in alsoPlotMean.vals){
							for(alsoPlotDefaultRun in alsoPlotDefaultRun.vals){
								yearsToPlot <- yearsToPlot.lst[[years.i]]
								location.plots.ci <- file.path(location.output,location.plots,'CI-plots',
																							 paste0(plotWeightType,'Weighted'),
																							 yearsToPlot.names[years.i],
																							 paste(uncertaintyType,
																							 			ifelse(alsoPlotMean,'withMean','withoutMean'),
																							 			ifelse(alsoPlotDefaultRun,'withDefaultRun','withoutDefaultRun'),sep='-'))
								dir.create(location.plots.ci,F,T)
								png(file.path(location.plots.ci,paste0(paste(varsToPlot[var.i],plotWeightType,'weighted',sep='-'),'.png')),
										width = plotWidth,height = plotHeight,units = plotUnit,res = plotRes)
								varsForExport.cleanNames.orig <- cleanNames(varsForExport.cleanNames.orig)
								varName <- varsForExport.fridaNames.orig[which(varsForExport.cleanNames.orig==varsToPlot[var.i])]
								layout(matrix(c(2,1),nrow=2),heights = c(0.9,0.1))
								par(mar=c(0,0,0,0))
								plot(0,0,type='n',axes=F)
								legend.text=c(
									if(alsoPlotMean&&uncertaintyType!='noise uncertainty'){'mean'},
									if(alsoPlotDefaultRun){'frida default'},
									'median',
									paste0(CIsToPlot[-1]*100,'% CI'),
									'Data')
								legend.lty=c(
									if(alsoPlotMean&&uncertaintyType!='noise uncertainty'){mean.lty},
									if(alsoPlotDefaultRun){def.lty},
									CIsToPlot.lty,
									NA)
								legend.lwd=c(
									if(alsoPlotMean&&uncertaintyType!='noise uncertainty'){mean.lwd},
									if(alsoPlotDefaultRun){def.lwd},
									CIsToPlot.lwd,
									NA)
								legend.pch=c(rep(NA,length(CIsToPlot)+sum(c(alsoPlotMean&uncertaintyType!='noise uncertainty',alsoPlotDefaultRun))),20)
								legend.col = c(
									if(alsoPlotMean&&uncertaintyType!='noise uncertainty'){mean.col},
									if(alsoPlotDefaultRun){def.col},
									CIsToPlot.lcol,
									calDat.col)
								legend('bottom',legend.text,lty=legend.lty,lwd=legend.lwd,pch=legend.pch,col=legend.col,
											 horiz=T,xpd=T)
								par(mar=c(3.1,4.1,4.1,2.1))
								plot(yearsToPlot,ciBounds[yearsToPlot,medianQIdx],
										 ylim=range(c(ciBounds[yearsToPlot,c(2,length(ciBoundQs)-1)],calDat[yearsToPlot,varsToPlot[var.i]]),na.rm=T),
										 xlim=range(as.numeric(yearsToPlot)),
										 xaxs='i',
										 type='n',
										 xlab='',
										 ylab=varName,
										 xaxt='n',
										 main=varName)
								mtext(paste('Samples',plotWeightType,'weighted. Ranges show ',uncertaintyType,'.'),3,0.5)
								xax <- axis(1,at=seq(as.numeric(yearsToPlot[1]),as.numeric(yearsToPlot[length(yearsToPlot)]),10))
								mtext('year',1,3)
								grid(nx=length(xax)-1,ny=NA)
								abline(h=axTicks(2),lty='dotted',col='gray')
								for(ci.i in length(CIsToPlot.col):1){
									if(CIsToPlot[ci.i]==0){
										#skip
									} else {
										if((length(ciBoundQs)%%2)==0){
											idxOfLowCiBounds1 <- length(ciBoundQs)/2
											ciBound.low <- ciBounds[,idxOfLowCiBounds1-ci.i+1]
											ciBound.high <- ciBounds[yearsToPlot,idxOfLowCiBounds1+ci.i]
										} else {
											idxOfLowCiBounds1 <- length(ciBoundQs)/2
											ciBound.low <- ciBounds[yearsToPlot,idxOfLowCiBounds1-ci.i+1.5]
											ciBound.high <- ciBounds[yearsToPlot,idxOfLowCiBounds1+ci.i-.5]
										}
										polygon(c(yearsToPlot,rev(yearsToPlot)),c(ciBound.low,rev(ciBound.high)),
														col=CIsToPlot.col[ci.i],lty = 0)
									}
								}
								for(q.i in 1:length(ciBoundQs)){
									lines(yearsToPlot,ciBounds[yearsToPlot,q.i],
												lty=ciBoundQs.lty[q.i],
												lwd=ciBoundQs.lwd[q.i],
												col=ciBoundQs.lcol[q.i])
								}
								if(alsoPlotMean){
									lines(yearsToPlot,means[yearsToPlot],lty=mean.lty,lwd=mean.lwd,col=mean.col)
								}
								if(alsoPlotDefaultRun){
									lines(yearsToPlot,defRun[yearsToPlot,varsToPlot[var.i]],lty=def.lty,lwd=def.lwd,col=def.col)
								}
								points(rownames(calDat),calDat[[varsToPlot[var.i]]],
											 col='red',pch=20)
								box()
								dev.off()
							}
						}
					}
					system(paste('touch',file.path(location.progress,varsToPlot[var.i])))
				}
				cat('done\n')
			}
			rm(runsData)
			gc()
		}
	}
}


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




