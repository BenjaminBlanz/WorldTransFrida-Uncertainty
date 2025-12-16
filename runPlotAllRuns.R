source('initialise.R')
source('config.R')
outputFolder <- file.path(location.output,'detectedParmSpace')
outputTypeFolders <- basename(list.dirs(outputFolder,recursive = F))
if(length(outputTypeFolders)==0){
	stop('PerVar output folders have not been created yet. Outdated output.\n')
}
if(length(grep('RDS',outputTypeFolders))>0){
	outputTypeFolder <- outputTypeFolders[grep('RDS',outputTypeFolders)]
} else {
	outputTypeFolder <- outputTypeFolders[1]
}
outputType <- strsplit(outputTypeFolder,'-')[[1]][2]
cat(sprintf('processing the results in\n %s\n',
						file.path(outputFolder,outputTypeFolder)))
calDat <- readRDS(file.path(location.output,'calDat.RDS'))$calDat
calDat.orig <- read.csv(file.path(location.frida,'Data','Calibration Data.csv'))
varsForExport.fridaNames.orig <- calDat.orig[,1]

resSigma <- readRDS(file.path(location.output,'sigma-indepParms.RDS'))
colnames(resSigma) <- rownames(resSigma) <- colnames(calDat)

# run files ####
allVarNames <- read.csv(file.path(location.frida.info,name.frida_extra_variables_to_export_list))$FRIDA.FQN
allVarNames <- allVarNames[nchar(allVarNames)>4]
allVarNames.orig <- c(varsForExport.fridaNames.orig,allVarNames)
allVarNames.orig <- gsub(' \\[(\\d)\\]','\\[\\1\\]',allVarNames.orig)
allVarNames.orig <- unique(allVarNames.orig)
allVarNames <- cleanNames(allVarNames.orig)

# read ####
cat('reading sampleParms...')
sampleParms <- readRDS(file.path(location.output,'sampleParmsParscaleRanged.RDS'))
cat('done\nreading samplePoints...')
# samplePoints <- as.data.frame(readRDS(file.path(location.output,'samplePoints.RDS')))
samplePoints <- data.frame(parValuesOmitted=rep(NA,numSample))
if(nrow(samplePoints)!=numSample){
	stop('Number of sample points and numSample do not match\n')
}
cat('done\n')

# plot vars ####
writeFRIDAExportSpec(varsForExport.fridaNames = allVarNames.orig,
										 location.frida)
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

dir.create(file.path(location.output,location.plots),recursive = T,showWarnings = F)
sink(file.path(location.output,location.plots,'pathToTheOutputUnderlyingTheseFigures.txt'))
cat(file.path(getwd(),location.output))
file.copy(file.path(location.output,'config.R'),file.path(location.output,location.plots,'configOfTheUnderlyingRuns.R'))
sink()

for(plotWeightType in plotWeightTypes){
	# weighting ####
	cat(sprintf('Plotting %s weighted\n',plotWeightType))
	if(plotWeightType %in% c('likelihood','logCutoff','linearly','logLikelihood','completeEqually')){
		# log like ####
		cat(' reading log likelihoods...\n')
		logLike <- readPerVarFile(file.path(outputFolder,outputTypeFolder,'logLike'),outputType)$logLike
		completeRunsSoFar <- sum(logLike > -.Machine$double.xmax+(1000*.Machine$double.eps))
		cat(sprintf('Collected %i sample log likes, %i runs in data where complete\n',
								numSample,completeRunsSoFar))
		samplePoints$logLike <- logLike
		logLike.ecdf <- ecdf(logLike)	
	}
	
	if(plotWeightType=='likelihood'){
		samplePoints$plotWeight <- exp(logLike)
	} else if(plotWeightType == 'logCutoff'){
		# somewhat wrong likelihood weighting
		samplePoints$plotWeight <- logLike.ecdf(logLike)
	} else if(plotWeightType == 'logLikelihood'){
		# somewhat wrong likelihood weighting
		rangeLogLike <- range(logLike,na.rm=T)
		samplePoints$plotWeight <- (logLike-rangeLogLike[1])/(rangeLogLike[2]-rangeLogLike[1])
	} else if(plotWeightType == 'equaly'){
		# equal weighting
		samplePoints$plotWeight <- rep(1,nrow(samplePoints))
	} else if(plotWeightType == 'completeEqually'){
		# equal weighting of completed runs
		samplePoints$plotWeight <- 0
		samplePoints$plotWeight[samplePoints$logLike > -.Machine$double.xmax+(1000*.Machine$double.eps)] <- 1
	}else if(plotWeightType == 'linearly'){
		samplePoints$plotWeight <- order(logLike)/nrow(samplePoints)
	} else {
		stop('unknown plotWeightType\n'	)
	}
	samplePoints$plotWeight[is.na(samplePoints$plotWeight)] <- 0
	for(varName.i in 1:length(allVarNames)){
		varName <- allVarNames[varName.i]
		varName.orig <- allVarNames.orig[which(allVarNames==varName)]
		if(length(varName.orig)>1){varName.orig<-varName.orig[1]}
		cat(sprintf('(%i of %i) Plotting %s...\n  ',
								varName.i,length(allVarNames),varName.orig))
		# assemble the list of all figures that will be created to check if they arleady
		# exit and we can skip
		figuresToBeCreated <- c()
		for(uncertaintyType in uncertaintiesToPlot){
			for(years.i in 1:length(yearsToPlot.lst)){
				for(alsoPlotMean in alsoPlotMean.vals){
					for(alsoPlotDefaultRun in alsoPlotDefaultRun.vals){
						for(alsoPlotRepSample in alsoPlotRepSample.vals){
							yearsToPlot <- yearsToPlot.lst[[years.i]]
							location.plots.ci <- file.path(location.output,location.plots,'CI-plots',
																						 paste0(plotWeightType,'Weighted'),
																						 yearsToPlot.names[years.i],
																						 paste(uncertaintyType,
																						 			ifelse(alsoPlotMean,'withMean','withoutMean'),
																						 			ifelse(alsoPlotDefaultRun,'withDefaultRun','withoutDefaultRun'),
																						 			ifelse(alsoPlotRepSample,'withRepSample','withoutRepSample'),sep='-'))
							figuresToBeCreated[length(figuresToBeCreated)+1] <- 
								file.path(location.plots.ci,paste0(paste(varName,plotWeightType,'weighted',sep='-'),'.png'))
						}
					}
				}
			}
		}
		if(sum(file.exists(figuresToBeCreated))==length(figuresToBeCreated)){
			cat('All figures already created, skipping.\n')
			next
		}
		cat('read...')
		if(!file.exists(file.path(outputFolder,outputTypeFolder,paste0(varName,'.',outputType)))){
			cat('missing\n')
			next
		}
		varData <- readPerVarFile(file.path(outputFolder,outputTypeFolder,varName),outputType)
		varData <- sort_by(varData,varData$id)
		varData <- t(varData[,-1])
		colnames(varData) <- gsub('(^X)([0-9]{4})','\\2',colnames(varData),perl = T)
		cat('CI borders...')
		CIsToPlot <- sort(CIsToPlot)
		ciBoundQs <- unique(c(rev((1-CIsToPlot)/2),1-(1-CIsToPlot)/2))
		ciBoundQs.lty <- c(rev(CIsToPlot.lty),CIsToPlot.lty[-1])
		ciBoundQs.lwd <- c(rev(CIsToPlot.lwd),CIsToPlot.lwd[-1])
		ciBoundQs.lcol <- c(rev(CIsToPlot.lcol),CIsToPlot.lcol[-1])
		medianQIdx <- which(ciBoundQs==0.5)
		ciBounds <- array(NA,dim=c(nrow(defRun),length(ciBoundQs)))
		means <- rep(NA,nrow(defRun))
		names(means) <- rownames(ciBounds) <- rownames(defRun)
		colnames(ciBounds) <- ciBoundQs
		for(year.i in 1:nrow(defRun)){
			means[year.i] <- weighted.mean(varData[year.i,],
																		 w = samplePoints$plotWeight,
																		 na.rm=T)
			ciBounds[year.i,] <- spatstat.univar::weighted.quantile(x=varData[year.i,],
																															w = samplePoints$plotWeight,
																															probs = ciBoundQs,
																															na.rm=T)
		}
		means.store <- means
		ciBounds.store <- ciBounds
		for(uncertaintyType in uncertaintiesToPlot){
			cat(sprintf('%s...',uncertaintyType))
			if(uncertaintyType %in% c('all uncertainty','noise uncertainty')
				 && ! varName %in% rownames(resSigma)){
				cat('skip (var not in Sigma)\n')
			} else {
				for(year.i in 1:nrow(defRun)){
					if(uncertaintyType=='all uncertainty'||uncertaintyType=='fit uncertainty'){
						means[year.i] <- means.store[year.i]
						ciBounds[year.i,] <- ciBounds.store[year.i,]
					} else if (uncertaintyType=='noise uncertainty'){
						means[year.i] <- defRun[year.i,varName]
						ciBounds[year.i,] <- defRun[year.i,varName]
					} else {
						stop('unkown uncertaintyType\n')
					}
					if(uncertaintyType=='noise uncertainty'||uncertaintyType=='all uncertainty'){
						if(treatVarsAsIndep){
							ciBounds[year.i,] <- ciBounds[year.i,] + qnorm(ciBoundQs,mean = 0, sd=sqrt(resSigma[varName,varName]))
						} else {
							stop('Noise uncertainty not implemented yet for dependent vars\n')
						}
					}
				}
				location.output.repSample <- file.path(location.output,'repSample',plotWeightType)
				if(file.exists(file.path(location.output.repSample,'subSampleParameterIndices.RDS'))){
					repSampleID <- unname(unlist(readRDS(file.path(location.output.repSample,'subSampleParameterIndices.RDS'))))
					repSampleIDExists <- T
					repSample <- varData[,repSampleID]
				} else {
					repSampleIDExists <- F
					alsoPlotRepSample.vals <- c(F)
					repSample <- NULL
				}
				## save data ####
				cat('saving...')
				plotData <- list()
				plotData$variable <- varName
				plotData$years <- rownames(defRun)
				plotData$CIsToPlot <- CIsToPlot
				plotData$ciBoundQs <- ciBoundQs
				plotData$ciBounds <- ciBounds
				plotData$uncertaintyType <- uncertaintyType
				plotData$means <- means
				plotData$defaultRun <- defRun[[varName]]
				plotData$repSample <- repSample
				plotData$calDat <- calDat[[varName]]
				plotData$varName.orig <- varName.orig
				plotData$plotWeightType <- plotWeightType
				plotData$ciBoundQs.lty <- ciBoundQs.lty
				plotData$ciBoundQs.lwd <- ciBoundQs.lwd
				plotData$ciBoundQs.lcol <- ciBoundQs.lcol
				dir.create(file.path(location.output,location.plots,'CI-plots',
														 paste0(plotWeightType,'Weighted'),'plotData'),F,T)
				saveRDS(plotData,file.path(location.output,location.plots,'CI-plots',
																	 paste0(plotWeightType,'Weighted'),'plotData',
																	 paste0(paste(varName,uncertaintyType,plotWeightType,'weighted',sep='-'),'.RDS')))
				csvExport.df <- data.frame(year=plotData$years,mean=means,defaultRun=defRun[[varName]])
				for(q.i in 1:length(ciBoundQs)){
					csvExport.df[[paste0('Quantile',ciBoundQs[q.i])]] <- ciBounds[,q.i]
				}
				if(repSampleIDExists){
					for(r.i in 1:length(repSampleID)){
						csvExport.df[[paste0('RepSample',r.i)]] <- repSample[r.i,]
					}
				}
				write.csv(csvExport.df,file.path(location.output,location.plots,'CI-plots',
																				 paste0(plotWeightType,'Weighted'),'plotData',
																				 paste0(paste(varName,uncertaintyType,plotWeightType,'weighted',sep='-'),'.csv')))
				cat('drawing...')
				for(years.i in 1:length(yearsToPlot.lst)){
					for(alsoPlotMean in alsoPlotMean.vals){
						for(alsoPlotDefaultRun in alsoPlotDefaultRun.vals){
							for(alsoPlotRepSample in alsoPlotRepSample.vals){
								yearsToPlot <- yearsToPlot.lst[[years.i]]
								location.plots.ci <- file.path(location.output,location.plots,'CI-plots',
																							 paste0(plotWeightType,'Weighted'),
																							 yearsToPlot.names[years.i],
																							 paste(uncertaintyType,
																							 			ifelse(alsoPlotMean,'withMean','withoutMean'),
																							 			ifelse(alsoPlotDefaultRun,'withDefaultRun','withoutDefaultRun'),
																							 			ifelse(alsoPlotRepSample,'withRepSample','withoutRepSample'),sep='-'))
								dir.create(location.plots.ci,F,T)
								png(file.path(location.plots.ci,paste0(paste(varName,plotWeightType,'weighted',sep='-'),'.png')),
										width = plotWidth,height = plotHeight,units = plotUnit,res = plotRes)
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
								legend.pch=c(rep(NA,length(CIsToPlot)+
																 	sum(c(alsoPlotMean&uncertaintyType!='noise uncertainty',alsoPlotDefaultRun))),20)
								legend.col = c(
									if(alsoPlotMean&&uncertaintyType!='noise uncertainty'){mean.col},
									if(alsoPlotDefaultRun){def.col},
									CIsToPlot.lcol,
									calDat.col)
								legend('bottom',legend.text,lty=legend.lty,lwd=legend.lwd,pch=legend.pch,col=legend.col,
											 horiz=T,xpd=T)
								par(mar=c(3.1,4.1,4.1,2.1))
								plot(yearsToPlot,ciBounds[yearsToPlot,medianQIdx],
										 ylim=range(c(ciBounds,calDat[yearsToPlot,varName]),na.rm=T),
										 xlim=range(as.numeric(yearsToPlot)),
										 xaxs='i',
										 type='n',
										 xlab='',
										 ylab=varName.orig,
										 xaxt='n',
										 main=varName.orig)
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
									lines(yearsToPlot,defRun[yearsToPlot,varName],lty=def.lty,lwd=def.lwd,col=def.col)
								}
								if(alsoPlotRepSample & repSampleIDExists){
									for(r.i in 1:length(repSampleID)){
										lines(yearsToPlot,repSample[r.i,yearsToPlot],lty=rs.lty,lwd=rs.lwd,col=rs.col)
									}
								}
								if(!is.null(calDat[[varName]])){
									points(rownames(calDat),calDat[[varName]],
												 col=calDat.col,pch=20)
								}
								box()
								dev.off()
							}
						}
					}
				}
			}
			cat('done\n')
		}
	}
}
