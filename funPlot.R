# funPlotDat ####
funPlotDat <- function(calDat,calDat.impExtrValue=NULL,defDat=NULL,yaxPad=0.04,
											 highlightConstrainingVars=F,shadowIncompleteYears=F){
	noImpExtrValues <- F
	if(is.null(calDat.impExtrValue)){
		calDat.impExtrValue <- calDat
		calDat.impExtrValue[,] <- F
		noImpExtrValues <- T
	}
	completeIdc <- which(complete.cases(calDat))
	incompleteIdc <- which(!complete.cases(calDat))
	incompleteYears <- as.numeric(rownames(calDat)[incompleteIdc])
	firstCompleteIdx <- completeIdc[1]
	lastCompleteIdx <- completeIdc[length(completeIdc)]
	years <- as.numeric(rownames(calDat))
	ncols <- ncol(calDat)
	sqrtNcols <- sqrt(ncols)
	plotCols <- round(sqrtNcols)
	plotRows <- ceiling(sqrtNcols)
	rm(sqrtNcols)
	# par(mfrow=c(plotRows,plotCols))
	heights <- c(rep(2,plotRows),0.5)
	heights <- heights/sum(heights)
	layout(
		rbind(
			matrix((plotCols+1):(plotCols+plotRows*plotCols),
						 byrow=T,nrow=plotRows),
			t(1:plotCols)),
		widths = 1,
		heights = heights)
	par(mar=c(0,1,1,0.1))
	xlims <- c(min(as.numeric(rownames(calDat)))-0.5,
						 max(as.numeric(rownames(calDat)))+0.5)
	for(i in 1:plotCols){
		oldMar <- par('mar')
		par(mar=c(0,oldMar[2],0,oldMar[4]))
		plot(0,0,type='n',axes=F,
				 xlim=xlims,
				 xaxs='i',
				 ylim=c(0,1))
		# abline(h=0)
		axis(1,pos=1)
		if(i == 1){
			mtext(sprintf(paste(ifelse(!shadowIncompleteYears,
																 'Gray areas do not have complete cases.',
																 ''),
													ifelse(highlightConstrainingVars&!shadowIncompleteYears,
																 'Red lines highlight the vars which most limit the complete cases window.',
																 ''),
													ifelse(!noImpExtrValues,
																 'Red crosses are imputed or extrapolated values.',
																 ''),
													ifelse(!is.null(defDat),
																 'Blue line is prior model calibration.',
																 ''),
													ifelse(!shadowIncompleteYears,
																 'Complete cases: %i',
																 ''))
										,nrow(calDat)-length(incompleteIdc)),
						3,line = -4,adj=0,cex=0.8)
		}
		par(mar=oldMar)
	}
	for(i in 1:ncol(calDat)){
		yrange <- range(calDat[[i]],na.rm=T)
		ylims <- c(yrange[1]-abs(diff(yrange))*yaxPad,
							 yrange[2]+abs(diff(yrange))*yaxPad)
		plot(rownames(calDat),calDat[[i]],type='n',
				 xaxt='n',yaxt='n',
				 xaxs='i',yaxs='i',
				 xlim=xlims,
				 ylim=ylims)
		if(shadowIncompleteYears){
			for(year in incompleteYears){
				year <- as.numeric(year)
				rect(year-0.5,par('usr')[3],year+0.5,par('usr')[4],
						 density=-1,col='gray',border=NA)
			}
			box()
		}
		if(highlightConstrainingVars){
			# highlight if this var has its last obs on the boundary of incomplete obs
			validRange <- funValidRange(calDat[[i]])
			if(firstCompleteIdx!=1 && validRange[1]==firstCompleteIdx){
				abline(v=years[validRange[1]-1]+0.5,col='red',lwd=3)
			}
			if(lastCompleteIdx!=nrow(calDat) && validRange[2]==lastCompleteIdx){
				abline(v=years[validRange[2]+1]-0.5,col='red',lwd=3)
			}
		}
		# add default model line
		if(!is.null(defDat)){
			lines(rownames(defDat),defDat[[i]],col='blue')
		}
		# add obs points
		points(rownames(calDat)[!calDat.impExtrValue[,i]],
					 calDat[[i]][!calDat.impExtrValue[,i]])
		# add imputed points in red
		if(imputeMissingVars||extrapolateMissingVarMethod!='n'){
			points(rownames(calDat)[calDat.impExtrValue[,i]],
						 calDat[[i]][calDat.impExtrValue[,i]],
						 col='red',pch=3)
		}
		# add axis ticks
		axis(1,labels=F,tcl=0.3)
		yAxVals <- axis(2,labels=F,tcl=0.3)
		axTextCex <- 1
		# min yax label
		adjVal <- c(0.5,0)
		yPosVal <- min(yAxVals)
		if(ydev2in(yPosVal-par('usr')[3]) < 
			 strwidth(as.character(min(yAxVals)),
			 				 # font=par('font'),
			 				 # family = 'Serif',
			 				 units='inch',
			 				 cex=axTextCex)){
			adjVal <- c(0,0)
			yPosVal <- par('usr')[3]
		}
		text(par('usr')[1]-0.01*diff(par('usr')[1:2]),
				 yPosVal,
				 min(yAxVals),
				 xpd=T,cex=axTextCex,srt=90,adj=adjVal)
		# max yax label
		adjVal <- c(0.5,0)
		yPosVal <- max(yAxVals)
		if(ydev2in(yPosVal-par('usr')[4]) < 
			 strwidth(as.character(max(yAxVals)),
			 				 font=par('font'),
			 				 units='inch',
			 				 cex=axTextCex)){
			#0.1*(par('usr')[4]-par('usr')[3])){
			adjVal <- c(1,0)
			yPosVal <- par('usr')[4]
		}
		text(par('usr')[1]-0.01*diff(par('usr')[1:2]),
				 yPosVal,
				 max(yAxVals),
				 xpd=T,cex=axTextCex,srt=90,adj=adjVal)
		# highlight the labeled points
		points(rep(par('usr')[1],2),range(yAxVals),col=1,xpd=F,pch=18,cex=2.1)
		# title
		text(mean(par('usr')[1:2]),par('usr')[4]+diff(par('usr')[3:4])*0.01,
				 colnames(calDat)[i],
				 cex=0.8,xpd=T,adj=c(0.5,0))
		# add label specifying var indext to top left
		text(years[1],yrange[2]+abs(diff(yrange))*yaxPad,i,adj=c(0,1))
	}
	padding <- rep(NA,(plotRows*plotCols)-ncol(calDat))
	plotLocations <- matrix(c(rep(NA,plotCols),1:ncol(calDat),padding),byrow=T,ncol = plotCols)
	return(plotLocations)
}


# funPlotParRanges ####
funPlotParRangesLikelihoods <- function(sampleParms,sampleParms.orig=NULL,
																				samplePoints=NULL,like=NULL,yaxPad=0.04,
																				outputFilePath=NULL,
																				outputFileName='parameterLogLikelihoods.png',
																				baseLike=NULL,
																				includeZeroYVal=FALSE,logY=F,ylim=NULL,
																				minLike=NULL, parallelPlot = F,
																				plotWidth = 10,
																				plotHight = 10,
																				plotRes = 150,
																				plotUnits = 'cm'){
	if(is.null(samplePoints)){
		samplePointBase <- seq(0,1,length.out=10)
		samplePoints <- array(rep(samplePointBase,nrow(sampleParms)),
													dim=c(length(samplePointBase),nrow(sampleParms)))
		samplePoints <- funStretchSamplePoints(samplePoints,sampleParms,F)
		colnames(samplePoints) <- sampleParms[,1]
	}
	if(is.null(sampleParms.orig)){
		sampleParms.orig <- sampleParms
	}
	if(is.null(like)){
		like <- array(rep(1,nrow(samplePoints)),dim=c(nrow(samplePoints),nrow(sampleParms)))
	}
	if(is.null(baseLike)){
		baseLike <- max(like,na.rm=T)
	}
	if(is.null(minLike)){
		minLike <- baseLike - abs(0.99*baseLike)
	}
	like.range <- range(c(like,baseLike),na.rm=T)
	
	sqrtNplots <- sqrt(nrow(sampleParms))
	plotCols <- round(sqrtNplots)
	plotRows <- ceiling(sqrtNplots)
	if(parallelPlot){
		if(!is.null(savePlotFilePath)) {
			stop('Parallel plotting requires a file path.\n')
		}
		if(is.null(cl)){
			stop('Parallel plotting requirs a cluster')
		}
		cl.local <- makeForkCluster(detectCores())
		parLapply(cl,funPlotParRangesLikelihoodsI,
							sampleParms=sampleParms,sampleParms.orig=sampleParms.orig,
							samplePoints=samplePoints,like=like,yaxPad=yaxPad,
							outputFilePath=outputFilePath,
							outputFileName=outputFileName,
							baseLike=baseLike,
							includeZeroYVal=includeZeroYVal,logY=logY,ylim=ylim,
							like.range=like.range,
							minLike=minLike, parallelPlot=parallelPlot,
							plotWidth=plotWidth,
							plotHight=plotHight,
							plotRes=plotRes,
							plotUnits=plotUnits)
	} else {
		if(!is.null(outputFilePath)){
			outputFileName <- tools::file_path_sans_ext(outputFileName)
			png(file.path(savePlotFilePathName,paste0(outputFileName,'.png')),
					width=plotCols*plotWidth,
					height=plotRows*plotHight,
					res=plotRes,units=plotUnits)
		}
		par(mfrow=c(plotRows,plotCols),mar=c(1,0.5,1,0.5),mgp=c(1,0.5,0))
		cat('\n')
		for(i in 1:nrow(sampleParms)){
			funPlotParRangesLikelihoodsI(i,sampleParms,sampleParms.orig,
																	 samplePoints,like,yaxPad,
																	 outputFilePath=NULL,
																	 outputFileName,
																	 baseLike,
																	 includeZeroYVal,logY,ylim,
																	 like.range,
																	 minLike, parallelPlot,
																	 plotWidth,
																	 plotHight,
																	 plotRes,
																	 plotUnits)
		}
		if(!is.null(savePlotFilePath)){
			dev.off()
		}
	}
	cat(sprintf('\r Plotted %i parms.               ',nrow(sampleParms)))
}

funPlotParRangesLikelihoodsI <- function(i,sampleParms,sampleParms.orig,
																				 samplePoints,like,yaxPad=0.04,
																				 outputFilePath=NULL,
																				 outputFileName='parameterLogLikelihoods.png',
																				 baseLike=NA,
																				 includeZeroYVal=FALSE,logY=F,ylim=NULL,
																				 like.range=NULL,
																				 minLike=NULL, parallelPlot = F,
																				 plotWidth = 10,
																				 plotHight = 10,
																				 plotRes = 150,
																				 plotUnits = 'cm'){
	cat(sprintf('\r Plotting parm %i of %i',i,nrow(sampleParms)))
	if(!is.null(outputFilePath)){
		dir.create(file.path(outputFilePath,outputFileName),F,T)
		outputFileName <- tools::file_path_sans_ext(outputFileName)
		if(parallelPlot){
			outputFile <- file.path(outputFilePath,outputFileName,
															paste0(sampleParms$Variable[i],'.png'))
		}else{
			outputFile <- file.path(outputFilePath,paste0(outputFileName,'.png'))
		}
		png(outputFile,
				width=plotWidth,
				height=plotHight,
				res=plotRes,units=plotUnits)
	}
	if(is.null(like.range)){
		like.range <- range(c(like,baseLike),na.rm=T)
	}
	if(is.null(minLike)){
		minLike <- min(like,na.rm=T)
	}
	if(!is.null(ylim)){
		yrange <- ylim
	} else {
		if(includeZeroYVal){
			yrange <- c(0,like.range[2])
		}else{
			yrange <- c(max(like.range[1],minLike),like.range[2])
		}
	}
	if(logY){
		yrange <- yrange-like.range[1]+1
		like <- like-like.range[1]
		baseLike <- baseLike-like.range[1]
	}
	plot(1,1,type='n',
			 axes=F,
			 xaxs='r',xlab='',ylab='',
			 log=if(logY){'y'}else{''},
			 xlim = c(sampleParms.orig$Min[i],sampleParms.orig$Max[i]),
			 xaxs = 'i',
			 ylim = yrange)
	axis(1,padj=-1,at = c(sampleParms.orig$Min[i]),hadj=0,cex.axis=0.8)
	axis(1,padj=-1,at = c(sampleParms.orig$Max[i]),hadj=1,cex.axis=0.8)
	box(col='gray')
	abline(v=c(sampleParms.orig$Min[i],sampleParms.orig$Max[i]),col='gray')
	axis(1,padj=-1,at = sampleParms$Value[i],labels = '',hadj=0,cex.axis=0.8,col='blue')
	axis(1,padj=-1,at = sampleParms.orig$Value[i],labels = '',hadj=0,cex.axis=0.8,col='red')
	if(sampleParms.orig$Min[i]<sampleParms$Min[i]){
		abline(v=sampleParms$Min[i],lty=2,col='red')
	}
	if(sampleParms$Max[i]<sampleParms.orig$Max[i]){
		abline(v=sampleParms$Max[i],lty=2,col='red')
	}
	points(samplePoints[which(!is.na(like)),i],like[which(!is.na(like))],pch=20,cex=0.5)
	points(sampleParms$Value[i],baseLike,pch=20,cex=1,col='red')
	# add label specifying var indext to top left
	text(mean(par('usr')[1:2]),par('usr')[4]+diff(par('usr')[3:4])*0.01,
			 sampleParms$Variable[i],
			 cex=0.8,xpd=T,adj=c(0.5,0))
	text(sampleParms.orig$Min[i],
			 yrange[2]+abs(diff(yrange))*yaxPad,i,adj=c(0,1),col='red')
}
