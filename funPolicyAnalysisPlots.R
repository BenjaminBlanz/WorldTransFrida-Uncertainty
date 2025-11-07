# plotType:
#  0: Boxplot of medians
#  1: Flat area of medians and flat area of sow extending beyond median
#  2: contourplot of medians
#  3: contourplot of medians overlayed with contour of sow
#  
parPlotPolResults<-function(i,varsFiles,polIDsToDrop,funFigFolder=NULL,colLevels=NULL,plotType=1,verbosity=0,
														baselinePlotProps=NULL,selPolID=NULL){
	if(!is.null(baselinePlotProps)){
		ylims <- baselinePlotProps[[i]]$ylims
	} else {
		ylims <- NULL
	}
	return(plotPolResults(varsFiles[i],polIDsToDrop=polIDsToDrop,
												funFigFolder=funFigFolder,
												plotType=plotType,
												colLevels=colLevels,
												verbosity=verbosity,
												ylims=ylims,
												selPolID=selPolID))
}
plotPolResults <- function(varFile,polIDsToDrop=NULL,funFigFolder=NULL,
													 plotType=1,
													 colLevels=NULL,
													 verbosity=9,
													 ylims=NULL,
													 selPolID=NULL){
	if(is.null(funFigFolder)){
		funFigFolder <- file.path(location.output,'figures',paste0('plotType',plotType))
	}
	dir.create(funFigFolder,showWarnings = F,recursive = T)
	varName <- tools::file_path_sans_ext(varFile)
	if(verbosity>0){cat(sprintf('reading %s...',varName))}
	varDat <- readPerVarFile(file.path(outputFolder,varFile))
	varUnit <- varsMeta$Unit[varsMeta$cleanName==varName]
	varFullName <- varsMeta$FRIDA.FQN[varsMeta$cleanName==varName]
	cat('applying filters...')
	varDat <- varDat[!varDat$polID %in% polIDsToDrop,]
	if(verbosity>0){cat('plotting...')}
	years <- as.numeric(colnames(varDat)[3:ncol(varDat)])
	png(file.path(funFigFolder,paste0(varName,'.png')),
			width = plotWidth,height = plotHeight,units = plotUnit,res = plotRes)
	if(is.null(ylims)){
		ylims <- quantile(varDat[,seq(3,ncol(varDat),3)],probs=plot.relyrange,na.rm=T)
	}
	if(plotType==0){
		plot(0,0,type='n',
				 xlim=range(years),
				 ylim=ylims,
				 xlab='year',
				 xaxs='i',
				 main=varFullName,ylab=varUnit)
		grid(col='gray',lwd=2)
		for(year in years){
			boxplot(varDat[[as.character(year)]][varDat$sowID==5],at=year,width = 0.8,add = T,pch='.',lwd=1,
							axes=F,range=0)
		}
	} else if(plotType==1){
		plot(0,0,type='n',
				 xlim=range(years),
				 ylim=ylims,
				 xlab='year',
				 xaxs='i',
				 main=varFullName,ylab=varUnit)
		grid(col='gray',lwd=2)
		maxBoundUp <- c() # area reached by the any SOW
		maxBoundDown <- c() # lower bound of area reached by any SOW
		medianBoundUp <- c() # upper bound of median results
		medianBoundDown <- c() # lower bound of median results
		for(year.i in 1:length(years)){
			year <- years[year.i]
			maxBoundUp[year.i] <- max(varDat[[as.character(year)]],na.rm=T)
			maxBoundDown[year.i] <- min(varDat[[as.character(year)]],na.rm=T)
			medianBoundUp[year.i] <- max(varDat[[as.character(year)]][varDat$sowID==5],na.rm=T)
			medianBoundDown[year.i] <- min(varDat[[as.character(year)]][varDat$sowID==5],na.rm=T)
			if(verbosity>1){cat('.')}
		}
		maxBoundUpLim <- maxBoundUp
		maxBoundUpLim[maxBoundUp>ylims[2]] <- ylims[2]
		maxBoundDownLim <- maxBoundDown
		maxBoundDownLim[maxBoundDown<ylims[1]] <- ylims[1]
		medianBoundUpLim <- medianBoundUp
		medianBoundUpLim[medianBoundUp>ylims[2]] <- ylims[2]
		medianBoundDownLim <- medianBoundDown
		medianBoundDownLim[medianBoundDown<ylims[1]] <- ylims[1]
		polygon(c(years,rev(years)),c(maxBoundUp,rev(maxBoundDown)),
						border = plot.lcol[1],col = plot.col[1],
						lty=plot.lty[1])
		polygon(c(years,rev(years)),c(medianBoundUp,rev(medianBoundDown)),
						border = plot.lcol[2],col = plot.col[2],
						lty=plot.lty[2])
	} else if(plotType==2 || plotType==3){
		# to build contours we need density values
		# create histogram at each year, use density of hist. 
		# ensure that all hists have the same breaks
		breaks <- c(-Inf,seq(ylims[1],ylims[2],length.out=1001),Inf)
		counts <- array(NA,dim=c(length(breaks)-3,length(years)))
		countsSOW <- counts
		for(year.i in 1:length(years)){
			year <- years[year.i]
			if(plotType==2){
				histDat <- hist(varDat[[as.character(year)]][varDat$sowID==5],
												breaks=breaks,plot=F)
			} else {
				histDat <- hist(varDat[[as.character(year)]],
												breaks=breaks,plot=F)
			}
			counts[,year.i] <- histDat$counts[-c(1,length(histDat$counts))]
			if(verbosity>1){cat('.')}
		}
		breaks <- breaks[-c(1,length(breaks))]
		breaksMids <- breaks[1:(length(breaks)-1)]+
			(breaks[2:(length(breaks))]-breaks[1:(length(breaks)-1)])/2
		if(is.null(colLevels)){
			logmax <- log(max(counts))
			levels <- exp(seq(0,logmax,length.out=plot.numColLevels))
		} else {
			levels <- colLevels
		}
		filled.contour(years,breaksMids,t(counts),
									 xlim=range(years),
									 ylim=ylims,
									 xlab='year',
									 xaxs='i',
									 levels=levels,
									 main=varFullName,ylab=varUnit,
									 col = rev(paletteer_c(ifelse(plotType==2,plot.palletteName,plot.palletteNameSOW),
									 											plot.numColLevels-1)),
									 key.title = 'Number of policies',
									 plot.axes = {
									 	grid(col=gray(0.5,0.8),lwd=2,lty=2)
									 	if(!is.null(selPolID)){
									 		lines(years,varDat[varDat$polID==selPolID & varDat$sowID==selectedRunSpec$sow,as.character(years)])
									 	}
									 	}
									 )
	}
	dev.off()
	if(verbosity>0){cat('done\n')}
	retlist <- list()
	retlist$ylims <- ylims
	return(retlist)
}
