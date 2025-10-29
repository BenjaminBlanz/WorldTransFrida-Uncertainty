
# plot function ####
# In the first run
# plot unfiltered and apply filter and report back the filtered out IDs
# In the second run 
# plot filtered runs
filterResults <- function(varFile,filterSpec,verbosity=0,polIDsToDrop=NULL,useCluster=T){
	varName <- tools::file_path_sans_ext(varFile)
	if(varName %in% names(filterSpec)){
		if(verbosity>0){cat(sprintf('reading %s...',varName))}
		varDat <- readPerVarFile(file.path(outputFolder,varFile))
		if(verbosity>0){cat('converting to data.table...')}
		varDat <- data.table(varDat)
		if(verbosity>0){cat('determining filtered...')}
		# if the last column is entirely NA this is probably
		# a generated variable, drop that col to not mess with the complete cases
		# filter
		if(sum(is.na(varDat[,ncol(varDat)]))==nrow(varDat)){
			varDat <- varDat[,-ncol(varDat)]
		}
		# if we have ex ante known polIDs to drop, we do not have to filter them
		if(is.null(polIDsToDrop)){
			if(verbosity>0){cat('applying prefilter...')}
			varDat <- varDat[!varDat$polID %in% polIDsToDrop,]
		}
		cat('dropping incomplete and inf...')
		polIDsToDrop <- unique(varDat$polID[!complete.cases(varDat) | !is.finite(varDat[,ncol(varDat)])])
		# years
		years <- colnames(varDat)[-c(1,2)]
		filterFun <- function(year.i){
			year <- years[year.i]
			if(filterSpec[[varName]][1] == 'ltabs'){
				polIDsToDrop <- varDat$polID[abs(varDat[[year]])<filterSpec[[varName]][2]]
			} else if(filterSpec[[varName]][1] == 'gtabs'){
				polIDsToDrop <- varDat$polID[abs(varDat[[year]])>filterSpec[[varName]][2]]
			} else if (filterSpec[[varName]][1] == 'ltval'){
				polIDsToDrop <- varDat$polID[varDat[[year]]<filterSpec[[varName]][2]]
			} else if (filterSpec[[varName]][1] == 'gtval'){
				polIDsToDrop <- varDat$polID[varDat[[year]]>filterSpec[[varName]][2]]
			} else {
				stop('unkown filter spec\n')
			}
			if(verbosity>0){cat('.')}
		}
		if(verbosity>0){cat('processing years')}
		if(!useCluster){
			if(verbosity>0){cat(' not using cluster...')}
			yearPolIDsToDrop <- lapply(1:length(years),filterFun)
		} else {
			if(verbosity>0){cat(' starting cluster...')}
			clFiltering <- makeForkCluster(min(10,detectCores()))
			if(verbosity>0){cat('running...')}
			yearPolIDsToDrop <- parLapply(clFiltering,1:length(years),filterFun)
			stopCluster(clFiltering)
		}
		polIDsToDrop <- unique(c(polIDsToDrop,unlist(yearPolIDsToDrop)))
	} else {
		polIDsToDrop <- c()
	}
	if(verbosity>0){cat('done\n')}
	return(polIDsToDrop)
}

# plotType:
#  0: Boxplot of medians
#  1: Flat area of medians and flat area of sow extending beyond median
#  2: contourplot of medians
#  3: contourplot of medians overlayed with contour of sow
plotPolResults <- function(varFile,polIDsToDrop=NULL,figuresFolder=NULL,
													 plotType=1,
													 verbosity=0){
	if(is.null(figuresFolder)){
		figuresFolder <- file.path(location.output,'figures',paste0('plotType',plotType))
	}
	dir.create(figuresFolder,showWarnings = F,recursive = T)
	varName <- tools::file_path_sans_ext(varFile)
	if(verbosity>0){cat(sprintf('reading %s...',varName))}
	varDat <- readPerVarFile(file.path(outputFolder,varFile))
	varUnit <- varsMeta$Unit[varsMeta$cleanName==varName]
	varFullName <- varsMeta$FRIDA.FQN[varsMeta$cleanName==varName]
	cat('applying filters...')
	varDat <- varDat[!varDat$polID %in% polIDsToDrop,]
	if(verbosity>0){cat('plotting...')}
	years <- as.numeric(colnames(varDat)[3:ncol(varDat)])
	png(file.path(figuresFolder,paste0(varName,'.png')),
			width = plotWidth,height = plotHeight,units = plotUnit,res = plotRes)
	ylims <- quantile(varDat[,seq(3,ncol(varDat),3)],probs=plot.relyrange,na.rm=T)
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
							axes=F)
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
			histDat <- hist(varDat[[as.character(year)]][varDat$sowID==5],
											breaks=breaks,plot=F)
			counts[,year.i] <- histDat$counts[-c(1,length(histDat$counts))]
			if(plotType==3){
				histDat <- hist(varDat[[as.character(year)]],
												breaks=breaks,plot=F)
				countsSOW[,year.i] <- histDat$counts[-c(1,length(histDat$counts))]
			}
			if(verbosity>1){cat('.')}
		}
		breaks <- breaks[-c(1,length(breaks))]
		breaksMids <- breaks[1:(length(breaks)-1)]+
			(breaks[2:(length(breaks))]-breaks[1:(length(breaks)-1)])/2
		logmax <- log(max(counts))
		levels <- exp(seq(0,logmax,length.out=plot.numColLevels))
		if(plotType==2){
			filled.contour(years,breaksMids,t(counts),
										 xlim=range(years),
										 ylim=ylims,
										 xlab='year',
										 xaxs='i',
										 levels=levels,
										 main=varFullName,ylab=varUnit,
										 col = rev(paletteer_c(plot.palletteName, plot.numColLevels-1))
										 key.title = 'Number of policies')
		} else if(plotType==3){
			filled.contour(years,breaksMids,t(counts),
										 xlim=range(years),
										 ylim=ylims,
										 xlab='year',
										 xaxs='i',
										 levels=levels,
										 main=varFullName,ylab=varUnit,
										 col = rev(paletteer_c(plot.palletteNameSOW, plot.numColLevels-1))
										 key.title = 'Number of policies',
										 plot.axes={
										 	filled.contour(years,breaksMids,t(counts),
										 								 xlim=range(years),
										 								 ylim=ylims,
										 								 xlab='year',
										 								 xaxs='i',
										 								 levels=levels,
										 								 main=varFullName,ylab=varUnit,
										 								 col = rev(paletteer_c(plot.palletteName, plot.numColLevels-1))
										 								 key.title = 'Number of policies')
										 })
		}
	}
	dev.off()
	if(verbosity>0){cat('done\n')}
}
parPlotPolResults<-function(i,varsFiles,polIDsToDrop,figuresFolder=NULL,plotType=1){
	plotPolResults(varsFiles[i],polIDsToDrop,figuresFolder,plotType)	
}
