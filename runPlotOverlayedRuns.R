# just for plot specification
plotWeightType <- 'equaly'
dataForOverlayedFiguresFolders <- c('workOutput/figures_JeffEMB/CI-plots/equalyWeighted/plotData/',
																		'workOutput/JeffGDP_figures/CI-plots/equalyWeighted/plotData/')
overlayNames <- c('EMB','GDP driven food demand')
overlayColors <- c('black','blue')

location.plots <- file.path('workOutput','jeffsOverlayedFigures',uncertaintyType)
dir.create(location.plots,F,T)
for(file in list.files(dataForOverlayedFiguresFolders[1])){
	if(sum(file.exists(file.path(dataForOverlayedFiguresFolders,file)))==length(dataForOverlayedFiguresFolders)){
		varName <- tools::file_path_sans_ext(file)
		png(file.path(location.plots,paste0(varName,'.png')),
				width = plotWidth,height = plotHeight,units = plotUnit,res = plotRes)
		plotData <- readRDS(file.path(dataForOverlayedFiguresFolders[1],file))
		yearsToPlot <- plotData$years
		uncertaintyType <- plotData$uncertaintyType
		ciBoundQs <- plotData$ciBoundQs
		ciBounds <- plotData$ciBounds
		means <- plotData$means
		defRun <- plotData$defaultRun
		varName.orig <- plotData$varName.orig
		calDat <- plotData$calDat
		ciBoundQs <- unique(c(rev((1-CIsToPlot)/2),1-(1-CIsToPlot)/2))
		ciBoundQs.lty <- c(rev(CIsToPlot.lty),CIsToPlot.lty[-1])
		ciBoundQs.lwd <- c(rev(CIsToPlot.lwd),CIsToPlot.lwd[-1])
		ciBoundQs.lcol <- c(rev(CIsToPlot.lcol),CIsToPlot.lcol[-1])
		plotData <- readRDS(file.path(dataForOverlayedFiguresFolders[o.i],file))
		yearsToPlot <- plotData$years
		uncertaintyType <- plotData$uncertaintyType
		ciBoundQs <- plotData$ciBoundQs
		ciBounds <- plotData$ciBounds
		means <- plotData$means
		defRun <- plotData$defaultRun
		varName.orig <- plotData$varName.orig
		calDat <- plotData$calDat
		medianQIdx <- which(ciBoundQs==0.5)
		layout(matrix(c(3,2,1),nrow=3),heights = c(0.8,0.1,0.1))
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
		par(mar=c(0,0,0,0))
		plot(0,0,type='n',axes=F)
		legend('bottom',overlayNames,pch=15,col=overlayColors,
					 horiz=T,xpd=T)
		par(mar=c(3.1,4.1,4.1,2.1))
		plot(yearsToPlot,ciBounds[yearsToPlot,medianQIdx],
				 ylim=range(c(ciBounds[yearsToPlot,c(2,length(ciBoundQs)-1)],calDat[yearsToPlot,varName]),na.rm=T),
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
		for(o.i in 1:length(dataForOverlayedFiguresFolders)){
			ciBoundQs.lcol <- rep(overlayColors[o.i],length(ciBoundQs.lcol))
			plotData <- readRDS(file.path(dataForOverlayedFiguresFolders[o.i],file))
			yearsToPlot <- plotData$years
			uncertaintyType <- plotData$uncertaintyType
			ciBoundQs <- plotData$ciBoundQs
			ciBounds <- plotData$ciBounds
			means <- plotData$means
			defRun <- plotData$defaultRun
			varName.orig <- plotData$varName.orig
			calDat <- plotData$calDat
			for(ci.i in length(CIsToPlot.col):1){
				CIsToPlot.col[ci.i] <- adjustcolor(overlayColors[o.i],0.4/ci.i)
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
		}
		for(o.i in 1:length(dataForOverlayedFiguresFolders)){
			ciBoundQs.lcol <- rep(overlayColors[o.i],length(ciBoundQs.lcol))
			plotData <- readRDS(file.path(dataForOverlayedFiguresFolders[o.i],file))
			yearsToPlot <- plotData$years
			uncertaintyType <- plotData$uncertaintyType
			ciBoundQs <- plotData$ciBoundQs
			ciBounds <- plotData$ciBounds
			means <- plotData$means
			defRun <- plotData$defaultRun
			varName.orig <- plotData$varName.orig
			calDat <- plotData$calDat
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
			if(!is.null(calDat[[varName]])){
				points(rownames(calDat),calDat[[varName]],
							 col=calDat.col,pch=20)
			}
		}
		box()
		dev.off()
	}
}
