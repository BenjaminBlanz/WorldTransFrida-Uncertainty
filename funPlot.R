# funPlotDat ####
funPlotDat <- function(calDat,calDat.impExtrValue=NULL,defDat=NULL){
	if(is.null(defDat)){
		defDat <- calDat
		defDat[,] <- NA
	}
	if(is.null(alDat.impExtrValue)){
		calDat.impExtrValue <- calDat
		calDat.impExtrValue[,] <- F
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
		for(i in 1:plotCols){
			oldMar <- par('mar')
			par(mar=c(0,oldMar[2],0,oldMar[4]))
			plot(0,0,type='n',axes=F,
					 xlim=as.numeric(range(rownames(calDat))),
					 ylim=c(0,1))
			# abline(h=0)
			axis(1,pos=1)
			if(i == 1){
				mtext(sprintf(paste('Gray areas do not have complete cases.',
														'Red lines highlight the vars which most limit the complete cases window.',
														'Red crosses are imputed or extrapolated values.',
														'Blue line is prior model calibration.',
														'Complete cases: %i'),nrow(calDat)-length(incompleteIdc)),
							3,line = -4,adj=0)
			}
			par(mar=oldMar)
		}
		for(i in 1:ncol(calDat)){
			plot(rownames(calDat),calDat[[i]],type='n',
					 xaxt='n',yaxt='n',
					 ylim=range(c(calDat[[i]],defDat[[i]]),na.rm=T))
			for(year in incompleteYears){
				year <- as.numeric(year)
				rect(year-0.5,par('usr')[3],year+0.5,par('usr')[4],
						 density=-1,col='gray',border=NA)
			}
			box()
			# add label specifying var indext to top left
			text(years[1],max(calDat[[i]],defDat[[i]],na.rm=T),i,adj=c(0,1))
			# highlight if this var has its last obs on the boundary of incomplete obs
			validRange <- funValidRange(calDat[[i]])
			if(firstCompleteIdx!=1 && validRange[1]==firstCompleteIdx){
				abline(v=years[validRange[1]-1]+0.5,col='red',lwd=3)
			}
			if(lastCompleteIdx!=nrow(calDat) && validRange[2]==lastCompleteIdx){
				abline(v=years[validRange[2]+1]-0.5,col='red',lwd=3)
			}
			# add default model line
			lines(rownames(defDat),defDat[[i]],col='blue')
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
		}
		padding <- rep(NA,(plotRows*plotCols)-ncol(calDat))
		plotLocations <- matrix(c(1:ncol(calDat),padding),byrow=T,nrow=plotRows)
		return(plotLocations)
}
