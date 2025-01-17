thing <- 'gdp_real_gdp_in_2021c-fit uncertainty-equaly-weighted'
#thing <- 'gdp_real_gdp_growth_rate-fit uncertainty-equaly-weighted'

# read Data
yesFB <- readRDS(file.path(
	'workOutput',
	'N-5e+05-ChS-100-LCR-1000-IgB-FALSE-FrB-FALSE-KcE-FALSE-Sym-Min-AAZ-FALSE-CFB-TRUE-Pol-policy_EMB',
	'figures','CI-plots','equalyWeighted','plotData',
	paste0(thing,'.RDS')))
yesFB.temp <- readRDS(file.path(
	'workOutput',
	'N-5e+05-ChS-100-LCR-1000-IgB-FALSE-FrB-FALSE-KcE-FALSE-Sym-Min-AAZ-FALSE-CFB-TRUE-Pol-policy_EMB',
	'figures','CI-plots','equalyWeighted','plotData',
	'energy_balance_model_surface_temperature_anomaly-fit uncertainty-equaly-weighted.RDS'))
noFB <- readRDS(file.path(
	'workOutput',
	'N-5e+05-ChS-100-LCR-1000-IgB-FALSE-FrB-FALSE-KcE-FALSE-Sym-Min-AAZ-FALSE-CFB-FALSE-Pol-policy_EMB',
	'figures','CI-plots','equalyWeighted','plotData',
	paste0(thing,'.RDS')))
noFB.temp <- readRDS(file.path(
	'workOutput',
	'N-5e+05-ChS-100-LCR-1000-IgB-FALSE-FrB-FALSE-KcE-FALSE-Sym-Min-AAZ-FALSE-CFB-FALSE-Pol-policy_EMB',
	'figures','CI-plots','equalyWeighted','plotData',
	'energy_balance_model_surface_temperature_anomaly-fit uncertainty-equaly-weighted.RDS'))

# y1: value at time 1
# y2: value at time 2
# N:  Number of time units between time 1 and time 2
# returns the average growth rate over that time
avgGrowthRate <- function(y1,y2,N){
	(y2/y1)^(1/(N-1)) - 1
}

thing <- gsub('-fit uncertainty-equaly-weighted','',thing)

# data post proc
temps <- list()
yesFB.avgGr <- yesFB
noFB.avgGr <- noFB
for(t in 1:(length(yesFB$years))){
	yesFB.avgGr$ciBounds[t,] <- avgGrowthRate(yesFB$ciBounds[1,],yesFB$ciBounds[t,],t)
	yesFB.avgGr$means[t] <- avgGrowthRate(yesFB$means[1],yesFB$means[t],t)
	yesFB.avgGr$defaultRun[t] <- avgGrowthRate(yesFB$defaultRun[1],yesFB$defaultRun[t],t)
	noFB.avgGr$ciBounds[t,] <- avgGrowthRate(noFB$ciBounds[1,],noFB$ciBounds[t,],t)
	noFB.avgGr$means[t] <- avgGrowthRate(noFB$means[1],noFB$means[t],t)
	noFB.avgGr$defaultRun[t] <- avgGrowthRate(noFB$defaultRun[1],noFB$defaultRun[t],t)
}
relDiff <- diff <- yesFB
diff.avgGr <- diff
for(field in c('ciBounds','means','defaultRun')){
	diff[[field]] <- noFB[[field]]-yesFB[[field]]
	relDiff[[field]] <- diff[[field]]/noFB[[field]]
	diff.avgGr[[field]] <- noFB.avgGr[[field]]-yesFB.avgGr[[field]]
	temps[[field]] <- data.frame(years=yesFB.temp$years,yesFB=yesFB.temp[[field]],noFB=noFB.temp[[field]]) 
}


# test the average growth rate
plot(yesFB.temp$defaultRun,diff.avgGr$defaultRun*100,
		 xlab='Surface Temperature Anomaly °C',
		 ylab='Difference in average growth rates')
plot(yesFB$defaultRun,ylim=range(yesFB$defaultRun,noFB$defaultRun))
points(noFB$defaultRun,col='red')
for(t in 2:length(yesFB$years)){
	color <- rep('gray',length(yesFB$years))
	color[1:t] <- 'black'
	ys <- c(yesFB$defaultRun[1])
	for(ti in 2:length(yesFB$years)){
		ys[ti] <- ys[ti-1]*(1+yesFB.avgGr$defaultRun[t])#+diff.avgGr$defaultRun[t])
		lines(c(ti-1,ti),c(ys[ti-1],ys[ti]),col=color[ti])
	}
}



# plot temp ####
plot(temps$defaultRun$years,temps$defaultRun$yesFB,type='l',
		 ylab='Surface Temperature Anomaly °C',
		 xlab='year',
		 ylim=range(temps$defaultRun[,-1]))
lines(temps$defaultRun$years,temps$defaultRun$noFB,col='red')
legend('topleft',
			 legend = c('with feedbacks','without feedbacks'),
			 lty=1,
			 col=c(1,'red'))

# plot diff ####
plot(noFB.temp$year,yesFB$defaultRun,type='l',
		 xlab='year',ylab=thing,
		 ylim=range(yesFB$defaultRun,noFB$defaultRun))
lines(noFB.temp$year,noFB$defaultRun,col='red')
legend('topright',
			 legend = c('with feedbacks','without feedbacks'),
			 lty=1,
			 col=c(1,'red'))

if(thing=="gdp_real_gdp_growth_rate"){
	plot(yesFB.temp$year,diff$ciBounds[,'0.5']*100,type='l',
			 xlab='year',ylab=thing,
			 main=sprintf('Difference in %s\n between yesFB and noFB',thing))
	
	
	plot(yesFB.temp$ciBounds[,'0.5'],diff$ciBounds[,'0.5']*100,type='l',
			 xlab='Surface Temperature Anomaly °C',ylab=paste(thing,'difference'),
			 main='Difference in thing between yesFB and noFB')
}

# annual loss ####

# annual gdp change ####
yesFB.growth <- yesFB
noFB.growth <- noFB
diff.growth <- diff
for(t in 1:(length(yesFB$years)-1)){
	yesFB.growth$ciBounds[t,] <- yesFB$ciBounds[t+1,] - yesFB$ciBounds[t,]
	yesFB.growth$means[t] <- yesFB$means[t+1] - yesFB$means[t]
	yesFB.growth$defaultRun[t] <- yesFB$defaultRun[t+1] - yesFB$defaultRun[t]
	noFB.growth$ciBounds[t,] <- noFB$ciBounds[t+1,] - noFB$ciBounds[t,]
	noFB.growth$means[t] <- noFB$means[t+1] - noFB$means[t]
	noFB.growth$defaultRun[t] <- noFB$defaultRun[t+1] - noFB$defaultRun[t]
	diff.growth$ciBounds[t,] <- diff$ciBounds[t+1,] - diff$ciBounds[t,]
	diff.growth$means[t] <- diff$means[t+1] - diff$means[t]
	diff.growth$defaultRun[t] <- diff$defaultRun[t+1] - diff$defaultRun[t]
}
t <- length(yesFB$years)
yesFB.growth$ciBounds[t,] <- yesFB.growth$means[t] <- yesFB.growth$defaultRun[t] <- 
	noFB.growth$ciBounds[t,] <- noFB.growth$means[t] <- noFB.growth$defaultRun[t] <-
	diff.growth$ciBounds[t,] <- diff.growth$means[t] <- diff.growth$defaultRun[t] <- NA
# points(yesFB.growth$defaultRun+diff.growth$defaultRun,col='blue')
par(pch=20)
plot(0,0,
		 xlab='Surface Temperature Anomaly °C',
		 ylab='Annual % GDP loss',
		 xlim=c(0,7),ylim=c(-10,80),yaxs='i')
ci.cols <- list('0.025'=hsv(0,1,1),
								'0.25'=hsv(0,1,0.5),
								'0.5'=hsv(0,0,0),
								'0.75'=hsv(240/359,1,0.5),
								'0.975'=hsv(240/359,1,1))
for(ci in c('0.025','0.25','0.5','0.75','0.975')){
	# points(yesFB.temp$ciBounds[,ci],100*diff.growth$ciBounds[,ci]/yesFB$ciBounds[,ci],
	points(yesFB.temp$ciBounds[,ci],100*(noFB.growth$ciBounds[,ci]-yesFB.growth$ciBounds[,ci])/noFB$ciBounds[,ci],
			 col=ci.cols[[ci]])
}
legend('topleft',
			 legend=paste('FRIDA',c('0.025','0.25','0.5','0.75','0.975'),'quantile'),
			 pch=par('pch'),
			 col=unlist(ci.cols))

# over time ####
par(pch=20)
plot(0,0,
		 xlab='Surface Temperature Anomaly °C',
		 ylab='Annual % GDP loss',
		 xlim=as.numeric(range(yesFB$years)),ylim=c(-10,80),yaxs='i')
ci.cols <- list('0.025'=hsv(0,1,1),
								'0.25'=hsv(0,1,0.5),
								'0.5'=hsv(0,0,0),
								'0.75'=hsv(240/359,1,0.5),
								'0.975'=hsv(240/359,1,1))
for(ci in c('0.025','0.25','0.5','0.75','0.975')){
	# points(yesFB.temp$ciBounds[,ci],100*diff.growth$ciBounds[,ci]/yesFB$ciBounds[,ci],
	points(yesFB$years,100*(noFB.growth$ciBounds[,ci]-yesFB.growth$ciBounds[,ci])/yesFB$ciBounds[,ci],
				 col=ci.cols[[ci]])
}
legend('topleft',
			 legend=paste('FRIDA',c('0.025','0.25','0.5','0.75','0.975'),'quantile'),
			 pch=par('pch'),
			 col=unlist(ci.cols))
# difference in average growth rate ####
plot(temps$defaultRun$yesFB,diff.avgGr$defaultRun*100,
		 xlab='Year',
		 ylab='R')


