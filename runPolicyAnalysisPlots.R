# plot policy analysis output
source('initialise.R')
source('configPolicyAnalysis.R')
figuresFolder <- file.path(location.output,'figures')
dir.create(figuresFolder,showWarnings = F,recursive = T)
plt.w <- 30
plt.h <- 30
plt.r <- 300
plt.u <- 'cm'

filters <- list()
filters$inf <- 0.4
filters$gdpgr <- 0.4
filters$sta1d <- 0.4

# read 
vars <- list()
vars[['sta']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','energy_balance_model_surface_temperature_anomaly.RDS'))
vars[['pop']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','demographics_population.RDS'))
vars[['gdp']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','gdp_real_gdp_in_2021c.RDS'))
vars[['inf']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','inflation_inflation_rate.RDS'))
vars[['rdu']] <- readPerVarFile(file = file.path(location.output,'detectedParmSpace','PerVarFiles-RDS','gdp_future_year_in_recession.RDS'))

vars$sta1d <- vars$sta
vars$sta1d[,3:(ncol(vars$sta1d)-1)] <- vars$sta1d[,4:(ncol(vars$sta1d))] - vars$sta1d[,3:(ncol(vars$sta1d)-1)]
vars$sta1d[,ncol(vars$sta1d)] <- NA

vars$gdpgr <- vars$gdp
vars$gdpgr[,3:(ncol(vars$gdpgr)-1)] <- (vars$gdpgr[,4:(ncol(vars$gdpgr))] - vars$gdpgr[,3:(ncol(vars$gdpgr)-1)]) /  vars$gdpgr[,3:(ncol(vars$gdpgr)-1)]
vars$gdpgr[,ncol(vars$gdpgr)] <- NA

#define filter
# the filter defines the things we keep
filter <- complete.cases(vars$sta)
for(var in names(filters)){
	filter <- filter & rowSums(abs(vars[[var]][,-c(1,2)]) > filters[[var]],na.rm = T) == 0 
}
# apply filters
numPolID <- length(unique(vars$sta[,'polID']))
polIDsToRemove <- unique(vars$sta[!filter,'polID'])
filterInclSow <- vars$sta[,'polID'] %in% polIDsToRemove
for(var in names(vars)){
	vars[[var]] <- vars[[var]][filterInclSow,]
}
gc(verbose=F)
cat(sprintf('Filters applied: %i removed , %i remain\n',length(polIDsToRemove),numPolID-length(polIDsToRemove)))

# selecting run
# highest GDP in 2150, that stays below 2deg, in median
# idxOfSubSel <- which.max(vars$gdp[['2150']][vars$sta[['2150']]<2 & vars$gdp$sowID==5])
# selRun <- rownames(vars$gdp[vars$sta[['2150']]<2 & vars$gdp$sowID==5,])[idxOfSubSel]

# highest GDP in 2150, that stays below 2deg and has no more than 10 years in recession, in median
idxOfSubSel <- which.max(vars$gdp[['2150']][vars$sta[['2150']]<2 & vars$rdu[['2150']]<10 & vars$gdp$sowID==5])
selRun <- rownames(vars$gdp[vars$sta[['2150']]<2 & vars$rdu[['2150']]<10 & vars$gdp$sowID==5,])[idxOfSubSel]

# highest real gdp in 2075
# idxOfSubSel <- which.max(vars$gdp[['2075']][vars$gdp$sowID==5])
# selRun <- rownames(vars$gdp[vars$gdp$sowID==5,])[idxOfSubSel]

selPolID <- vars$gdp[selRun,1]

# get information on the selPolID
pdp.lst <- readRDS(file.path(location.output,'pdp.lst.RDS'))
pdpMeta <- readRDS(file.path(location.output,'pdpMeta.RDS'))
jointPolicies <- readRDS(file.path(location.output,'jointPolicies.RDS'))
samplePoints <- readRDS(file.path(location.output,'samplePoints.RDS'))

selPolDescStrs <- c()
i <- 0
for(domID in colnames(samplePoints)){
	if(!is.na(samplePoints[selPolID,domID])){
		i <- i+1
		pdpName <- pdpMeta$domain[pdpMeta$domID==domID][1]
		sdmID <- jointPolicies$sdmID[jointPolicies$domID==domID & jointPolicies$dplID==samplePoints[selPolID,domID]]
		if(!is.na(pdpMeta$subdomain[pdpMeta$domID==domID & pdpMeta$sdmID==sdmID])){
			pdpName <- paste0(pdpName,'+',pdpMeta$subdomain[pdpMeta$domID==domID & pdpMeta$sdmID==sdmID])
		} 
		sdpID <- jointPolicies$sdpID[jointPolicies$domID==domID & jointPolicies$dplID==samplePoints[selPolID,domID]]
		selPolDescStrs <-c(selPolDescStrs, pdp.lst[[pdpName]][pdp.lst[[pdpName]]$polID==sdpID,2])
	}
}
write(selPolDescStrs,file.path(figuresFolder,'selRunPolDesc.csv'))

# plot ####
plotLims <- list()
plotLims$sta <- c(0,7)
plotLims$gdp <- c(0,2e6)
plotLims$pop <- c(0,1.5e4)
plotLims$inf <- c(-0.4,0.4)
plotLims$sta1d <- c(-0.4,0.4)
plotLims$gdpgr <- c(-0.4,0.4)
plotLims$rdu <- c(0,100)
varUnits <- list()
varUnits$sta <- 'degC'
varUnits$gdp <- 'billion constant 2021 $'
varUnits$pop <- 'million people'
varUnits$sta1d <- 'degC'
varUnits$inf <- 'rate'
varUnits$gdpgr <- 'rate'
varUnits$rdu <- 'years'
varNames <- list()
varNames$sta <- 'Surface Temperature Anomaly'
varNames$gdp <- 'Real GDP'
varNames$pop <- 'Population'
varNames$inf <- 'Inflation'
varNames$sta1d <- 'Change in Surface Temperature Anomaly'
varNames$gdpgr <- 'Real GDP Growth Rate'
varNames$rdu <- 'Duration spent in Recession'
years <- seq(2150,1980,-2)
for(var in names(vars)){
	cat(sprintf('plotting %s ...',var))
	png(file.path(figuresFolder,paste0(var,'.png')),width = plt.w,height = plt.h,units = plt.u,res = plt.r)
	plot(0,0,type='n',xlim=range(years),ylim=plotLims[[var]],
			 xlab='year',
			 main=varNames[[var]],ylab=varUnits[[var]])
	grid(col='gray',lwd=2)
	for(year in years){
		boxplot(vars[[var]][[as.character(year)]][vars[[var]]$sowID==5],at=year,width = 0.8,add = T,pch='.',lwd=1,
						axes=F)
	}
	lines(as.numeric(colnames(vars[[var]][,-c(1,2)])),
				vars[[var]][selRun,-c(1,2)])
	dev.off()
	cat('\n')
}


