# testCleanNames.R ####
#
# cleanNames names the per var files, the plot files and the columns of the run
# data. It must give every name the same result it always has, except where the
# old version dropped a _1 other than the run index, renamed a time other than
# Stella's Time column, or kept characters that do not belong in a file name.
#
# Run from the repository root:
#   Rscript developmentTools/testCleanNames.R

source('funRunFRIDA.R')

# the version before the rewrite, the reference for everything that must not change
cleanNamesBefore <- function(colNames){
	gsub('time','year',
			 gsub('_+$','',
			 		 gsub('_+','_',
			 		 		 gsub(',','_',
			 		 		 		 gsub('\\$','',
			 		 		 		 		 gsub('_1','',
			 		 		 		 		 		 gsub('\\]','_',
			 		 		 		 		 		 		 gsub('\\[\\*','_',
			 		 		 		 		 		 		 		 gsub('\\[\\d+','_',
			 		 		 		 		 		 		 		 		 gsub('[. ]','_',
			 		 		 		 		 		 		 		 		 		 tolower(colNames)))))))))))
}

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),'')))
		fail <<- fail+1}
}

exportList <- read.csv(file.path('FRIDA-info','frida_extra_variables_to_export_list.csv'))$FRIDA.FQN
exportList <- unique(exportList[nchar(exportList)>4])
# the columns Stella writes for them, the run index filled in. An array exported
# whole comes out per element, see testExpandArrayVarNames.R
headers <- sub('\\[\\*','[1',exportList[!grepl('\\[\\*,',exportList)])
# the names the old version got wrong: a 1 after a separator, time, or a
# character that is neither letter, digit nor part of the array notation
affected <- function(x){
	grepl('[ ._,[]1',x) | grepl('time',tolower(x)) | grepl('[^A-Za-z0-9 ._,*$[]',sub('\\]$','',x))
}

cat('names that were right stay as they were\n')
same <- !affected(exportList)
check(sprintf('all %i export entries without a known problem',sum(same)),
			identical(cleanNames(exportList[same]),cleanNamesBefore(exportList[same])),
			paste(exportList[same][cleanNames(exportList[same])!=cleanNamesBefore(exportList[same])],collapse='; '))
sameHeaders <- headers[!affected(headers)]
check('and their Stella columns',
			identical(cleanNames(sameHeaders),cleanNamesBefore(sameHeaders)))
check('read unmangled they get what they got read through read.csv\'s make.names',
			identical(cleanNames(sameHeaders),cleanNamesBefore(make.names(sameHeaders))))

cat('the file name and the column name agree\n')
check('the merge strips the run index before cleaning, the columns do not',
			identical(cleanNames(gsub('\\[\\d+','',headers)),cleanNames(headers)))

cat('the run index and Time\n')
check('X[1], X[*] and X] are all x',
			all(cleanNames(c('M.X[1]','M.X[*]','M.X]'))=='m_x'))
check('Time is year',cleanNames('Time')=='year')
check('a $ at the end is dropped',cleanNames('GDP.Real GDP in 2021c$[1]')=='gdp_real_gdp_in_2021c')

cat('what the old version got wrong\n')
fixed <- c('Demographics.Aged 1 to 20 Years[1]'='demographics_aged_1_to_20_years',
					 'Finance.Regulatory Tier 1 capital to assets[1]'='finance_regulatory_tier_1_capital_to_assets',
					 'GDP.future time in recession[1]'='gdp_future_time_in_recession',
					 'GDP.current recession time counter[1]'='gdp_current_recession_time_counter',
					 'Emissions.MtN2O emission per PCal animal products 1980[1]'='emissions_mtn2o_emission_per_pcal_animal_products_1980',
					 'policy_100DollarCarbonTax'='policy_100dollarcarbontax')
for(n in names(fixed)){
	check(sprintf('%s is %s',n,fixed[[n]]),cleanNames(n)==fixed[[n]],cleanNames(n))
}
numbered <- cleanNames(paste0('Ages.Population[1, ',c(1,2,10,11,12,100),']'))
check('numbered elements stay apart',!any(duplicated(numbered)),paste(numbered,collapse=', '))
check('and none of them looks like the whole variable',!'ages_population'%in%numbered)
odd <- cleanNames('Energy.Output (fossil+renewable)/capita[1]')
check('characters that do not belong in a file name become _',
			odd=='energy_output_fossil_renewable_capita',odd)
check('the * of a further wildcard dimension is kept',
			cleanNames('Coastal Assets.Coastal Assets[*,*]')=='coastal_assets_coastal_assets_*')

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('cleanNames does not give the names it should\n')
}
