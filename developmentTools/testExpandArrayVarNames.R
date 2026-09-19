# testExpandArrayVarNames.R ####
#
# An array exported whole ([*,*]) has one per var file per element, and the
# plots have to find each of them from the one export entry. Everything else
# has to come through as it is.
#
# Run from the repository root:
#   Rscript developmentTools/testExpandArrayVarNames.R

source('funRunFRIDA.R')

ok <- 0
fail <- 0
check <- function(label,cond,detail=''){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s%s\n',label,ifelse(nchar(detail)>0,paste0('\n       ',detail),'')))
		fail <<- fail+1}
}

orig <- c('Coastal Assets.Coastal Assets[*,*]',
					'Sea level rise costs and impacts.Total fatalities due to coastal floods[*]',
					'GDP.Real GDP in 2021c$[1]',
					'Coastal Assets.Coastal Assets growth',
					'Energy.Output (fossil+renewable)/capita[*,*]',
					'Mix.Stock[*,*,Coal]',
					'Mix.Stock[*,*,Gas]',
					'Cube.Stock[*,*,*]',
					'Ages.Population[*,*]',
					'Absent.Nowhere[*,*]')
varNames <- cleanNames(orig)
# the columns of a run as Stella writes them, and cleaned as the per var files
# are named
headers <- c('Coastal Assets.Coastal Assets[1, insufficient]',
						 'Coastal Assets.Coastal Assets[1, well_protected]',
						 'Sea level rise costs and impacts.Total fatalities due to coastal floods[1]',
						 'GDP.Real GDP in 2021c$[1]',
						 'Coastal Assets.Coastal Assets growth',
						 'Energy.Output (fossil+renewable)/capita[1, oil]',
						 'Mix.Stock[1, insufficient, Coal]',
						 'Mix.Stock[1, well_protected, Coal]',
						 'Mix.Stock[1, insufficient, Gas]',
						 'Mix.Stock[1, well_protected, Gas]',
						 'Cube.Stock[1, insufficient, Coal]',
						 'Cube.Stock[1, insufficient, Gas]',
						 'Cube.Stock[1, well_protected, Coal]',
						 'Cube.Stock[1, well_protected, Gas]',
						 paste0('Ages.Population[1, ',c(1,2,10,12),']'))
available <- cleanNames(headers)

res <- funExpandArrayVarNames(varNames,orig,available,headers)
rowsOf <- function(entryOrig){res[res$entry==cleanNames(entryOrig),]}

cat('an array exported whole\n')
coastal <- rowsOf(orig[1])
check('becomes one row per element',
			setequal(coastal$name,c('coastal_assets_coastal_assets_insufficient',
															'coastal_assets_coastal_assets_well_protected')),
			paste(coastal$name,collapse=', '))
check('titled as Stella names the element, with the run as *',
			setequal(coastal$orig,c('Coastal Assets.Coastal Assets[*, insufficient]',
															'Coastal Assets.Coastal Assets[*, well_protected]')),
			paste(coastal$orig,collapse=', '))
check('without taking a variable of its own that shares the prefix',
			!'coastal_assets_coastal_assets_growth'%in%coastal$name)
energy <- rowsOf(orig[5])
check('with brackets, + and / in its name',
			identical(energy$name,available[6]) &&
				identical(energy$orig,'Energy.Output (fossil+renewable)/capita[*, oil]'),
			paste(energy$name,energy$orig))

cat('more dimensions\n')
coal <- rowsOf(orig[6])
gas <- rowsOf(orig[7])
check('an element fixed in the entry keeps the other entry\'s elements out',
			setequal(coal$name,available[7:8]) && setequal(gas$name,available[9:10]),
			paste(c(coal$name,gas$name),collapse=', '))
cube <- rowsOf(orig[8])
check('three wildcard dimensions give every combination',
			setequal(cube$name,available[11:14]),paste(cube$name,collapse=', '))
check('titled with every element',
			'Cube.Stock[*, well_protected, Gas]'%in%cube$orig,paste(cube$orig,collapse=', '))
ages <- rowsOf(orig[9])
check('numbered elements stay apart',
			nrow(ages)==4 && !any(duplicated(ages$name)),paste(ages$name,collapse=', '))

cat('everything else\n')
for(i in c(2,3,4)){
	r <- rowsOf(orig[i])
	check(sprintf('%s is kept as it is',orig[i]),
				nrow(r)==1 && r$name==varNames[i] && r$orig==orig[i],
				paste(r$name,r$orig))
}
absent <- rowsOf(orig[10])
check('an array with no element in the run is kept, so it is reported missing',
			nrow(absent)==1 && absent$name==varNames[10],paste(absent$name))
check('in the order of the entries',
			identical(unique(res$entry),varNames))

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('arrays exported whole are not expanded to their elements\n')
}
