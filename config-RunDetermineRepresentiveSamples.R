# for DPIS rep sample
preexsistingBaselineFolder <- NA
numSample <- 2e4
specFileForScenaro <- 'policy_EMB.csv'
# spec files for policies to be evaluated in the chosen scenario
specFilesForPols <- c('policy_300DollarCarbonTax.csv','policy_100DollarCarbonTax.csv')

discountFactor <- 0.05

# representative subsample ####
subSample.NumSamplePerVar <- 11
subSample.Ps <- seq(0.5/subSample.NumSamplePerVar,1-0.5/subSample.NumSamplePerVar,
										length.out=subSample.NumSamplePerVar)
subSample.TargetVars <- c('demographics_real_gdp_per_person')
