#
# Prepares a csv of sample points for a run, called by
# submit_UncertaintyAnalysisLevante.sh --samplePointsFile.
#
#   Rscript runPrepareSamplePoints.R <csv> <parameter info> <out dir>
#
# writes to <out dir>
#   suppliedSamplePoints.RDS  the points, runs numbered by csv row
#   runCases.csv.gz           id,Case of every run, if the csv has a Case column
#   samplePointsParms.csv     frida_info of the csv's parameters: the model value
#                             from <parameter info>, e.g. Parameter Info.csv of
#                             the model, and the sampled range
# and prints numSample=<number of runs> as its last line.
#

args <- commandArgs(TRUE)
if(length(args)!=3){
	stop('usage: Rscript runPrepareSamplePoints.R <csv> <parameter info> <out dir>\n')
}
suppressMessages(library(data.table))
source('funRunFRIDA.R')
source('funParmSpace.R')

cat(sprintf('Reading sample points from %s...',args[1]))
sp <- funReadSamplePointsFile(args[1])
cat(sprintf('%i runs of %i parameters\n',nrow(sp$samplePoints),ncol(sp$samplePoints)))
parameterInfo <- funReadParameterFile(args[2],check.names=FALSE)
parms <- funSamplePointsParms(sp$samplePoints,parameterInfo)

dir.create(args[3],showWarnings=FALSE,recursive=TRUE)
saveRDS(sp$samplePoints,file.path(args[3],'suppliedSamplePoints.RDS'))
write.csv(parms,file.path(args[3],'samplePointsParms.csv'),row.names=FALSE)
if(!is.null(sp$cases)){
	fwrite(data.table(id=seq_along(sp$cases),Case=sp$cases),file.path(args[3],'runCases.csv.gz'))
}
cat(sprintf('numSample=%i\n',nrow(sp$samplePoints)))
