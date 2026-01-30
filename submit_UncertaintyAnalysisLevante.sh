#!/bin/bash

##############################################################################
########                  Change input here                           ########
##############################################################################
# If you want to change other things, these need to be changed in the
# respective .R files directly! Or, ideally, add the option to change them also
# here in the script. Any change to e.g. config.R or runMLEandParmSpace.R will
# always be used in the future, as these are the templates for this script!

numWorkers=250
numSample=10000
chunkSizePerWorker=100
likeCutoffRatio=1000

#### experiment id custom addition
# experiment id at bottom of input section to include policy name by default
expIDPreString='UA'

### SLURM settings
# How long will it take approximately? Job will be killed after this time!
# But max 8:00 hours and the shorter the run is, the earlier the job gets run
hours=2
minutes=00

# Does the job need larger memory?
# (the larger the memory you request, the longer the job might sit in the queue, if the machine is full)
memorySize='256G' # can be ['256G' '512G', '1024G'], '256G' sometimes fails with 100,000 samples


# Use a different group account for ressources? Which partition?
account=REPLACE_ME 
partition=compute

# Enter Email here in case you want to receive a mail, when the job failed
email=REPLACE_ME

# Copy ranges from pre-existing work?
# This should be the name of workOutput from which the files will be copied!
copyID='UA_EMB_nS100000'
# What to copy?
copyParmRangesAndScales="false"
copySamplePoints="false" # for run-by-run comparisons

# Which baseline configuration file to use as the template
configFile=config.R

# where to store outputs
baseOutputDir=workOutput

outputType='RDS' # can be 'both', 'csv' or 'RDS'
plotting='true' # avoiding the plotting can save quit some compute time
detRepSample='true' # needed for the plots that show the representative subsample
#and for generating the rep sample


#### OPTIONAL: USING DIFFERENT INPUT FILES?
# FRIDA folder
FRIDA='FRIDAforUncertaintyAnalysis'

# These need to be located in the FRIDA-configs/ folder!
policyFile='policy_EMB.csv'
climateFeedbackFile='ClimateFeedback_On.csv'
climateSTAOverrideFile='ClimateSTAOverride_Off.csv'
climateSTAOverrideFileTS='ClimateSTAOverrideTS.csv'
baselineParmFile=''

# These need to be located in the FRIDA-info/ folder!
infoFile='frida_info.csv'
externalRangesFile='frida_external_ranges.csv'
excludeParmFile='frida_parameter_exclusion_list.csv'
excludeVarFile='frida_variable_exclusion_list.csv'
extraExportFile='frida_extra_variables_to_export_list.csv'

##############################################################################
########                 End of input section                      ###########
##############################################################################

##############################################################################
########            Overrides from command line                    ###########
##############################################################################
while [ $# -gt 0 ]; do
  case "$1" in
    -s|--expIDPreString)
      expIDPreString="$2"
      ;;
    -w|--numWorkers)
      numWorkers="$2"
      ;;
    -n|--numSample)
      numSample="$2"
      ;;
    -k|--chunkSizePerWorker)
      chunkSizePerWorker="$2"
      ;;
    -h|--hours)
      hours="$2"
      ;;
    -m|--minutes)
      minutes="$2"
      ;;
    -M|--memorySize)
      memorySize="$2"
      ;;
    -a|--account)
      account="$2"
      ;;
    -p|--partition)
      partition="$2"
      ;;
    -e|--email)
      partition="$2"
      ;;
    -c|--cid|--copyID)
      copyID="$2"
      ;;
    --cpps|--copyParmRangesAndScales)
      copyParmRangesAndScales="$2"
      ;;
    --cpsp|--copySamplePoints)
      copySamplePoints="$2"
      ;;
    --frida|--FRIDA)
      FRIDA="$2"
      ;;
    --pol|--policyFile)
      policyFile="$2"
      ;;
    --cfb|--climateFeedbackFile)
      climateFeedbackFile="$2"
      ;;
    --sta|--climateSTAOverrideFile)
      climateSTAOverrideFile="$2"
      ;;
    --stats|--climateSTAOverrideFileTS)
      climateSTAOverrideFileTS="$2"
      ;;
    --blp|--baselineParmFile)
      baselineParmFile="$2"
      ;;
    --plt|--plotting)
    	plotting="$2"
    	;;
    --drs|--detRepSample)
    	detRepSample="$2"
    	;;
  	--lcr|--likeCutoffRatio)
  		likeCutoffRatio="$2"
  		;;
    --configFile)
      configFile="$2"
      ;;
    --baseOutputDir)
      baseOutputDir="$2"
      ;;
    --outputType)
      outputType="$2"
      ;;
    --extraExportFile)
    	extraExportFile="$2"
    	;;
    *)
      printf "Error: Invalid argument: $1 \n"
      exit 1
  esac
  shift
  shift
done

expID=${expIDPreString}-S${numSample}-${policyFile%.*}-${climateFeedbackFile%.*}-${climateSTAOverrideFile%.*}
if [ -n "${baselineParmFile%.*}" ]; then
	expID=${expID}-${baselineParmFile%.*}
fi
#############################################################################
########     Preparing the R-scripts and the runscript             ##########
#############################################################################

########
# Modify the config file according to the settings above
config=${expID}_config.R
cp $configFile $config

sed -i "s/^numWorkers <-.*$/numWorkers <- ${numWorkers}/" $config
sed -i "s/^numSample <-.*$/numSample <- ${numSample}/" $config
sed -i "s/^likeCutoffRatio <-.*$/likeCutoffRatio <- ${likeCutoffRatio}/" $config
sed -i "s/^chunkSizePerWorker <-.*$/chunkSizePerWorker <- ${chunkSizePerWorker}/" $config
sed -i "s/^name.output <-.*$/name.output <- '${expID}'/" $config
sed -i "s/config.R/${config}/g" $config
sed -i "s/FRIDAforUncertaintyAnalysis/${FRIDA}/" $config
sed -i "s/^name.frida_info <-.*$/name.frida_info <- '${infoFile}'/" $config
sed -i "s/^policyFileName <-.*$/policyFileName <- '${policyFile}'/" $config
sed -i "s/^climateFeedbackSpecFile <-.*$/climateFeedbackSpecFile <- '${climateFeedbackFile}'/" $config
sed -i "s/^climateOverrideSpecFile <-.*$/climateOverrideSpecFile <- '${climateSTAOverrideFile}'/" $config
sed -i "s/^climateOverrideSpecFileTS <- .*$/climateOverrideSpecFileTS <- '${climateSTAOverrideFileTS}'/" $config
sed -i "s/^name.baselineParmFile <-.*$/name.baselineParmFile <- '${baselineParmFile}'/" $config
sed -i "s/^name.frida_external_ranges <-.*$/name.frida_external_ranges <- '${externalRangesFile}'/" $config
sed -i "s/^name.frida_parameter_exclusion_list <-.*$/name.frida_parameter_exclusion_list <- '${excludeParmFile}'/" $config
sed -i "s/^name.frida_variable_exclusion_list <-.*$/name.frida_variable_exclusion_list <- '${excludeVarFile}'/" $config
sed -i "s/^name.frida_extra_variables_to_export_list <-.*$/name.frida_extra_variables_to_export_list <- '${extraExportFile}'/" $config
sed -i "s@'workOutput'@'${baseOutputDir}'@" $config # use @ as field boundaries as the path cotains /

if [ "$outputType" = "csv" ]; then
	sed -i "/^perVarOutputTypes/c\perVarOutputTypes <- c('csv')" $config
elif [ "$outputType" = "RDS" ]; then
	sed -i "/^perVarOutputTypes/c\perVarOutputTypes <- c('RDS')" $config
elif [ "$outputType" = "both" ]; then
	sed -i "/^perVarOutputTypes/c\perVarOutputTypes <- c('RDS','csv')" $config		
else
	echo "Output Type ${OutputType} not defined, exiting"
	exit 1
fi



#########
# Need to also modify the runMLE and runInitData, to avoid the breaks in R 
# that wait for someone to hit enter and to use the updated config

# modify runInitialise
runInit=${expID}_runInitialiseData.R
cp runInitialiseData.R $runInit

sed -i "/^continue <- readline/d" $runInit

# modify clusterHelp to use correct config
clusterHelp=${expID}_clusterHelp.R
cp clusterHelp.R $clusterHelp

sed -i "s/config.R/${config}/" $clusterHelp

# modify runMLE
runMLE=${expID}_runMLE.R
cp runMLEandParmSpace.R $runMLE

sed -i "s/config.R/${config}/g" $runMLE
sed -i "s/runInitialiseData.R/${runInit}/g" $runMLE
sed -i "s/clusterHelp.R/${clusterHelp}/g" $runMLE
sed -i "/^continue <- readline/d" $runMLE

# modify runPlotAllRuns
if [ "${plotting}" = "true" ]; then
	runPlot=${expID}_runPlotAllRuns.R
	cp runPlotAllRuns.R $runPlot
	sed -i "s/config.R/${config}/g" $runPlot
fi

# modify runDetermineRepresentativeSamples
if [ "${detRepSample}" = "true" ]; then
	runDetRepSample=${expID}_runDetermineRepresentativeSamples.R
	cp runDetermineRepresentativeSamples.R $runDetRepSample
	sed -i "s/config.R/${config}/g" $runDetRepSample
	sed -i "s/runInitialiseData.R/${runInit}/g" $runDetRepSample
fi


# Create the runscript from the template
template='UncertaintyAnalysis.run'
runscript="${expID}.run"
#runscript="${expID}_numW=${numWorkers}_numS=${numSample}_maxC=${maxConnections}.run"

cp $template $runscript

# Modify the template
sed -i "s/time=01:00:00/time=0${hours}:${minutes}:00/" $runscript
sed -i "s/account=mh0033/account=${account}/" $runscript
sed -i "s/partition=compute/partition=${partition}/" $runscript
sed -i "s/ntasks-per-node=128/ntasks-per-node=${numWorkers}/" $runscript
sed -i "s/constraint=256G/constraint=${memorySize}/" $runscript
sed -i "s/expID/${expID}/g" $runscript
sed -i "s/jdoe@mail.com/${email}/" $runscript
sed -i "s/config.R/${config}/" $runscript
sed -i "s/runMLEandParmSpace.R/${runMLE}/" $runscript
if [ "${plotting}" = "true" ]; then
	 sed -i "s/runPlotAllRuns.R/${runPlot}/" $runscript
else
	sed -i "s/source('runPlotAllRuns.R'); //" $runscript
fi
if [ "${detRepSample}" = "true" ]; then
	 sed -i "s/runDetermineRepresentativeSamples.R/${runDetRepSample}/" $runscript
else
	sed -i "s/source('runDetermineRepresentativeSamples.R'); //" $runscript
fi

#############################################################################
############ Copying parameter ranges and scales ############################
#############################################################################

expDir="workOutput/${expID}"
if [ -d "$expDir" ]; then
	echo $expDir already exists       
else
	mkdir $expDir
fi


if [ "$copyParmRangesAndScales" = "true" ]; then
	for filename in 'parscale.RDS' 'sampleParmsParscale.csv' 'sampleParmsParscaleRanged.csv' 'sampleParmsParscaleRanged.RDS' 'sigma-indepParms.RDS' 'sigma.RDS'; do
		file=workOutput/${copyID}/${filename}
		if [ -e $file ]; then
			cp ${file}  ${expDir}/
		else
			echo "File ${file} not found, exiting"
			exit 1
		fi
	done
fi


if [ "$copySamplePoints" = "true" ]; then
	for filename in 'samplePoints.RDS'; do
		file=workOutput/${copyID}/${filename}
		if [ -e $file ]; then
			cp ${file}  ${expDir}/
		else
			echo "File ${file} not found, exiting"
			exit 1
		fi
	done
fi

#############################################################################
#######                    Submit the Job                            ########
#############################################################################
echo submitted > $expDir/status
echo $runscript
sbatch $runscript

