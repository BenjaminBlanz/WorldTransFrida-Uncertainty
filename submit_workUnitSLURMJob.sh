#!/bin/bash

##############################################################################
########                  Change config here                          ########
##############################################################################

# These settings are only for running this file manually,
# if using SLURMrunner, these values are taken from 
# respective config.

# How long will it take approximately? Job will be killed after this time!
# But max 8:00 hours and the shorter the run is, the earlier the job gets run
hours=8
minutes=00

# Does the job need larger memory?
# (the larger the memory you request, the longer the job might sit in the queue, if the machine is full)
memorySize='256G' # can be ['256G' '512G', '1024G'], '256G' sometimes fails with 100,000 samples

# Use a different group account for ressources? Which partition?
account=REPLACEME 
partition=compute

# Enter Email here in case you want to receive a mail, when the job failed
email=na@na.na

# where to store outputs
baseOutputDir=workOutput
outputType='both' # can be 'both', 'csv' or 'RDS'
plotting='true' # avoiding the plotting can save quit some compute time
detRepSample='true' # needed for the plots that show the representative subsample
#and for generating the rep sample

numWorkers=250
numSample=10000
chunkSizePerWorker=100
likeCutoffRatio=1000

#### OPTIONAL: USING DIFFERENT INPUT FILES?
# baseline configuration
configFile=config.R
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
    -w|--numWorkers)
      numWorkers="$2"
      ;;
    -s|--expID)
      expID="$2"
      ;;
    -u|--workUnitID)
      workUnitID="$2"
      ;;
    -o|--locationOutput)
      locationOutput="$2"
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
      email="$2"
      ;;
    --configFile)
      configFile="$2"
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
    --baseOutputDir)
      baseOutputDir="$2"
      ;;
    *)
      printf "Error: Invalid argument: $1 \n"
      exit 1
  esac
  shift
  shift
done

#############################################################################
########     Preparing the R-scripts and the runscript             ##########
#############################################################################

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
sed -i "s/'workOutput'/'${baseOutputDir}'/" $config

# Create the runscript from the template
template='workUnit.run'
runscript="${locationOutput}/workUnits/workUnit-${workUnitID}/runFileForSLURM.run"
cp $template $runscript

# Modify the template
sed -i "s/time=01:00:00/time=0${hours}:${minutes}:00/" $runscript
sed -i "s/account=REPLACEME/account=${account}/" $runscript
sed -i "s/partition=compute/partition=${partition}/" $runscript
sed -i "s/ntasks-per-node=128/ntasks-per-node=${numWorkers}/" $runscript
sed -i "s/constraint=256G/constraint=${memorySize}/" $runscript
# locationOutput contains slashes so use @ as the field limit symbol for sed
sed -i "s@locationOutput@${locationOutput}@g" $runscript
sed -i "s/max-connections=1024/max-connections=$(($(($numWorkers+5)) > 128 ? $(($numWorkers+5)) : 128))/g" $runscript
sed -i "s/workUnitID/${workUnitID}/g" $runscript
sed -i "s/expID/${expID}/g" $runscript
sed -i "s/jdoe@mail.com/${email}/" $runscript
sed -i "s/config.R/${configFile}/" $runscript

#############################################################################
#######                    Submit the Job                            ########
#############################################################################
echo submitted > $locationOutput/workUnits/workUnit-${workUnitID}/status.txt
#echo $runscript
if ! sbatch $runscript; then
	# if the sbatch failed (e.g. bad username or whatever) set the status back to prepared
	echo prepared > $locationOutput/workUnits/workUnit-${workUnitID}/status.txt
fi

