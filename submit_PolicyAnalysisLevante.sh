#!/bin/bash

##############################################################################
########                  Change input here                           ########
##############################################################################
# If you want to change other things, these need to be changed in the
# respective .R files directly! Or, ideally, add the option to change them also
# here in the script. Any change to e.g. config.R or runMLEandParmSpace.R will
# always be used in the future, as these are the templates for this script!

# number of joint policy scenarios
numInitialJointPol=100000

# probability of status quo policy per policy domain
# in sampling joint policies
# This ensures that not every joint policy run includes 
# variation in all policies.
# Higher values cause more sparse joint policies
nullPolProb=0.5

chunkSizePerWorker=100
#### experiment id custom addition
# experiment id at bottom of input section to include policy name by default
expIDPreString='PA'

### SLURM settings
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

outputType='RDS' # can be 'both', 'csv' or 'RDS'

#### OPTIONAL: USING DIFFERENT INPUT FILES?
# FRIDA folder
FRIDA='FRIDAforPolicyAnalysis'

# These need to be located in the FRIDA-info/ folder!
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
    -n|--numJointPol)
      numInitialJointPol="$2"
      ;;
    -P|--nullPolProb)
      nullPolProb="$2"
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
    --frida|--FRIDA|--FRIDAlocation)
      FRIDA="$2"
      ;;
    *)
      printf "Error: Invalid argument: $1 \n"
      exit 1
  esac
  shift
  shift
done

expID=${expIDPreString}-S${numJointPol}-P${nullPolProb}

#############################################################################
########     Preparing the R-scripts and the runscript             ##########
#############################################################################

########
# Modify the config file according to the settings above
config=${expID}_config.R
cp configPolicyAnalysis.R $config

sed -i "s/numWorkers <- parallel::detectCores()/numWorkers <- ${numWorkers}/" $config
sed -i "s/numInitialJointPol <- 1e5/numInitialJointPol <- ${numInitialJointPol}/" $config
sed -i "s/nullPolProb <- 0.5/nullPolProb <- ${nullPolProb}/" $config
sed -i "s/chunkSizePerWorker <- 100/chunkSizePerWorker <- ${chunkSizePerWorker}/" $config
sed -i "s/name.output <- gsub('\\\\\\.','_',paste0('N-',numInitialJointPol,'-nPr-',nullPolProb))/name.output <- '${expID}'/" $config
sed -i "s/config.R/${config}/g" $config
sed -i "s/FRIDAforPolicyAnalysis/${FRIDA}/" $config
sed -i "s/frida_extra_variables_to_export_list.csv/${extraExportFile}/" $config

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

# modify clusterHelp to use correct config
clusterHelp=${expID}_clusterHelp.R
cp clusterHelp.R $clusterHelp

sed -i "s/config.R/${config}/" $clusterHelp

# modify runPolicyAnalysis
runPolicyAnalysis=${expID}_runPolicyAnalysis.R
cp runPolicyAnalysis.R $runPolicyAnalysis

sed -i "s/config.R/${config}/g" $runMLE
sed -i "s/clusterHelp.R/${clusterHelp}/g" $runMLE

# Create the runscript from the template
template='PolicyAnalysis.run'
runscript="${expID}.run"

cp $template $runscript

# Modify the template
sed -i "s/time=01:00:00/time=0${hours}:${minutes}:00/" $runscript
sed -i "s/account=REPLACACCNAME/account=${account}/" $runscript
sed -i "s/partition=compute/partition=${partition}/" $runscript
sed -i "s/ntasks-per-node=128/ntasks-per-node=${numWorkers}/" $runscript
sed -i "s/constraint=256G/constraint=${memorySize}/" $runscript
sed -i "s/expID/${expID}/g" $runscript
sed -i "s/jdoe@mail.com/${email}/" $runscript
sed -i "s/config.R/${config}/" $runscript
sed -i "s/runPolicyAnalysis.R/${runMLE}/" $runscript

#############################################################################
#######                    Submit the Job                            ########
#############################################################################
echo submitted > $expDir/status
echo $runscript
sbatch $runscript

