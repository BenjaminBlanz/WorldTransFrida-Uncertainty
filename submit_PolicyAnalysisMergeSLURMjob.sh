#!/bin/bash

##############################################################################
########                  Change config here                          ########
##############################################################################
# Change properties of the SLURM job here. Change run specificaiton in 
# configPolicyAnalysis.R

### SLURM settings
# How long will it take approximately? Job will be killed after this time!
# But max 8:00 hours and the shorter the run is, the earlier the job gets run
hours=8
minutes=00

# Does the job need larger memory?
# (the larger the memory you request, the longer the job might sit in the queue, if the machine is full)
memorySize='256G' # can be ['256G' '512G', '1024G'], '256G' sometimes fails with 100,000 samples

# Use a different group account for ressources? Which partition?
account=uc1275
partition=compute

# Enter Email here in case you want to receive a mail, when the job failed
email=benjamin.blanz@uni-hamburg.de

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
      partition="$2"
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

# Create the runscript from the template
template='PolicyAnalysisMerge.run'
runscript="${locationOutput}/runFileForMergeSLURM.run"
cp $template $runscript

# Modify the template
sed -i "s/time=01:00:00/time=0${hours}:${minutes}:00/" $runscript
sed -i "s/account=REPLACEME/account=${account}/" $runscript
sed -i "s/partition=compute/partition=${partition}/" $runscript
sed -i "s/ntasks-per-node=128/ntasks-per-node=${numWorkers}/" $runscript
sed -i "s/constraint=256G/constraint=${memorySize}/" $runscript
# locationOutput contains slashes so use @ as the field limit symbol for sed
sed -i "s/max-connections=1024/max-connections=$(($(($numWorkers+5)) > 128 ? $(($numWorkers+5)) : 128))/g" $runscript
sed -i "s/runPolicyAnalysisMerger.R 128/runPolicyAnalysisMerger.R ${numWorkers}/g" $runscript
sed -i "s/expID/${expID}/g" $runscript
sed -i "s/jdoe@mail.com/${email}/" $runscript

#############################################################################
#######                    Submit the Job                            ########
#############################################################################
echo submitted > $locationOutput/mergeStatus.txt
#echo $runscript
if ! sbatch $runscript; then
	echo 'sbatch failed' > $locationOutput/mergeStatus.txt
fi

