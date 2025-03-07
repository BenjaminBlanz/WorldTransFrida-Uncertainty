#!/bin/bash

##############################################################################
########                  Change input here                           ########
##############################################################################
# If you want to change other things, these need to be changed in the
# respective .R files directly! Or, ideally, add the option to change them also
# here in the script. Any change to e.g. config.R or runMLEandParmSpace.R will
# always be used in the future, as these are the templates for this script!

expID=UA_nW120_nS10000
numWorkers=120
numSamples=1.0e4

### SLURM settings
# How long will it take approximately? Job will be killed after this time!
# But max 8:00 hours and the shorter the run is, the earlier the job gets run
hours=1
minutes=00

# Use a different group account for ressources? Which partition?
account=mh0033 
partition=compute

# Enter Email here in case you want to receive a mail, when the job failed
email=Testmailforscript@UncertaintyAnalysis.abc

# Copy ranges from pre-existing work?
copyParmRangesAndScales="true"
# This should be the name of workOutput from which the files will be copied!
copyID=UA_nW250_nS10000_maxC1024_no2



#### OPTIONAL: USING DIFFERENT INPUT FILES?
# FRIDA folder
FRIDA='FRIDAforSLRpaper'

# These need to be located in the FRIDA-configs/ folder!
policyFile='policy_EMB.csv'

# These need to be located in the FRIDA-info/ folder!
infoFile='frida_info.csv'
externalRangesFile='frida_external_ranges.csv'
excludeParmFile='frida_parameter_exclusion_list.csv'
excludeVarFile='frida_variable_exclusion_list.csv'
extraExportFile='frida_extra_variables_to_export_list.csv'
##############################################################################
########                 End of input section                      ###########
##############################################################################


#############################################################################
############ Preparing the R-scripts and the runscript ######################
#############################################################################

# Modify the config file according to the settings above
config=${expID}_config.R
cp TEMPLATE_config.R $config

sed -i "s/numWorkers <- 120/numWorkers <- ${numWorkers}/" $config
sed -i "s/numSample <- 1e4/numSample <- ${numSamples}/" $config
sed -i "s/file.path('workOutput',name.output)/file.path('workOutput','${expID}')/" $config
sed -i "s/config.R/${config}/g" $config
sed -i "s/FRIDAforUncertaintyAnalysis/${FRIDA}/" $config
sed -i "s/frida_info.csv/${infoFile}/" $config
sed -i "s/frida_external_ranges.csv/${externalRangesFile}/" $config
sed -i "s/frida_parameter_exclusion_list.csv/${excludeParmFile}/" $config
sed -i "s/frida_variable_exclusion_list.csv/${excludeVarFile}/" $config
sed -i "s/frida_extra_variables_to_export_list.csv/${extraExportFile}/" $config


# Need to also modify the runMLE and runInitData, to avoid the breaks in R 
# that wait for someone to hit enter and to use the updated config

# modify runInitialise
runInit=${expID}_runInitialiseData.R
cp TEMPLATE_runInitialiseData.R $runInit

sed -i "/^continue <- readline/d" $runInit

# modify clusterHelp to use correct config
clusterHelp=${expID}_clusterHelp.R
cp TEMPLATE_clusterHelp.R $clusterHelp

sed -i "s/config.R/${config}/" $clusterHelp

# modify runMLE
runMLE=${expID}_runMLE.R
cp TEMPLATE_runMLEandParmSpace.R $runMLE

sed -i "s/config.R/${config}/g" $runMLE
sed -i "s/runInitialiseData.R/${runInit}/g" $runMLE
sed -i "s/clusterHelp.R/${clusterHelp}/g" $runMLE
sed -i "/^continue <- readline/d" $runMLE

# Create the runscript from the template
template='UncertaintyAnalysis.run'
runscript="${expID}.run"
#runscript="${expID}_numW=${numWorkers}_numS=${numSamples}_maxC=${maxConnections}.run"

cp $template $runscript

# Modify the template
sed -i "s/time=01:00:00/time=0${hours}:${minutes}:00/" $runscript
sed -i "s/account=mh0033/account=${account}/" $runscript
sed -i "s/partition=compute/partition=${partition}/" $runscript
sed -i "s/LOG.UncertaintyAnalysis/LOG.${expID}/g" $runscript
sed -i "s/jdoe@mail.com/${email}/" $runscript
sed -i "s/job-name=runUncertaintyAnalysis/job-name=${expID}/" $runscript
sed -i "s/runMLEandParmSpace.R/${runMLE}/" $runscript

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
	for filename in 'parscale.RDS' 'sampleParmsParscale.csv' 'sampleParmsParscaleRanged.csv' 'sampleParmsParscaleRanged.RDS' 'samplePointsBase.RDS' 'samplePoints.RDS' 'sigma-indepParms.RDS' 'sigma.RDS'; do
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
sbatch $runscript
