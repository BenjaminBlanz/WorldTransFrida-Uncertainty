#!/bin/bash

if [ "$#" -eq 0 ]; then
  echo "Error: No input provided! Need to give expID of crashed run" >&2
  exit 1
fi

expID=$1

# This will only remove them, if they are links and not directories
rm FRIDAforUncertaintyAnalysis Stella_Simulator_Linux workerDirs

mv FRIDAforUncertaintyAnalysis-store FRIDAforUncertaintyAnalysis
mv Stella_Simulator_Linux-store Stella_Simulator_Linux

runScripts=workOutput/${expID}/runScriptsAndConfiguration
if [ ! -d ${runScripts} ]; then
	# output folders without runScriptsAndConfiguration keep the scripts at the top level
	runScripts=workOutput/${expID}
fi
cp ${runScripts}/*.R .
# a job that did not reach its cleanup leaves the .run in place
[ -f ${expID}.run ] || cp ${runScripts}/${expID}.run .
sbatch ${expID}.run
