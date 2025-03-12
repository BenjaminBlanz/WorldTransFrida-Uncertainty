#!/bin/bash

expID=$1

# This will only remove them, if they are links and not directories
rm FRIDAforUncertaintyAnalysis Stella_Simulator_Linux workersDir

mv FRIDAforUncertaintyAnalysis-store FRIDAforUncertaintyAnalysis
mv Stella_Simulator_Linux-store Stella_Simulator_Linux

cp workOutput/${expID}/*.Run .
sbatch ${expID}.run
