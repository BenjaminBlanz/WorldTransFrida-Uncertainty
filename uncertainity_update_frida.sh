#!/bin/bash

RED='\033[0;31m'
NOCOLOR='\033[0m'

function cleanup_uncertainty_repo() {
	cd ./FRIDAforUncertaintyAnalysisGit

	#remove all the stuff we don't need
	rm -rf "./Data Processing"
	rm -r FRIDA.isdb
	rm -f LICENSE
	rm -f ReadMe.md

	#save what we need to save from the /Data folder
	mv "./Data/Calibration Data.csv" "Calibration Data.csv"
	mv "./Data/frida_input_data.csv" "frida_input_data.csv"

	#recreate the data folder
	rm -rf ./Data
	mkdir ./Data

	#put back the stuff we need
	mv "Calibration Data.csv" "./Data/Calibration Data.csv"
	mv "frida_input_data.csv" "./Data/frida_input_data.csv"
	sed -i 's/<sim_specs isee:sim_duration="0" isee:run_prefix="Run" isee:simulation_delay="0" isee:restore_on_start="false" isee:save_interval="1" method="RK4" time_units="Year" isee:instantaneous_flows="false" isee:ignore_module_errors="false" isee:strict_units="false" isee:loop_scores="false" isee:loop_exhaustive_allowed="1000">/<sim_specs isee:sim_duration="0" isee:run_prefix="Run" isee:simulation_delay="0" isee:restore_on_start="false" isee:save_interval="1" method="RK4" time_units="Year" isee:instantaneous_flows="false" isee:ignore_module_errors="false" isee:strict_units="false" isee:loop_scores="false" isee:loop_exhaustive_allowed="1000" isee:analyze_mode="false">/' FRIDA.stmx
	git apply $PWD/../FRIDAforAnalysis.patch
	cd ..

	rm -rf ./FRIDAforUncertaintyAnalysis
	rsync -av --exclude=".*" ./FRIDAforUncertaintyAnalysisGit/ ./FRIDAforUncertaintyAnalysis

	#bring in anything from the template
	cp -R ./FRIDAforUncertaintyAnalysis_Template/* ./FRIDAforUncertaintyAnalysis/

}

if [ ! -d ./FRIDAforUncertaintyAnalysisGit ] ; then
	#clone the FRIDA model so you always have the latest
	echo "Cloning FRIDA from the main metno/WorldTransFRIDA repo"
	git clone https://github.com/metno/WorldTransFRIDA.git FRIDAforUncertaintyAnalysisGit

	cleanup_uncertainty_repo

else
	# reset your model to be the main model
	echo "Resetting your FRIDA to be the latest..."
	cd ./FRIDAforUncertaintyAnalysisGit
	git fetch origin main
	git reset --hard v2.1

	cd ..

	cleanup_uncertainty_repo
fi



echo -e ${RED}"Make sure you have the latest ./FRIDA-Info/frida_info.csv before running an uncertainity analysis!"${NOCOLOR}
read -n 1 -s -r -p $'Press any key to continue...\n'


