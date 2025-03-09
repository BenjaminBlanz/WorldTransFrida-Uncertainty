#!/bin/bash

RED='\033[0;31m'
NOCOLOR='\033[0m'

function cleanup_repo() {
	cd ./FRIDAforUncertaintyAnalysisGit

	#remove all the stuff we don't need
	rm -rf "./Data Processing"
	rm -rf "./FRIDA Dashboard"
	rm -rf "./FRIDA ILE"
	rm -rf "./FRIDA Uncertainity Dashboard"
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

	git apply $PWD/../FRIDAforUncertaintyAnalysis.patch
	cd ..

	rm -rf ./FRIDAforUncertaintyAnalysis
	rsync -av --exclude=".*" ./FRIDAforUncertaintyAnalysisGit/ ./FRIDAforUncertaintyAnalysis

	#bring in anything from the template
	cp -R ./FRIDAforUncertaintyAnalysis_Template/* ./FRIDAforUncertaintyAnalysis/

}

if [ ! -d ./FRIDAforUncertaintyAnalysisGit ] ; then
	#clone the FRIDA model so you always have the latest
	echo "Cloning FRIDA from the main metno/WorldTransFRIDA repo"
	git clone git@github.com:metno/WorldTransFRIDA.git FRIDAforUncertaintyAnalysisGit

	cleanup_repo

else
	# reset your model to be the main model
	echo "Resetting your FRIDA to be the latest..."
	cd ./FRIDAforUncertaintyAnalysisGit
	git fetch origin main
	#git reset --hard origin/main
	
	cd ..

	cleanup_repo
fi



echo -e ${RED}"Make sure you have the latest ./FRIDA-Info/frida_info.csv before running an uncertainity analysis!"${NOCOLOR}
read -n 1 -s -r -p "Press any key to continue"


