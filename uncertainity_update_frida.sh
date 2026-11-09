#!/bin/bash

RED='\033[0;31m'
NOCOLOR='\033[0m'

function fail() {
	echo -e "${RED}Error: $1, FRIDAforUncertaintyAnalysis is left as it was${NOCOLOR}"
	exit 1
}

if ! command -v git > /dev/null; then
	fail "git not found (on levante: source /sw/etc/profile.levante; module load git)"
fi

frida_source=$(grep name.frida_source config.R | sed 's/name.frida_source <- //' | sed "s/['\"]//g")
frida_branch=$(grep name.frida_branch config.R | sed 's/name.frida_branch <- //' | sed "s/['\"]//g")

function cleanup_uncertainty_repo() {
	cd ./FRIDAforUncertaintyAnalysisGit || fail "FRIDAforUncertaintyAnalysisGit not found"

	#remove all the stuff we don't need
	rm -rf "./Data Processing"
	rm -rf "./docs"
	rm -f FRIDA.isdb
	rm -f LICENSE
	rm -f README.md

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
	git apply "$PWD/../FRIDAforAnalysis.patch" || fail "FRIDAforAnalysis.patch does not apply"
	cd ..

	rm -rf ./FRIDAforUncertaintyAnalysis
	if ! rsync -av --exclude=".*" ./FRIDAforUncertaintyAnalysisGit/ ./FRIDAforUncertaintyAnalysis; then
		echo -e "${RED}Error: copying FRIDAforUncertaintyAnalysisGit to FRIDAforUncertaintyAnalysis failed${NOCOLOR}"
		exit 1
	fi

	#bring in anything from the template
	cp -R ./FRIDAforUncertaintyAnalysis_Template/* ./FRIDAforUncertaintyAnalysis/

}

if [ ! -d ./FRIDAforUncertaintyAnalysisGit ] ; then
	#clone the FRIDA model so you always have the latest
	echo "Cloning FRIDA from the ${frida_branch} metno/WorldTransFRIDA repo"
	git clone -b ${frida_branch} --depth 1 ${frida_source} FRIDAforUncertaintyAnalysisGit ||
		fail "cloning ${frida_branch} of ${frida_source} failed"

	cleanup_uncertainty_repo

else
	# reset your model to be the ${frida_branch} model
	echo "Resetting your FRIDA to be the latest..."
	cd ./FRIDAforUncertaintyAnalysisGit
	git reset --hard || fail "resetting FRIDAforUncertaintyAnalysisGit failed"
	# explicit refspec so this works even when the existing clone was made
	# single-branch (--depth 1 -b X) against a *different* branch
	git fetch --depth 1 origin "+refs/heads/${frida_branch}:refs/remotes/origin/${frida_branch}" ||
		fail "fetching ${frida_branch} failed"
	git checkout -B ${frida_branch} "origin/${frida_branch}" || fail "checking out ${frida_branch} failed"
	git reset --hard "origin/${frida_branch}" || fail "resetting to origin/${frida_branch} failed"

	cd ..

	cleanup_uncertainty_repo
fi



echo -e ${RED}"IGNORE if you're v3.1 or newer: Make sure you have the latest ./FRIDA-Info/frida_info.csv before running an uncertainity analysis!"${NOCOLOR}
if [ -t 0 ]; then
	read -n 1 -s -r -p $'Press any key to continue...\n'
fi


