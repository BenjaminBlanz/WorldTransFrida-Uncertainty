#!/bin/bash

RED='\033[0;31m'
NOCOLOR='\033[0m'

function cleanup_policy_repo() {
	cd ./FRIDAforPolicyAnalysisGit

	#remove all the stuff we don't need
	rm -rf "./Graphics"
	rm -f LICENSE
	rm -f README.md
	sed -i 's/<sim_specs isee:sim_duration="0" isee:run_prefix="Run" isee:simulation_delay="0" isee:restore_on_start="false" isee:save_interval="1" method="RK4" time_units="Year" isee:instantaneous_flows="false" isee:ignore_module_errors="false" isee:strict_units="false" isee:loop_scores="false" isee:loop_exhaustive_allowed="1000">/<sim_specs isee:sim_duration="0" isee:run_prefix="Run" isee:simulation_delay="0" isee:restore_on_start="false" isee:save_interval="1" method="RK4" time_units="Year" isee:instantaneous_flows="false" isee:ignore_module_errors="false" isee:strict_units="false" isee:loop_scores="false" isee:loop_exhaustive_allowed="1000" isee:analyze_mode="false">/' FRIDA.stmx
	sed -i 's/<import frequency="on_demand" isee:overwrite="true" resource="r..\/Data\/subSampleParameterValues.csv"\/>/<import frequency="on_demand" isee:overwrite="true" resource="r..\/Data\/uncertainty_parameters.csv"\/>/' FRIDA.stmx
	git apply --verbose $PWD/../FRIDAforPolicyAnalysis.patch
	cd ..

	rm -rf ./FRIDAforPolicyAnalysis
	rsync -av --exclude=".*" "./FRIDAforPolicyAnalysisGit/" ./FRIDAforPolicyAnalysis

	#bring in anything from the template
	cp -R ./FRIDAforPolicyAnalysis_Template/* ./FRIDAforPolicyAnalysis/
}


if [ ! -d ./FRIDAforPolicyAnalysisGit ] ; then
	#clone the FRIDA model so you always have the latest
	echo "Cloning FRIDA from the main metno/WorldTransFRIDA repo"
	git clone https://github.com/BenjaminBlanz/WorldTransFRIDA-SimpleDashboard.git FRIDAforPolicyAnalysisGit

	cleanup_policy_repo

else
	# reset your model to be the main model
	echo "Resetting your FRIDA to be the latest..."
	cd ./FRIDAforPolicyAnalysisGit
	git fetch origin main
	git reset --hard origin/main

	cd ..

	cleanup_policy_repo
fi

echo -e ${RED}"Make sure the FRIDA Uncertainity Dashboard has been updated with the latest reduced uncertainity ensemble"${NOCOLOR}
read -n 1 -s -r -p $'Press any key to continue...\n'


