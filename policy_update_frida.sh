#!/bin/bash

RED='\033[0;31m'
NOCOLOR='\033[0m'

function cleanup_policy_repo() {
	cd ./FRIDAforPolicyAnalysisGit

	#remove all the stuff we don't need
	rm -rf "./Data"
	rm -rf "./Data Processing"
	rm -rf "./FRIDA Dashboard"
	rm -rf "./FRIDA ILE"
	rm -rf "./FRIDA_Modules"
	rm -rf "./FRIDA Uncertainity Dashboard/Graphics"
	rm -r FRIDA.isdb
	rm -r FRIDA.sxtmx
	rm -f LICENSE
	rm -f ReadMe.md

	git apply $PWD/../FRIDAforAnalysis.patch
	cd ..

	rm -rf ./FRIDAforPolicyAnalysis
	rsync -av --exclude=".*" "./FRIDAforPolicyAnalysisGit/FRIDA Uncertainity Dashboard/" ./FRIDAforPolicyAnalysis

	#bring in anything from the template
	cp -R ./FRIDAforPolicyAnalysis_Template/* ./FRIDAforPolicyAnalysis/
}


if [ ! -d ./FRIDAforPolicyAnalysisGit ] ; then
	#clone the FRIDA model so you always have the latest
	echo "Cloning FRIDA from the main metno/WorldTransFRIDA repo"
	git clone git@github.com:metno/WorldTransFRIDA.git FRIDAforPolicyAnalysisGit

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


