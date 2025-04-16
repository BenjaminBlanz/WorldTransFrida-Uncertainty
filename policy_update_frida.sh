#!/bin/bash

RED='\033[0;31m'
NOCOLOR='\033[0m'

function cleanup_policy_repo() {
	cd ./FRIDAforPolicyAnalysisGit

	#remove all the stuff we don't need
	rm -rf "./Graphics"
	rm -f LICENSE
	rm -f README.md

	git apply --verbose $PWD/../FRIDAforPolicyAnalysis.patch
	cd ..

	rm -rf ./FRIDAforPolicyAnalysis
	rsync -av --exclude=".*" "./FRIDAforPolicyAnalysisGit/FRIDA Uncertainity Dashboard/" ./FRIDAforPolicyAnalysis

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


