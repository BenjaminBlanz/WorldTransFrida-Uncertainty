#!/bin/bash
# testSubmitSamplePoints.sh ####
#
# submit_UncertaintyAnalysisLevante.sh --samplePointsFile runs the points of a
# csv: the row count sets numSample, a frida_info file of its parameters goes to
# FRIDA-info/, the points and the Case labels to the output folder. Resubmitting
# the same command resumes the run, a different csv under the same name stops.
#
# Runs in a copy of the repository with a stand-in model, stella and sbatch, as
# the submit script writes its scripts into the repository and runs the model.
#
# Run from the repository root:
#   bash developmentTools/testSubmitSamplePoints.sh

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
repo="$tmp/repo"
mkdir -p "$repo" "$tmp/bin"
git ls-files -z --cached --others --exclude-standard | tar --null -T - -cf - | tar -C "$repo" -xf -

mkdir -p "$repo/fixtureModel"
echo 'a model' > "$repo/fixtureModel/FRIDA.stmx"
printf 'Variable,Value,Min,Max\n"p.a[1]",2,0,10\n"p.b[1]",5,0,10\n' > "$repo/fixtureModel/Parameter Info.csv"
mkdir -p "$repo/Stella_Simulator_Linux"
printf '#!/bin/sh\nexit 0\n' > "$repo/Stella_Simulator_Linux/stella_simulator"
chmod +x "$repo/Stella_Simulator_Linux/stella_simulator"
# the package check of the submit script needs every package of a run, which
# runPrepareSamplePoints.R does not
echo 'invisible(TRUE)' > "$repo/initialise.R"
printf '#!/bin/sh\necho Submitted batch job 1\n' > "$tmp/bin/sbatch"
chmod +x "$tmp/bin/sbatch"
export PATH="$tmp/bin:$PATH"

printf 'Case,p.a[1],p.b[1]\nA:1,1,2\nB:1,3,4\nA_B[p.a[1]]:1,3,2\n' > "$repo/pts.csv"
printf 'p.a[1],p.b[1]\n1,2\n3,4\n' > "$repo/ptsNoCase.csv"

ok=0
fail=0
check(){
	if eval "$2"; then echo "  ok   $1"; ok=$((ok+1)); else echo "  FAIL $1"; fail=$((fail+1)); fi
}
submit(){
	(cd "$repo" && ./submit_UncertaintyAnalysisLevante.sh --frida fixtureModel --plt false --drs false "$@") > "$tmp/log" 2>&1
	local status=$?
	{ echo "=== submit $*"; cat "$tmp/log"; } >> "$tmp/allLogs"
	return $status
}

expID=spTest-S3-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_Off
out="$repo/workOutput/$expID"

echo 'a csv with Case'
submit --samplePointsFile pts.csv -n 99 -s spTest
check 'the submit succeeds' "[ $? -eq 0 ]"
check 'numSample comes from the csv, not -n' "[ -d '$out' ] && grep -q 'numSample 3' '$tmp/log'"
check 'the points are in the output folder' "[ -f '$out/suppliedSamplePoints.RDS' ]"
check 'the Case labels too, one per run' \
	"[ \"\$(zcat '$out/runCases.csv.gz' | head -1)\" = 'id,Case' ] && [ \$(zcat '$out/runCases.csv.gz' | wc -l) -eq 4 ]"
check 'the frida_info of the points is in FRIDA-info' "[ -f '$repo/FRIDA-info/$expID-samplePointsParms.csv' ]"
config="$repo/${expID}_config.R"
check 'the config runs 3 samples' "grep -q '^numSample <- 3$' '$config'"
check 'from the csv' "grep -qF \"name.samplePointsFile <- 'pts.csv'\" '$config'"
check 'with its frida_info as info file' "grep -qF \"name.frida_info <- '$expID-samplePointsParms.csv'\" '$config'"
check 'and as external ranges' "grep -qF \"name.frida_external_ranges <- '$expID-samplePointsParms.csv'\" '$config'"
check 'nothing staged is left behind' "[ -z \"\$(ls -A '$repo/workOutput' | grep samplePoints)\" ]"

echo 'resubmitting'
submit --samplePointsFile pts.csv -s spTest
check 'the same csv again is accepted' "[ $? -eq 0 ]"
sed -i 's/^B:1,3,4$/B:1,3,5/' "$repo/pts.csv"
submit --samplePointsFile pts.csv -s spTest
check 'a changed csv under the same name stops' "[ $? -ne 0 ] && grep -q 'differ' '$tmp/log'"

echo 'a csv without Case'
submit --samplePointsFile ptsNoCase.csv -s spNoCase
check 'the submit succeeds' "[ $? -eq 0 ]"
outNoCase="$repo/workOutput/spNoCase-S2-policy_EMB-ClimateFeedback_On-ClimateSTAOverride_Off"
check 'the points are there' "[ -f '$outNoCase/suppliedSamplePoints.RDS' ]"
check 'without labels' "[ ! -e '$outNoCase/runCases.csv.gz' ]"

echo 'options the csv replaces'
submit --samplePointsFile ptsNoCase.csv -s spConflict --infoFile x.csv
check '--infoFile stops' "[ $? -ne 0 ] && grep -q 'samplePointsFile' '$tmp/log'"
submit --samplePointsFile ptsNoCase.csv -s spConflict --cpsp true
check '--cpsp true stops' "[ $? -ne 0 ] && grep -q 'samplePointsFile' '$tmp/log'"
grep -v '^name.samplePointsFile <-' "$repo/config.R" > "$repo/olderConfig.R"
submit --samplePointsFile ptsNoCase.csv -s spOlderConfig --conf olderConfig.R
check 'a config template without name.samplePointsFile stops' \
	"[ $? -ne 0 ] && grep -q 'name.samplePointsFile' '$tmp/log'"

echo
echo "$ok ok, $fail failed"
if [ $fail -gt 0 ]; then
	echo 'output of the submits, without the package loading:'
	grep -v -e '^Loading' -e '^Attaching' -e '^The following' -e '^    ' "$tmp/allLogs"
	exit 1
fi
