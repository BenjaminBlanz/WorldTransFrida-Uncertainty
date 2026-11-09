#!/bin/bash
# testUpdateFrida.sh ####
#
# uncertainity_update_frida.sh replaces FRIDAforUncertaintyAnalysis with a fresh,
# patched checkout of name.frida_source. When git is missing or a git step or the
# patch fails it stops, and the model in FRIDAforUncertaintyAnalysis stays. It
# asks for a key at the end only when there is a terminal to press it on.
#
# Runs in a copy of the repository against a local bare repository with a
# stand-in model, so no network is needed.
#
# Run from the repository root:
#   bash developmentTools/testUpdateFrida.sh

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
# one level down, so that a cd .. out of the repository stays inside $tmp
repo="$tmp/work/repo"
mkdir -p "$repo"
git ls-files -z --cached --others --exclude-standard | tar --null -T - -cf - | tar -C "$repo" -xf -
gitAs(){ git -c user.name=test -c user.email=test@example.org "$@"; }

model="$tmp/model"
mkdir -p "$model/Data" "$model/docs" "$model/Data Processing"
echo 'a model' > "$model/FRIDA.stmx"
echo 'calibration' > "$model/Data/Calibration Data.csv"
echo 'input v1' > "$model/Data/frida_input_data.csv"
echo 'other data' > "$model/Data/other.csv"
echo 'docs' > "$model/docs/index.md"
echo 'processing' > "$model/Data Processing/run.py"
echo 'model readme' > "$model/README.md"
echo 'model licence' > "$model/LICENSE"
gitAs -C "$model" init -q -b main
gitAs -C "$model" add -A
gitAs -C "$model" commit -q -m 'v1'
gitAs clone -q --bare "$model" "$tmp/frida.git"

gitAs clone -q "$tmp/frida.git" "$tmp/patchWork"
echo 'a patched model' > "$tmp/patchWork/FRIDA.stmx"
git -C "$tmp/patchWork" diff > "$tmp/good.patch"
sed 's/^-a model$/-another model/' "$tmp/good.patch" > "$tmp/bad.patch"
cp "$tmp/good.patch" "$repo/FRIDAforAnalysis.patch"

setSource(){
	sed -i "s|^name.frida_source <-.*|name.frida_source <- '$1'|;s|^name.frida_branch <-.*|name.frida_branch <- '$2'|" "$repo/config.R"
}
setSource "$tmp/frida.git" main

# PATH without git
mkdir "$tmp/nogit"
for t in bash env grep sed rm mv mkdir rsync cp; do
	ln -s "$(command -v $t)" "$tmp/nogit/$t"
done

ok=0
fail=0
check(){
	if eval "$2"; then echo "  ok   $1"; ok=$((ok+1)); else echo "  FAIL $1"; fail=$((fail+1)); fi
}
# stdin is a pipe that stays open, as in a non-interactive shell; arguments are
# environment settings for the script
update(){
	(cd "$repo" && timeout 20 env "$@" bash uncertainity_update_frida.sh < <(sleep 25)) > "$tmp/log" 2>&1
	local status=$?
	{ echo "=== update $*"; cat "$tmp/log"; } >> "$tmp/allLogs"
	return $status
}
installed="$repo/FRIDAforUncertaintyAnalysis"
kept="$installed/kept"
# the update stopped with an error naming the cause, before anything was deleted
stopped(){
	[ $1 -ne 0 ] && [ $1 -ne 124 ] && grep -q "$2" "$tmp/log" && [ -f "$kept" ] &&
		grep -q 'a patched model' "$installed/FRIDA.stmx"
}

echo 'a first install'
update
status=$?
check 'succeeds' "[ $status -eq 0 ]"
check 'without waiting for a key' "[ $status -ne 124 ]"
check 'installs the patched model' "grep -q 'a patched model' '$installed/FRIDA.stmx'"
check 'with the data it needs' "[ -f '$installed/Data/frida_input_data.csv' ] && [ ! -e '$installed/Data/other.csv' ]"
check 'and the template' "[ -f '$installed/Data/policyParameters.csv' ]"
touch "$kept"

echo 'failures leave the installed model in place'
update PATH="$tmp/nogit"
check 'git not on PATH stops' "stopped $? 'git not found'"
setSource "$tmp/frida.git" noSuchBranch
update
check 'a fetch that fails stops' "stopped $? 'fetch'"
setSource "$tmp/frida.git" main
touch "$repo/FRIDAforUncertaintyAnalysisGit/.git/index.lock"
update
check 'a reset that fails stops' "stopped $? 'reset'"
rm "$repo/FRIDAforUncertaintyAnalysisGit/.git/index.lock"
cp "$tmp/bad.patch" "$repo/FRIDAforAnalysis.patch"
update
check 'a patch that does not apply stops' "stopped $? 'FRIDAforAnalysis.patch'"
cp "$tmp/good.patch" "$repo/FRIDAforAnalysis.patch"
mv "$repo/FRIDAforUncertaintyAnalysisGit" "$tmp/gitAway"
setSource "$tmp/noSuchSource.git" main
update
check 'a clone that fails stops' "stopped $? 'clon'"
check 'and leaves the files of this repository alone' "[ -f '$repo/README.md' ]"
mv "$tmp/gitAway" "$repo/FRIDAforUncertaintyAnalysisGit"
setSource "$tmp/frida.git" main

echo 'an update'
echo 'input v2' > "$model/Data/frida_input_data.csv"
gitAs -C "$model" commit -q -am 'v2'
gitAs -C "$model" push -q "$tmp/frida.git" main
update
check 'succeeds' "[ $? -eq 0 ]"
check 'installs the new version' "grep -q 'input v2' '$installed/Data/frida_input_data.csv'"
check 'patched' "grep -q 'a patched model' '$installed/FRIDA.stmx'"
check 'in place of the old one' "[ ! -e '$kept' ]"

echo
echo "$ok ok, $fail failed"
if [ $fail -gt 0 ]; then
	echo 'output of the updates:'
	cat "$tmp/allLogs"
	exit 1
fi
