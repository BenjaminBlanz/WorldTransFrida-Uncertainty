# Shared helpers for the local SLURM stand-ins (sbatch, squeue, scancel).
#
# These scripts let the submit_*.sh scripts of this repository run on a machine
# without a SLURM installation. Each of them first looks for the real binary on
# the PATH and hands over to it if it exists, so the same checkout works both on
# Levante and locally.
#
# This file is sourced, not executed.

# location of the shim directory itself, so that the registry and the sibling
# shims are found no matter what the working directory is
localSlurmDir="$(cd -- "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")" && pwd)"

# the registry lives in workOutput, which is in .gitignore
localSlurmRegistry="${LOCAL_SLURM_REGISTRY:-${localSlurmDir}/../workOutput/localSlurmRegistry}"
localSlurmJobDir="${localSlurmRegistry}/jobs"

localSlurmUser="${USER:-$(id -un)}"
localSlurmHost="$(hostname -s 2>/dev/null || hostname)"

# hand over to a real slurm binary if one is installed ########################

# Walks the PATH and execs the first executable of the given name that does not
# live in our own directory. Returns (instead of execing) if there is none.
delegateToRealBinary(){
	local name="$1"; shift
	local dir resolved
	local IFS=:
	for dir in $PATH; do
		# an empty PATH entry means the current directory
		[ -n "$dir" ] || dir=.
		[ -x "${dir}/${name}" ] || continue
		resolved="$(cd -- "$dir" 2>/dev/null && pwd)" || continue
		# skip ourselves, otherwise we would call ourselves forever
		[ "$resolved" = "$localSlurmDir" ] && continue
		exec "${dir}/${name}" "$@"
	done
}

# tell the user that we are the stand-in and not the real thing
localSlurmNote(){
	printf 'localSlurm: %s\n' "$*" >&2
}

# the job registry ###########################################################

localSlurmInitRegistry(){
	mkdir -p "$localSlurmJobDir"
}

localSlurmJobFile(){
	printf '%s/%s.job\n' "$localSlurmJobDir" "$1"
}

# Hands out the next job id. Under flock, as the work unit submitters submit
# many jobs in a loop and two of them must never get the same id.
localSlurmNextJobId(){
	localSlurmInitRegistry
	local counter="${localSlurmRegistry}/nextJobId"
	local id
	{
		flock 9
		id=$(cat "$counter" 2>/dev/null)
		# start at 1 if the counter does not exist yet or got clobbered
		case "$id" in
			''|*[!0-9]*) id=1 ;;
		esac
		printf '%s\n' "$((id + 1))" > "$counter"
	} 9>"${counter}.lock"
	printf '%s\n' "$id"
}

# Reads a record into the localSlurmJob associative array. The record is a
# plain key=value file, one entry per line, values are single line by
# construction (paths and names, no newlines).
localSlurmReadJob(){
	local file
	file="$(localSlurmJobFile "$1")"
	[ -f "$file" ] || return 1
	unset localSlurmJob
	declare -gA localSlurmJob=()
	local line key value
	while IFS= read -r line || [ -n "$line" ]; do
		key="${line%%=*}"
		value="${line#*=}"
		[ -n "$key" ] && [ "$key" != "$line" ] || continue
		localSlurmJob["$key"]="$value"
	done < "$file"
	return 0
}

# Writes the localSlurmJob array back. Atomically, so that a squeue reading at
# the same moment never sees half a record.
localSlurmWriteJob(){
	local file tmp key
	file="$(localSlurmJobFile "${localSlurmJob[JobId]}")"
	localSlurmInitRegistry
	tmp="${file}.$$.tmp"
	for key in "${!localSlurmJob[@]}"; do
		printf '%s=%s\n' "$key" "${localSlurmJob[$key]}"
	done > "$tmp"
	mv -f "$tmp" "$file"
}

# Sets a field of a stored record without disturbing the rest of it.
localSlurmSetJobField(){
	local jobId="$1" key="$2" value="$3"
	localSlurmReadJob "$jobId" || return 1
	localSlurmJob["$key"]="$value"
	localSlurmWriteJob
}

localSlurmListJobIds(){
	local file id
	for file in "$localSlurmJobDir"/*.job; do
		[ -e "$file" ] || continue
		id="${file##*/}"
		printf '%s\n' "${id%.job}"
	done | sort -n
}

# liveness ###################################################################

# Field 22 of /proc/<pid>/stat is the process start time in clock ticks since
# boot. Recording it alongside the pid makes the liveness check immune to pid
# reuse: a recycled pid has a different start time.
localSlurmProcStartTime(){
	local stat
	stat=$(cat "/proc/$1/stat" 2>/dev/null) || return 1
	# the second field is the command name in parentheses and may itself
	# contain spaces, so cut everything up to the closing parenthesis first
	local rest="${stat##*) }"
	# shellcheck disable=SC2086 # word splitting is what we want here
	set -- $rest
	# after dropping "pid (comm) " the start time is the 20th field
	printf '%s\n' "${20}"
}

# Is the process group leader of this job still alive?
localSlurmJobIsAlive(){
	local pgid="$1" startTime="$2"
	[ -n "$pgid" ] || return 1
	[ -d "/proc/$pgid" ] || return 1
	local current
	current=$(localSlurmProcStartTime "$pgid") || return 1
	[ -n "$startTime" ] && [ "$current" != "$startTime" ] && return 1
	return 0
}

# Refreshes the state of a record that claims to be running but whose process
# is gone, e.g. after a reboot or a kill -9. Leaves the record in place so that
# squeue -t all can still show what happened.
localSlurmReapJob(){
	local jobId="$1"
	localSlurmReadJob "$jobId" || return 1
	[ "${localSlurmJob[State]}" = "RUNNING" ] || return 0
	if localSlurmJobIsAlive "${localSlurmJob[PGid]}" "${localSlurmJob[PGidStartTime]}"; then
		return 0
	fi
	localSlurmJob[State]=FAILED
	localSlurmJob[ExitCode]=${localSlurmJob[ExitCode]:--1}
	localSlurmJob[EndTime]=${localSlurmJob[EndTime]:-$(date +%s)}
	localSlurmJob[Reason]='process vanished'
	localSlurmWriteJob
}

# formatting #################################################################

# SLURM state names to the two letter codes squeue prints
localSlurmStateCode(){
	case "$1" in
		RUNNING)   printf 'R\n' ;;
		PENDING)   printf 'PD\n' ;;
		COMPLETED) printf 'CD\n' ;;
		FAILED)    printf 'F\n' ;;
		CANCELLED) printf 'CA\n' ;;
		*)         printf '%s\n' "$1" ;;
	esac
}

# seconds to the M:SS / H:MM:SS / D-HH:MM:SS progression squeue uses
localSlurmFormatElapsed(){
	local seconds="$1"
	[ -n "$seconds" ] && [ "$seconds" -ge 0 ] 2>/dev/null || seconds=0
	local days=$((seconds / 86400))
	local hours=$(((seconds % 86400) / 3600))
	local minutes=$(((seconds % 3600) / 60))
	local secs=$((seconds % 60))
	if [ "$days" -gt 0 ]; then
		printf '%d-%02d:%02d:%02d\n' "$days" "$hours" "$minutes" "$secs"
	elif [ "$hours" -gt 0 ]; then
		printf '%d:%02d:%02d\n' "$hours" "$minutes" "$secs"
	else
		printf '%d:%02d\n' "$minutes" "$secs"
	fi
}

# seconds to HH:MM:SS as used in the notification mail subjects
localSlurmFormatRunTime(){
	local seconds="$1"
	[ -n "$seconds" ] && [ "$seconds" -ge 0 ] 2>/dev/null || seconds=0
	printf '%02d:%02d:%02d\n' \
		$((seconds / 3600)) $(((seconds % 3600) / 60)) $((seconds % 60))
}

# How long has this job been running, or how long did it run?
localSlurmJobElapsed(){
	local startTime="$1" endTime="$2"
	[ -n "$startTime" ] || { printf '0\n'; return; }
	[ -n "$endTime" ] || endTime=$(date +%s)
	printf '%s\n' "$((endTime - startTime))"
}

# notification ###############################################################
# The .run templates all carry --mail-user/--mail-type. Without a mail setup a
# detached local job would finish in complete silence, so the same directives
# also drive a ping to the terminal the job was submitted from.

# Which channels to use. tty and mail by default, desktop has to be asked for
# as a popup per work unit would be unbearable in the work unit loops.
localSlurmNotifyChannels(){
	printf '%s\n' "${LOCAL_SLURM_NOTIFY:-tty,mail}"
}

localSlurmWantsChannel(){
	local channel="$1" list
	list=",$(localSlurmNotifyChannels),"
	case "$list" in
		*,none,*) return 1 ;;
		*,"$channel",*) return 0 ;;
	esac
	return 1
}

# Does the --mail-type of the job ask for this event? Same semantics as slurm:
# ALL covers everything and END fires whether the job succeeded or not.
localSlurmWantsEvent(){
	local mailType="$1" event="$2" list
	list=",$(printf '%s' "$mailType" | tr '[:lower:]' '[:upper:]'),"
	case "$list" in
		*,NONE,*) return 1 ;;
		*,ALL,*) return 0 ;;
		*,"$event",*) return 0 ;;
	esac
	return 1
}

# A one line ping written straight to the terminal the job was submitted from.
# This is a plain write to the device, it only draws on the screen and cannot
# put anything into the shell's input. Never fatal: the terminal may be long
# gone by the time a job finishes.
localSlurmNotifyTty(){
	local tty="$1" message="$2"
	[ -n "$tty" ] || return 0
	[ -c "$tty" ] && [ -w "$tty" ] || return 0
	printf '\n\a[localSlurm] %s\n' "$message" > "$tty" 2>/dev/null || true
}

# Mimics the mail slurm would send, using whatever is installed. If nothing is,
# say so in the job log once and carry on.
localSlurmNotifyMail(){
	local address="$1" subject="$2"
	# the submit scripts ship with a REPLACE_ME placeholder, say why no mail
	# went out rather than failing to send to it
	case "$address" in
		*REPLACE_ME*|'')
			printf 'localSlurm: no usable --mail-user (%s), not sending mail\n' "${address:-unset}"
			return 0
			;;
		*@*) ;;
		*)
			printf 'localSlurm: --mail-user %s is not an address, not sending mail\n' "$address"
			return 0
			;;
	esac
	if ! command -v mail >/dev/null 2>&1 && ! command -v sendmail >/dev/null 2>&1; then
		printf 'localSlurm: no mail or sendmail found, not notifying %s\n' "$address"
		return 0
	fi
	if command -v mail >/dev/null 2>&1 \
		&& mail -s "$subject" "$address" < /dev/null >/dev/null 2>&1; then
		return 0
	fi
	# mail is often a thin layer over the MTA, so try the MTA directly before
	# giving up. Whether anything is actually delivered depends on the mail
	# setup of this machine, which is none of our business.
	if command -v sendmail >/dev/null 2>&1 && {
			printf 'To: %s\n' "$address"
			printf 'Subject: %s\n' "$subject"
			printf '\n'
		} | sendmail -t >/dev/null 2>&1; then
		return 0
	fi
	printf 'localSlurm: sending mail to %s failed, is an MTA set up on this machine?\n' \
		"$address"
}

localSlurmNotifyDesktop(){
	local summary="$1" body="$2"
	command -v notify-send >/dev/null 2>&1 || return 0
	notify-send -a localSlurm "$summary" "$body" >/dev/null 2>&1 || true
}

# Sends the notification for one event (BEGIN, END or FAIL) of one job over all
# enabled channels. Reads the job record itself so that the runner only has to
# say which job and which event.
localSlurmNotify(){
	local jobId="$1" event="$2"
	localSlurmReadJob "$jobId" || return 0
	localSlurmWantsEvent "${localSlurmJob[MailType]:-NONE}" "$event" || return 0

	local name="${localSlurmJob[JobName]}"
	local state="${localSlurmJob[State]}"
	local exitCode="${localSlurmJob[ExitCode]:-0}"
	local elapsed runTime subject message
	elapsed=$(localSlurmJobElapsed "${localSlurmJob[StartTime]}" "${localSlurmJob[EndTime]}")
	runTime=$(localSlurmFormatRunTime "$elapsed")

	if [ "$event" = BEGIN ]; then
		subject="SLURM Job_id=${jobId} Name=${name} Began, Queued time 00:00:00"
		message="job ${jobId} ${name} started - log: ${localSlurmJob[StdOut]}"
	elif [ "$state" = COMPLETED ]; then
		subject="SLURM Job_id=${jobId} Name=${name} Ended, Run time ${runTime}, COMPLETED, ExitCode ${exitCode}"
		message="job ${jobId} ${name} COMPLETED in $(localSlurmFormatElapsed "$elapsed") - log: ${localSlurmJob[StdOut]}"
	else
		subject="SLURM Job_id=${jobId} Name=${name} Failed, Run time ${runTime}, ${state}, ExitCode ${exitCode}"
		message="job ${jobId} ${name} ${state} (exit ${exitCode}) after $(localSlurmFormatElapsed "$elapsed") - log: ${localSlurmJob[StdOut]}"
	fi

	localSlurmWantsChannel tty && localSlurmNotifyTty "${localSlurmJob[NotifyTty]}" "$message"
	localSlurmWantsChannel mail && localSlurmNotifyMail "${localSlurmJob[MailUser]}" "$subject"
	localSlurmWantsChannel desktop && localSlurmNotifyDesktop "SLURM job ${jobId} ${name}" "$message"
	return 0
}

# running a job ##############################################################

# The body of a local job. sbatch starts this detached in its own session, with
# stdout and stderr already redirected to the job's log files, so everything
# printed here lands in the log exactly as slurm's would.
localSlurmRunJob(){
	local jobId="$1"
	localSlurmReadJob "$jobId" || return 1

	# We are the session leader created by setsid, so our own process group is
	# the one scancel has to signal to reach the R workers as well.
	local pgid
	pgid=$(ps -o pgid= -p $$ 2>/dev/null | tr -d ' ')
	[ -n "$pgid" ] || pgid=$$
	localSlurmJob[PGid]="$pgid"
	localSlurmJob[PGidStartTime]="$(localSlurmProcStartTime "$pgid")"
	localSlurmJob[StartTime]="$(date +%s)"
	localSlurmJob[State]=RUNNING
	localSlurmWriteJob

	printf 'localSlurm: job %s (%s) started at %s on %s\n' \
		"$jobId" "${localSlurmJob[JobName]}" "$(date '+%F %T')" "$localSlurmHost"
	if [ -n "${localSlurmJob[TimeLimit]}" ]; then
		printf 'localSlurm: the time limit of %s is not enforced here\n' \
			"${localSlurmJob[TimeLimit]}"
	fi

	localSlurmNotify "$jobId" BEGIN

	# The .run templates start with a module load, which does not exist off the
	# cluster. Standing in for it keeps the log free of command not found and
	# makes the job use the system R.
	if ! command -v module >/dev/null 2>&1; then
		module(){
			printf 'localSlurm: ignoring "module %s"\n' "$*"
			return 0
		}
		export -f module
	fi

	# the part of a job's environment that a job script may reasonably look at
	export SLURM_JOB_ID="$jobId"
	export SLURM_JOBID="$jobId"
	export SLURM_JOB_NAME="${localSlurmJob[JobName]}"
	export SLURM_JOB_USER="${localSlurmJob[User]}"
	export SLURM_JOB_PARTITION="${localSlurmJob[Partition]}"
	export SLURM_JOB_ACCOUNT="${localSlurmJob[Account]}"
	export SLURM_SUBMIT_DIR="${localSlurmJob[WorkDir]}"
	export SLURM_SUBMIT_HOST="$localSlurmHost"
	export SLURM_JOB_NUM_NODES=1
	export SLURM_NNODES=1
	export SLURM_NTASKS_PER_NODE="${localSlurmJob[NTasks]}"
	export SLURM_CPUS_ON_NODE="$(nproc 2>/dev/null || echo 1)"
	export SLURM_JOB_NODELIST="$localSlurmHost"
	export SLURM_NODELIST="$localSlurmHost"
	export SLURM_CLUSTER_NAME=local
	export SLURM_PROCID=0
	export SLURM_LOCALID=0

	# The --max-connections of the .run templates makes parallelly's .onLoad
	# believe it is running under R CMD check, after which it sets
	# _R_CHECK_LIMIT_CORES_ and parallel refuses clusters of more than two
	# workers. Say no to that here, unless the environment already has an
	# opinion.
	if [ -z "${_R_CHECK_LIMIT_CORES_}" ]; then
		export _R_CHECK_LIMIT_CORES_=false
		printf 'localSlurm: setting _R_CHECK_LIMIT_CORES_=false so that R can use more than two workers\n'
	fi

	cd -- "${localSlurmJob[WorkDir]}" || return 1

	local rc=0
	# shellcheck disable=SC2086 # job arguments are intentionally split
	bash "${localSlurmJob[Script]}" ${localSlurmJob[ScriptArgs]} || rc=$?

	localSlurmReadJob "$jobId"
	localSlurmJob[EndTime]="$(date +%s)"
	localSlurmJob[ExitCode]="$rc"
	if [ "${localSlurmJob[State]}" = CANCELLED ]; then
		: # scancel already had the last word on this one
	elif [ "$rc" -eq 0 ]; then
		localSlurmJob[State]=COMPLETED
	else
		localSlurmJob[State]=FAILED
	fi
	localSlurmWriteJob

	printf 'localSlurm: job %s exited with status %s (%s) at %s\n' \
		"$jobId" "$rc" "${localSlurmJob[State]}" "$(date '+%F %T')"

	localSlurmNotify "$jobId" END
	[ "$rc" -eq 0 ] || localSlurmNotify "$jobId" FAIL
	return "$rc"
}

# Drops finished records older than a week so that the registry does not grow
# without bound. Jobs that are still going are never touched, however old.
localSlurmPruneJobs(){
	local file jobId
	for file in "$localSlurmJobDir"/*.job; do
		[ -e "$file" ] || continue
		[ -n "$(find "$file" -maxdepth 0 -mtime +7 2>/dev/null)" ] || continue
		jobId="${file##*/}"; jobId="${jobId%.job}"
		localSlurmReadJob "$jobId" || continue
		case "${localSlurmJob[State]}" in
			RUNNING|CONFIGURING|PENDING) continue ;;
		esac
		rm -f "$file"
	done
}
