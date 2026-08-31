# localSlurm

Stand-ins for `sbatch`, `squeue` and `scancel` for machines without a SLURM
installation, so that the `submit_*.sh` scripts of this repository work locally
without any change to them.

Put this directory in front of your PATH:

```sh
export PATH="$PWD/localSlurm:$PATH"
./submit_UncertaintyAnalysisLevante.sh -w 12 -n 100 -k 10
```

Each of the three scripts looks for the real SLURM binary on the PATH first and
hands over to it if it finds one, so having this directory on the PATH on
Levante changes nothing.

## What they do instead

`sbatch` runs the job script right away, detached, and returns immediately with
the usual `Submitted batch job <id>`, exactly as the real one does. It reads the
`#SBATCH` directives of the job script for the log files (expanding `%j`, `%x`,
`%u`, `%N`), the job name, the working directory and the mail settings, so the
job writes to the same `workOutput/<expID>/LOG.<expID>_<jobid>.log` it would on
the cluster. It also exports the `SLURM_JOB_ID`, `SLURM_JOB_NAME`,
`SLURM_SUBMIT_DIR` and friends a job script may look at, and stands in for the
`module load` line of the `.run` templates so the job uses the system R.

`squeue` prints the local jobs in the same table as the real one, and `scancel`
stops a job by signalling its whole process group, which is what reaches the
`Rscript` and its workers and not just the shell that started them.

The time limit from `--time` is *not* enforced: a wall clock that fits a Levante
node means something else on a laptop, and being killed at the Levante limit
would be worse than running long. The job log says so at the start.

## Knowing when a job is done

A detached job would otherwise finish in silence, so the `--mail-user` and
`--mail-type` directives that the `.run` templates already carry are honoured:
on the events they ask for (`BEGIN`, `END`, `FAIL`, `ALL`) the job reports back.

`LOCAL_SLURM_NOTIFY` picks the channels, `tty,mail` by default:

- `tty` writes one line plus a beep to the terminal the job was submitted from.
  This always works, and is silently skipped if that terminal is gone.
- `mail` sends the mail SLURM would send, via `mail` or `sendmail`. Whether it
  is delivered depends on the mail setup of the machine; if there is none, the
  job log says so and nothing else happens.
- `desktop` additionally raises a `notify-send` popup. Not on by default, as one
  popup per work unit would be unbearable in the work unit loops.
- `none` turns all of it off.

## Where the state lives

In `workOutput/localSlurmRegistry/`, one `key=value` file per job plus the job
id counter. `workOutput` is in `.gitignore`, so none of it is committed. Set
`LOCAL_SLURM_REGISTRY` to keep it elsewhere. Finished records are dropped after
a week; running jobs are never touched.

## Running FRIDA locally

The shims do not touch the contents of the job script. In particular they do not
reduce `numWorkers`, which the submit scripts default to a value that fits a
Levante node. Pass a number your machine can host, e.g. `-w 12`.

Two things bite when running the uncertainty analysis off the cluster:

- The `#SBATCH --time` of the runscript is ignored, see above, but a job sized
  for 128 Levante cores still takes correspondingly longer on a laptop.
- The `Rscript --max-connections=1024` in the `.run` templates makes
  `parallelly` (loaded through `caret` in `initialise.R`) think R is running
  under `R CMD check`, and it then sets `_R_CHECK_LIMIT_CORES_=TRUE`.
  `makePSOCKcluster` refuses more than two workers after that, and the run dies
  at `cluster setup...start cluster...` with
  `Error in .check_ncores(length(names)) : N simultaneous processes spawned`.
  `sbatch` therefore starts every local job with `_R_CHECK_LIMIT_CORES_=false`
  and says so in the log. Set that variable yourself before submitting if you
  want a different value; the shim leaves an existing setting alone.
