This repository tracks scripts for working on the uncertainty analysis of FRIDA.

To run the stella simulator is needed which needs to be obtained seperately, and unzipped into Stella_Simulator_Linux (or whatever you change the location.stella variable too), you could conceivably also use another version of Stella but that is not recommended.


## Notes on the FRIDA uncertainty analysis


### Preparing files

- R is case sensitive, i.e. when preparing parameter or policy files, make sure to use the same capitalisation as in FRIDA

- Use uncertainity_update_frida.sh to ensure you have access to the latest model and to prep the model for use with this program.  This bash script takes care of all sub-bullets here
    - The uncertainty analysis directory has to include a folder of FRIDA.
        - The FRIDA.stmx has to include additional settings for importing and exporting data from/to files. Check out the diff to an unmodified version to find out what has to be changed
        - A few files have to be added to the Data subdirectory inside the FRIDA folder. Check for additonal files in current version compared to normal FRIDA.
        - The end time for the analysis can be controlled via the "Simulation End" variable in FRIDA

- If additional output variables are required, add them into the file ```frida_extra_variables_to_export_list.csv```. This should be done in the WorldTrans Google Drive WP3 folder for the uncertainty analysis, and then the file has to be exported as csv. The csv file then has to be stored under FRIDA-info/ (The file in the WP3 folder has a slightly different name)
- The model needs a ```frida_info.csv``` (inside FRIDA-info/), which includes all FRIDA model parameters and their taggs and ranges. Currently, Billy creates these files by parsing the FRIDA model for all constant parameters that have a range attached to it.

    - All model parameters that have a range, but should not be varied in the uncertainty analysis need to be tagged with "no sensi"!
    - The frida_info.csv has to be generated from the FRIDA version used for the uncertainty analyis, i.e. the one specified in the config file

- There is a list called ```frida_parameter_exclusion_list.csv```, which contains variables that are excluded in the analysis. In the future, this list should be reduced to 0 lines, as the exclusion of parameters is supposed to work via the "No Sensi" tag in FRIDA. Variables which are excluded by the ranging algorithim are included here as a part of program operation.  Variables fail ranging generally because they have no impact on any variable which has calibration data.  This can happen if its a parameter which control a future policy assumption.

- If output variables are not needed, specify them in ```frida_variable_exclusion_list.csv``` by default all variables in the calibration data set are included.

- If there are model parameters, for which the a fixed range should be used (i.e. the uncertainty analysis should not search for a smaller range on its own), these parameters and the allowed min and max values should be listed in ```FRIDA-info/frida_external_ranges.csv```

- For specific policy settings, create additional files inside the ```FRIDA-configs``` folder. The policy files should be named ```policyXYZ``` for consistency.

### Running via a Slurm job on Levante
Using the uncertainty analysis on a Levante compute node generally functions the same way as in an interactive session described below. However, here everything is setup in a "runscript", which is then submitted to the cluster via SLURM.

IMPORTANT: The shell script that creates and submits the runscript is typically the only thing that needs to be changed. However, there might be things that cannot be inserted there. Then these have to be changed directly in the respective R scripts. It might be useful to then adapt the shell script to also allow these changes there too.

The shell script that creates the scripts is currently called **```submit_UncertaintyAnalysisLevante.sh```**. The things that should be changed in there are the following:
- Set an experiment ID. This will be the name used for the scripts and the workOutput
- Define the number of Workers and the number of Samples.
- Set the time limit of the slurm job in hours and minutes. If new parameter scales and ranges should be defined, set this to several hours (max: 8:00 h). If only the FRIDA sampling shall be done, this can be less (??)
- If you want to use the scales and ranges from a different experiment, a respective switch can be set to true and the reference experiment ID has to be given.
- If you designed specific new input files for the analysis (e.g. a policy file), you can define them here too.


### Running in an interactive session

#### R setup and packages
- Ideally, use an R version >4.4, and then use the option ```--max-connections=1024``` when starting R
- add the option ```--no-site-file``` when starting R to avoid the package load messages like ```Loading required package: graphics```. These are related to packages 


- It is possible that the packages listed in ```initialise.R``` need to be installed first (only once). Also the packages ```imputeTs```, ```matrixcalc``` and ```caret``` need to be installed. To do that, in an R sessions type ```install.packages('packagename')```. 
    - The ```imputeTs``` package is currently (04. March 2025) not loadable on Levante with R v4.4. Can comment this package out in runInitialiseData.R, because it is only needed for likelihood analysis, which we currently don't do.

- prepare the ```config.R``` file
    -  set the correct name for the FRIDA folder (l. 122)
    -  for a test run, reduce the number of samples by setting ```numSample <- 1e4``` (l. 51)
    -  For specific policy runs, set the policy file name (l. 117)

- set the number of parallel threads in the R session by setting ```numWorkers <- 250``` (250 should work on a Levante node that has 128 cores (128x2, minus a few, which are always occupied). When using R version < 4.4 or on a smaller node, use 120). This can also be set inside ```config.R``` instead.

#### Running the runMLEandParmSpace.R
- Before running a new analysis (e.g. with different policy), where you want to use the same parameter ranges as in a previous analysis, some files can be copied to the new directory and then the script automatically skips the range determination. The files that have to be copied are listed below. This potentially has to be done during the first pause after having already started the runMLE script, as only then the new directory is already created.
    - all ```sample*``` files
    - all ```sigma*``` files
    - ```parscale.RDS```

- execute ```source('runMLEandParmSpace.R')```, this automatically sources the following files:
    - ```config.R, initialiseData.R, setupTMPFS.R``` -> These can also be exectuded one by one beforehand for testing
    - At the start of the execution one has to check for Errors in the log output and hit enter twice.

- (execute ```source('runPlotAllRuns.R')```)

# Notes on the FRIDA policy analysis

## Before First Run
### Preparation of the environment
- Place Stella Simulator for Linux files into the ```Stella_Simulator_Linux``` folder
- Run ```policy_update_frida.sh``` to check out the version of FRIDA that includes the rpresentative member subsample

### Technical Configuration
- ```PolicyAnalysis.run``` and ```PolicyAnalysisMerge.run```
    - The headers of these files (the lines starting with ```$SBATCH```) will be set automatically, but remove any that are incompatible with your installation of SLURM.
    - Change the module load command to load a version of R available on your machine, or remove it if R is always available. An R version greater than 4.4 is required if you intend to run with more than 125 workers per batch job.
- ```submit_PolicyAnalysisSLURMjob.sh```
    - Modify the account, partition, email, etc. in the ```Change config here``` section.

### Prepare Single Domain Policy Files
- Execute ```python3 createFilesForPolicyAnalysis.py``` in ```policyFileGeneration```
- This will create a folder called ```PolicyFolder``` which contains the single domain policy files. 

## For each Run
### Run Configuration
- Copy the single domain policy files from which policies should be sampled from ```policyFileGeneration/PolicyFolder``` to ```policy-singleDomainPolicyMatrices```.
- edit ```configPolicyAnalysis.R```
    - Specify the run configuration. This file contains documentation in comments
    - Note that policies to sample are specified toward the bottom of the file.
### Run
- Execute ```Rscript runPolicyAnalysisSLURMsetup.R```
    - This splits the total ensemble into work units suitable for SLURM batch jobs and creates corresponding run files and output folder structures.
- Execute ```Rscript runPolicyAnalysisSLURMrunner.R```
    - This submitts the work unit batch jobs to SLURM and reports the state of the chunks.
    - This can be run repeatedly, it will report on the current state of the batch jobs, and submit further work units as the queue empties, until all jobs have been submitted
    - Once this script reports all jobs have been completed proceed to the next step.
- Execute ```Rscript runPolicyAnalysisMerger.R```
    - This will merge the outputs from the work units into single files per variable. (Variables to export can be specidied as in the uncertainty analysis above)
- Exceute ```Rscript runPolicyAnalysisPlots.R```
    - Creates some figures, amend as desired, or simply work on the per variable exported files yourself.
