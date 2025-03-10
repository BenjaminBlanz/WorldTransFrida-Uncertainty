This repository tracks scripts for working on the uncertainty analysis of FRIDA.

To run the stella simulator is needed which needs to be obtained seperately, and unzipped into Stella_Simulator_Linux (or whatever you change the location.stella variable too), you could conceivably also use another version of Stella but that is not recommended.


## Notes on the FRIDA uncertainty analysis


### Preparing files

- R is case sensitive, i.e. when preparing parameter or policy files, make sure to use the same capitalisation as in FRIDA

- The uncertainty analysis directory has to include a folder of FRIDA.
    - The FRIDA.stmx has to include additional settings for importing and exporting data from/to files. Check out the diff to an unmodified version to find out what has to be changed
    - A few files have to be added to the Data subdirectory inside the FRIDA folder. Check for additonal viles in current version compared to normal FRIDA.
    - Set the FRIDA end year to 2100 in the stmx

- If additional output variables are required, add them into the file ```frida_extra_variables_to_export_list.csv```. This should be done in the WorldTrans Google Drive WP3 folder for the uncertainty analysis, and then the file has to be exported as csv. The csv file then has to be stored under FRIDA-info/ (The file in the WP3 folder has a slightly different name)
- The model needs a ```frida_info.csv``` (inside FRIDA-info/), which includes all FRIDA model parameters and their taggs and ranges. Currently, Billy creates these files by parsing the FRIDA model for all constant parameters that have a range attached to it.

    - All model parameters that have a range, but should not be varied in the uncertainty analysis need to be tagged with "no sensi"!
    - The frida_info.csv has to be generated from the FRIDA version used for the uncertainty analyis, i.e. the one specified in the config file

- There is a list called ```frida_parameter_exclusion_list.csv```, which contains variables that are excluded in the analysis. In the future, this list should be reduced to 0 lines, as the exclusion of parameters is supposed to work via the "No Sensi" tag in FRIDA.

- If output variables are not needed, specify them in ```frida_variable_exclusion_list.csv```

- If there are model parameters, for which the a fixed range should be used (i.e. the uncertainty analysis should not search for a smaller range on its own), these parameters and the allowed min and max values should be listed in ```FRIDA-info/frida_external_ranges.csv```

- For specific policy settings, create additional files inside the ```FRIDA-configs``` folder. The policy files should be named ```policyXYZ``` for consistency.

### Running via a Slurm job on Levante
Using the uncertainty analysis on a Levante compute node generally functions the same way as in an interactive session described below. However, here everything is setup in a "runscript", which is then submitted to the cluster via SLURM.

IMPORTANT: The shell script that creates and submits the runscript is the only thing that needs to be changed. BUT: it uses the templates of the ```config.R```, ```runInitialiseData.R``` ```clusterHelp.R``` and ```runMLEandParmSpace.R``` (indicated by a "TEMPLATE_" prefix). Therefore, if one of these scripts is changed and you want this change also to be included in all future runs, the respective TEMPLATE has to be changed too (if there are no other unwanted changes to the original file, it is sufficient to e.g. do ```cp config.R TEMPLATE_config.R``` )

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
