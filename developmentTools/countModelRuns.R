# countModelRuns.R ####
#
# Counts stella runs, per section of runMLEandParmSpace.R.
#
# Wall clock in the parscale determination and the range finding is the number of
# stella runs divided by the number of workers, so that count is the only honest
# measure of whether a change to either has helped. Timings taken on a machine
# that is also running a job measure the other job.
#
# Everything here is development only. The one thing that has to live in the
# production code is a guarded increment at the top of runFRIDASpecParms in
# funRunFRIDA.R, which does nothing unless devTools.countModelRuns is TRUE.
#
# ---- how to use it ----
#
# In an R session that has already sourced config.R and clusterHelp.R, so that cl
# exists:
#
#   source('developmentTools/countModelRuns.R')
#   devToolsCountModelRuns(cl)            # switch counting on, here and on workers
#   devToolsMarkSection('parscale')       # before the section of interest
#   ... run the section ...
#   devToolsMarkSection('range finding')  # closes the previous section
#   ... run the section ...
#   devToolsReportModelRuns(cl)           # closes the last one and prints
#
# Or, to instrument a whole run without editing it, put the first two calls after
# the source('clusterHelp.R') line in runMLEandParmSpace.R and the report at the
# end. Those three lines are the whole footprint.

devTools.sectionCounts <- list()
devTools.currentSection <- NULL

# counting ####
# Switch counting on in this session and on every worker.
devToolsCountModelRuns <- function(cl=NULL){
	devTools.countModelRuns <<- TRUE
	devTools.modelRunCount <<- 0
	if(!is.null(cl)){
		parallel::clusterEvalQ(cl,{
			devTools.countModelRuns <- TRUE
			devTools.modelRunCount <- 0
		})
	}
	cat('counting stella runs\n')
	invisible(TRUE)
}

# Runs since the last reset, this session's plus every worker's. The workers hold
# their own counts, so this has to go and ask them.
devToolsModelRunCount <- function(cl=NULL){
	own <- if(exists('devTools.modelRunCount')){devTools.modelRunCount}else{0}
	workers <- if(is.null(cl)){
		0
	} else {
		sum(unlist(parallel::clusterEvalQ(cl,
			if(exists('devTools.modelRunCount')){devTools.modelRunCount}else{0})))
	}
	return(own+workers)
}

devToolsResetModelRunCount <- function(cl=NULL){
	devTools.modelRunCount <<- 0
	if(!is.null(cl)){
		parallel::clusterEvalQ(cl,devTools.modelRunCount <- 0)
	}
	invisible(TRUE)
}

# sections and reporting ####
# Close the section that was running, if any, and open a new one. Call with NULL
# to close the last section without opening another.
devToolsMarkSection <- function(name,cl=NULL){
	if(!is.null(devTools.currentSection)){
		n <- devToolsModelRunCount(cl)
		devTools.sectionCounts[[devTools.currentSection]] <<-
			(if(is.null(devTools.sectionCounts[[devTools.currentSection]])){0}
			 else{devTools.sectionCounts[[devTools.currentSection]]})+n
		cat(sprintf('  %-28s %8d stella runs\n',devTools.currentSection,n))
	}
	devTools.currentSection <<- name
	devToolsResetModelRunCount(cl)
	invisible(TRUE)
}

# Close the open section and print the tally. Writes it beside the run's other
# output when location.output is available, so two runs can be compared later.
devToolsReportModelRuns <- function(cl=NULL,file=NULL){
	devToolsMarkSection(NULL,cl)
	if(length(devTools.sectionCounts)==0){
		cat('no sections were marked\n')
		return(invisible(NULL))
	}
	counts <- unlist(devTools.sectionCounts)
	cat('\nstella runs by section\n')
	cat(sprintf('%s\n',strrep('-',44)))
	for(nm in names(counts)){
		cat(sprintf('%-28s %8d  %5.1f%%\n',nm,counts[[nm]],100*counts[[nm]]/sum(counts)))
	}
	cat(sprintf('%-28s %8d\n','TOTAL',sum(counts)))
	if(is.null(file)&&exists('location.output')){
		file <- file.path(location.output,'devTools-modelRunCounts.csv')
	}
	if(!is.null(file)){
		utils::write.csv(data.frame(section=names(counts),
																stellaRuns=as.integer(counts)),
										 file,row.names=FALSE)
		cat(sprintf('\nwritten to %s\n',file))
	}
	invisible(counts)
}

# Two of those csv files side by side, which is how a change is shown to have
# helped: same config, before and after.
devToolsCompareModelRuns <- function(before,after){
	b <- utils::read.csv(before)
	a <- utils::read.csv(after)
	m <- merge(b,a,by='section',suffixes=c('.before','.after'),all=TRUE)
	m$stellaRuns.before[is.na(m$stellaRuns.before)] <- 0
	m$stellaRuns.after[is.na(m$stellaRuns.after)] <- 0
	m$saving <- ifelse(m$stellaRuns.after>0,
										 m$stellaRuns.before/m$stellaRuns.after,NA)
	cat(sprintf('%-28s %10s %10s %8s\n','section','before','after','saving'))
	for(r.i in seq_len(nrow(m))){
		cat(sprintf('%-28s %10d %10d %7.2fx\n',m$section[r.i],
								m$stellaRuns.before[r.i],m$stellaRuns.after[r.i],m$saving[r.i]))
	}
	cat(sprintf('%-28s %10d %10d %7.2fx\n','TOTAL',
							sum(m$stellaRuns.before),sum(m$stellaRuns.after),
							sum(m$stellaRuns.before)/max(sum(m$stellaRuns.after),1)))
	invisible(m)
}
