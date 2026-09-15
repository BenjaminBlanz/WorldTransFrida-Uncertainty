# runAllTests.R ####
#
# Every test in developmentTools, in one go. None of them needs stella, a cluster
# or calibration data, so the whole suite runs in seconds. Run it after any change
# to the parscale determination or the border search.
#
# Run from the repository root:  Rscript developmentTools/runAllTests.R

testFiles <- sort(list.files('developmentTools',pattern='^test.*\\.R$',full.names=TRUE))
if(length(testFiles)==0){
	stop('no tests found in developmentTools\n')
}

results <- data.frame(test=basename(testFiles),passed=NA,stringsAsFactors=FALSE)
for(t.i in seq_along(testFiles)){
	cat(sprintf('\n%s\n%s\n',testFiles[t.i],strrep('-',nchar(testFiles[t.i]))))
	# a separate process each, so one test's stand-in negLLike cannot leak into
	# the next one's environment
	status <- system2('Rscript',testFiles[t.i])
	results$passed[t.i] <- status==0
}

cat(sprintf('\n%s\n',strrep('=',60)))
for(r.i in seq_len(nrow(results))){
	cat(sprintf('%-40s %s\n',results$test[r.i],
							ifelse(results$passed[r.i],'passed','FAILED')))
}
cat(sprintf('%d of %d passed\n',sum(results$passed),nrow(results)))
if(!all(results$passed)){
	quit(status=1)
}
