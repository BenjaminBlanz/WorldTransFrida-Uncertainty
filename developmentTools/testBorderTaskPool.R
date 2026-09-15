# testBorderTaskPool.R ####
#
# Min and Max run as one worker pool, with the tasks dispatched longest-first.
# The risk is in the bookkeeping: every result has to land on the parameter and
# direction it was computed for, and reordering the task list must not move it.
#
# Drives the pooling and scattering logic against a stand-in for the border
# search, so it needs neither stella nor a cluster. Run from the repository root:
#   Rscript developmentTools/testBorderTaskPool.R

source('funParmSpace.R')

nPar <- 12
parVect <- setNames(seq_len(nPar)*1.0,paste0('p',seq_len(nPar)))
parBounds <- cbind(Min=parVect-(seq_len(nPar)*2),
									 Max=parVect+(seq_len(nPar)*3))
rownames(parBounds) <- names(parVect)
parscale.parvect <- rep(1,nPar)
# make the cost spread wide and uneven
parscale.parvect[c(3,7,11)] <- 1e-3

rangeDetSkip <- rep(FALSE,nPar)
rangeDetSkip[c(5,9)] <- TRUE
notDeterminedBorders <- array(TRUE,dim=c(nPar,2),
															dimnames=list(NULL,c('Min','Max')))
# a border that is already known does not go back into the pool
notDeterminedBorders[2,'Min'] <- FALSE
notDeterminedBorders[8,'Max'] <- FALSE

# ---- the code under test, mirroring runMLEandParmSpace.R
borderTasks <- list()
for(direction in c('Min','Max')){
	for(td in which(notDeterminedBorders[,direction]&!rangeDetSkip)){
		borderTasks[[length(borderTasks)+1]] <-
			list(parIdx=td,max=(direction=='Max'),direction=direction)
	}
}
borderTaskCost <- function(tsk){
	scale <- parscale.parvect[tsk$parIdx]
	span <- abs(parBounds[tsk$parIdx,ifelse(tsk$max,2,1)]-parVect[tsk$parIdx])
	if(!is.finite(scale)||scale==0||!is.finite(span)){
		return(Inf)
	}
	return(span/scale)
}
costsBefore <- sapply(borderTasks,borderTaskCost)
borderTasks <- borderTasks[order(costsBefore,decreasing=TRUE)]

# a stand-in for findDensValBorder that encodes which task produced it
fakeSearch <- function(task,...){
	(if(task$max){1e6}else{-1e6})+task$parIdx
}
borderResults <- lapply(borderTasks,fakeSearch)

border.coefs <- array(NA_real_,dim=c(nPar,2),
											dimnames=list(NULL,c('Min','Max')))
for(task.i in seq_along(borderTasks)){
	border.coefs[borderTasks[[task.i]]$parIdx,borderTasks[[task.i]]$direction] <-
		as.numeric(borderResults[[task.i]])
}

# ---- what should have happened
expected <- array(NA_real_,dim=c(nPar,2),dimnames=list(NULL,c('Min','Max')))
for(direction in c('Min','Max')){
	for(td in which(notDeterminedBorders[,direction]&!rangeDetSkip)){
		expected[td,direction] <- (if(direction=='Max'){1e6}else{-1e6})+td
	}
}

ok <- 0
fail <- 0
check <- function(label,cond){
	if(isTRUE(cond)){cat(sprintf('  ok   %s\n',label)); ok <<- ok+1}
	else{cat(sprintf('  FAIL %s\n',label)); fail <<- fail+1}
}

nMin <- sum(notDeterminedBorders[,'Min']&!rangeDetSkip)
nMax <- sum(notDeterminedBorders[,'Max']&!rangeDetSkip)
check(sprintf('one pool holds both directions (%d min + %d max = %d tasks)',
							nMin,nMax,nMin+nMax),
			length(borderTasks)==nMin+nMax)
check('skipped parameters are not in the pool',
			!any(sapply(borderTasks,function(t){t$parIdx})%in%which(rangeDetSkip)))
check('an already known border is not in the pool',
			!any(sapply(borderTasks,function(t){t$parIdx==2&&t$direction=='Min'})))
check('every result lands on its own parameter and direction',
			identical(border.coefs,expected))
check('borders never searched stay NA',
			is.na(border.coefs[5,'Min'])&&is.na(border.coefs[9,'Max'])&&
				is.na(border.coefs[2,'Min'])&&is.na(border.coefs[8,'Max']))
check('tasks are dispatched most expensive first',
			!is.unsorted(rev(sapply(borderTasks,borderTaskCost))))
check('the expensive parameters lead the pool',
			all(sapply(borderTasks[1:6],function(t){t$parIdx})%in%c(3,7,11)))

cat(sprintf('\n%d ok, %d failed\n',ok,fail))
if(fail>0){
	stop('the border task pool does not place its results correctly\n')
}
