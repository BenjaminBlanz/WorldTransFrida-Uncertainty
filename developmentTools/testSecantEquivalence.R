# testSecantEquivalence.R ####
#
# secant visits the same points as the three-evaluations-per-iteration form while
# evaluating the objective far fewer times. Any deterministic function shows it:
# the root must be unchanged and the call count must drop.
#
# Run from the repository root:  Rscript developmentTools/testSecantEquivalence.R

source('funParmSpace.R')
secant.new <- secant

# the three-evaluations-per-iteration form, as the reference.
secant.old <- function(fun, x0, x1, tol=1e-07, niter=1e4, doWarn=T, trace=0,
											 bound=NULL,hasToBePositive=FALSE,...){
	if(is.null(bound)){
		bound <- sign(x1-x0)*Inf
	}
	for ( i in 1:niter ){
		f0 <- fun(x0,...)
		f1 <- fun(x1,...)
		x2 <- x1-f1*(x1-x0)/(f1-f0)
		if(is.infinite(x2)||is.nan(x2)){
			return(bound)
		}
		if(hasToBePositive && x2 < 0){
			return(NA)
		}
		if(abs(fun(x2,...)) < tol || abs(x2)>abs(bound)){
			return(x2)
		}
		if(x0==x2){
			if(doWarn){
				warning("In secant cycle detected\n")
			}
			return(x2)
		}
		x0 <- x1
		x1 <- x2
	}
	if(doWarn){
		warning("In secant exceeded allowed number of iterations\n")
	}
	return(x2)
}

# a counting wrapper, standing in for the stella run
counted <- function(f){
	n <- 0
	list(fun=function(x,...){n <<- n+1; f(x,...)},
			 count=function(){n})
}

cases <- list(
	list(name='cubic, well behaved',      f=function(x){x^3-2*x-5},    x0=1,    x1=3),
	list(name='exponential',              f=function(x){exp(x)-4},     x0=0,    x1=3),
	list(name='shallow, many iterations', f=function(x){x^5-1e-3},     x0=0.01, x1=2),
	list(name='order-of-magnitude shape', f=function(x){abs(x*1e3)-1}, x0=0,    x1=1e-2),
	list(name='root far from start',      f=function(x){x-1e4},        x0=1,    x1=2),
	list(name='flat, hits niter',         f=function(x){1e-9*x+1},     x0=1,    x1=2)
)

cat(sprintf('%-28s %12s %12s %7s %7s %8s\n',
						'case','root old','root new','n old','n new','saved'))
allSame <- TRUE
totOld <- 0
totNew <- 0
for(cs in cases){
	co <- counted(cs$f)
	ro <- suppressWarnings(secant.old(co$fun,cs$x0,cs$x1,niter=100,doWarn=FALSE))
	cn <- counted(cs$f)
	rn <- suppressWarnings(secant.new(cn$fun,cs$x0,cs$x1,niter=100,doWarn=FALSE))
	# the fval attribute is carried separately; compare the roots themselves
	rn.bare <- as.numeric(rn)
	same <- isTRUE(all.equal(as.numeric(ro),rn.bare)) ||
		(is.na(ro)&&is.na(rn.bare))
	allSame <- allSame && same
	totOld <- totOld + co$count()
	totNew <- totNew + cn$count()
	cat(sprintf('%-28s %12.6g %12.6g %7d %7d %7.2fx%s\n',
							cs$name,ro,rn.bare,co$count(),cn$count(),
							co$count()/max(cn$count(),1),
							ifelse(same,'',' <-- ROOTS DIFFER')))
	# where secant evaluated the point it returns, it must report that value
	fval <- attr(rn,'fval')
	if(!is.null(fval)){
		stopifnot(isTRUE(all.equal(as.numeric(fval),cs$f(rn.bare))))
	}
}
cat(sprintf('\n%-28s %12s %12s %7d %7d %7.2fx\n',
						'TOTAL','','',totOld,totNew,totOld/totNew))
if(!allSame){
	stop('secant.new does not reproduce the roots secant.old found\n')
}
cat('\nAll roots identical; fval attribute agrees with the function where set.\n')
