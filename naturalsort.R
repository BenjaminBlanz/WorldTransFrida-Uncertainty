
# copy of the naturalsort package code
# copied here to skip installation step, which failed due to not explictly supported R 
# version, but code runs fine on modern R.

# Package: naturalsort
# Type: Package
# Title: Natural Ordering
# Version: 0.1.3
# Suggests: testthat
# Date: 2016-08-30
# Author: Kosei Abe
# Maintainer: Kosei Abe <mail@recyclebin.jp>
# 	Description: Provides functions related to human natural ordering.
# It handles adjacent digits in a character sequence as a number so that
# natural sort function arranges a character vector by their numbers, not digit
# characters. It is typically seen when operating systems lists file names. For
# example, a sequence a-1.png, a-2.png, a-10.png looks naturally ordered because 1
# < 2 < 10 and natural sort algorithm arranges so whereas general sort algorithms
# arrange it into a-1.png, a-10.png, a-2.png owing to their third and fourth
# characters.
# License: BSD_3_clause + file LICENSE
# BugReports: https://github.com/kos59125/naturalsort/issues
# RoxygenNote: 5.0.1
# NeedsCompilation: no
# Packaged: 2016-08-30 05:50:59 UTC; abe
# Repository: CRAN
# Date/Publication: 2016-08-30 12:48:28
# 
# LICENSE:
# YEAR: 2013
# COPYRIGHT HOLDER: Kosei ABE
# ORGANIZATION: RecycleBin



#' @rdname naturalsort
#' @export
naturalsort <- function(text, decreasing=FALSE, na.last=NA) {
	text[naturalorder(text, decreasing=decreasing, na.last=na.last)]
}

#' Natural Ordering Factor
#' 
#' \code{naturalfactor} creates a factor with levels in natural order.
#' 
#' @param x
#' a character vector.
#' @param levels
#' a character vector whose elements might be appeared in \code{x}.
#' @param ordered
#' logical flag that determines whether the factor is ordered.
#' @param ...
#' arguments that are passed to \code{factor} function.
#' 
#' @rdname naturalfactor
#' @export
naturalfactor <- function(x, levels, ordered = TRUE, ...) {
	text <- as.character(x)
	if (missing(levels)) {
		levels <- unique(text)
	}
	levels <- naturalsort(levels)
	factor(text, levels = levels, ordered = ordered, ...)
}

#' Natural Ordering Sort
#' 
#' Natural ordering is a kind of alphanumerical ordering.
#' \code{naturalorder} returns the order of the argument character #' vector in human natural ascending or descending order.
#' \code{naturalsort} returns the sorted vector.
#' 
#' @param text
#' a character vector to sort.
#' @param decreasing
#' logical.
#' @param na.last
#' logical. If \code{NA}, \code{NA}s will be removed of the result.
#' 
#' @return
#' For \code{naturalorder}, the results are indices of vector elements in natural order.
#' For \code{naturalsort}, the results are sorted vectors.
#' 
#' @examples
#' text <- c("a-1.png", "a-2.png", "a-10.png")
#' print(sort(text))
#' print(naturalsort(text))
#'
#' @rdname naturalsort
#' @export
naturalorder <- function(text, decreasing=FALSE, na.last=TRUE) {  # different with base::order in order or arguments
	if (!is.logical(decreasing) || length(decreasing) != 1) {
		decreasing <- as.logical(decreasing)[1]
	}
	if (is.na(decreasing)) {
		stop("'decreasing' must be either TRUE or FALSE")
	}
	if (!is.logical(na.last) || length(na.last) != 1) {
		na.last <- as.logical(na.last)[1]
	}
	if (!is.character(text)) {
		text <- as.character(text)
	}
	if (length(text) == 0L) {
		return(integer(0L))
	}
	
	sign <- (-1L) ^ decreasing
	removingNA <- is.na(na.last)  # used at last to remove NAs
	## na.last | (is.na(na.last) || na.last)
	## --------+----------------------------
	## NA      | TRUE
	## TRUE    | TRUE
	## FALSE   | FALSE
	na.last <- xor(is.na(na.last) || na.last, decreasing)
	
	## If strsplit is applied to an empty character, an empty character vector is returned.
	## Therefore, if all elements in 'text' are empty, 'maxLength' will be 0.
	## Otherwise, when there is at least one ordinal value or NA in 'text', 'maxLength' will be greater than 0.
	tokenList <- strsplit(text, "(?<=\\d)(?=\\D)|(?<=\\D)(?=\\d)", perl=TRUE)
	maxLength <- max(sapply(tokenList, length))
	if (maxLength == 0L) {  # all elements are empty ("").
		return(seq_along(text))
	}
	tokenList <- lapply(tokenList, function(tokens) c(tokens, rep("", maxLength - length(tokens))))
	tokenList <- Reduce(rbind, tokenList, matrix(, 0, maxLength))
	tokenList <- as.data.frame(tokenList, stringsAsFactors=FALSE)
	
	ranks <- lapply(tokenList, function(tokens) {
		isInteger <- grepl("^\\d+$", tokens, useBytes=TRUE)
		## stability for zero-padding equivalent values
		## (e.g. 1, 01, 001, ...)
		integers <- ifelse(isInteger, tokens, "")
		zeroPaddingRank <- sign * rank(integers, na.last=TRUE)
		## sort as integer values
		## string values sholud be ranked identically
		integers <- rep(-1, length(text))
		integers[isInteger] <- as.integer(tokens[isInteger])
		integerRank <- sign * rank(integers, na.last=TRUE)
		## sort as string values
		## integer values sholud be ranked identically
		strings <- ifelse(isInteger, "0", tokens)
		stringRank <- sign * rank(strings, na.last=na.last)
		list(stringRank, integerRank, zeroPaddingRank)
	})
	ranks <- unlist(ranks, recursive=FALSE)
	orderFunction <- sprintf("order(%s)", paste(names(ranks), collapse=","))
	result <- with(ranks, eval(parse(text=orderFunction)))
	if (removingNA) {
		result <- result[!(result %in% which(is.na(text)))]
	}
	result
}
