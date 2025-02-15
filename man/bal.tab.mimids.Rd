\name{bal.tab.mimids}
\alias{bal.tab.mimids}
\alias{bal.tab.wimids}
\title{
Balance Statistics for MatchIt.mice Objects
}
\description{
Generates balance statistics for \code{mimids} and \code{wimids} objects from \pkg{MatchIt.mice}.
}
\usage{
\method{bal.tab}{mimids}(mimids, 
    int = FALSE, 
    poly = 1, 
    distance = NULL, 
    addl = NULL, 
    data = NULL, 
    continuous,  
    binary, 
    s.d.denom, 
    m.threshold = NULL, 
    v.threshold = NULL, 
    ks.threshold = NULL,
    cluster = NULL, 
    abs = FALSE,
    subset = NULL,
    quick = TRUE, 
    ...)
}

\arguments{
  \item{mimids}{
a \code{mimids} or \code{wimids} object; the output of a call to \code{matchitmice()} or \code{weightitmice()} from the \pkg{MatchIt.mice} package.
}
  \item{int}{
\code{logical} or \code{numeric}; whether or not to include 2-way interactions of covariates included in \code{covs} and in \code{addl}. If \code{numeric}, will be passed to \code{poly} as well. In older versions of \pkg{cobalt}, setting \code{int = TRUE} displayed squares of covariates; to replicate this behavior, set \code{int = 2}.
}
  \item{poly}{
\code{numeric}; the highest polynomial of each continuous covariate to display. For example, if 2, squares of each continuous covariate will be displayed (in addition to the covariate itself); if 3, squares and cubes of each continuous covariate will be displayed, etc. If 1, the default, only the base covariate will be displayed. If \code{int} is numeric, \code{poly} will take on the value of \code{int}.
}
  \item{distance}{
Optional; either a vector or data.frame containing distance values (e.g., propensity scores) for each unit or a string containing the name of the distance variable in \code{data}. Note that the distance measure generated by \code{matchitmice()} is automatically included.
}
  \item{addl}{
an optional data frame or the quoted names of additional covariates for which to present balance. These may be covariates included in the original dataset but not included in the call to \code{matchit()}. If variable names are specified, \code{bal.tab()} will look first in the argument to \code{data}, if specified, and next in the \code{matchit} object.
}
  \item{data}{
an optional data frame containing variables that might be named in arguments to \code{distance}, \code{addl}, and \code{cluster}. If \code{method = "exact"} or \code{approach = "across"} were used in the call to \code{matchitmice()}, an argument to \code{data} should be specified. See Examples.
}
  \item{continuous}{
whether mean differences for continuous variables should be standardized ("std") or raw ("raw"). Default "std". Abbreviations allowed. This option can be set globally using \code{\link{set.cobalt.options}}.
}
  \item{binary}{
whether mean differences for binary variables (i.e., difference in proportion) should be standardized ("std") or raw ("raw"). Default "raw". Abbreviations allowed. This option can be set globally using \code{\link{set.cobalt.options}}.
}
  \item{s.d.denom}{
whether the denominator for standardized differences (if any are calculated) should be the standard deviation of the treated group ("treated"), the standard deviation of the control group ("control"), or the pooled standard deviation ("pooled"), computed as the square root of the mean of the group variances. Abbreviations allowed. The default is "treated" for \code{mimids} objects and "pooled" for \code{wimids} objects.
}
  \item{m.threshold}{
a numeric value for the threshold for mean differences. .1 is recommended. 
}
  \item{v.threshold}{
a numeric value for the threshold for variance ratios. Will automatically convert to the inverse if less than 1.
}
  \item{ks.threshold}{
a numeric value for the threshold for Kolmogorov-Smirnov statistics. Must be between 0 and 1. 
}
  \item{cluster}{
either a vector containing cluster membership for each unit or a string containing the name of the cluster membership variable in \code{data} or the \code{mimids} or \code{wimids} object. See \code{\link{bal.tab.cluster}} for details.
}
  \item{abs}{
\code{logical}; whether displayed balance statistics should be in absolute value or not. 
}
  \item{subset}{
a \code{logical} vector denoting whether each observation should be included. It should be the same length as the variables in the original call to \code{matchitmice()} or \code{weightitmice()} (i.e., one for each individual for each imputation). \code{NA}s will be treated as \code{FALSE}. This can be used as an alternative to \code{cluster} to examine balance on subsets of the data.
}
  \item{quick}{
\code{logical}; if \code{TRUE}, will not compute any values that will not be displayed. Set to \code{FALSE} if computed values not displayed will be used later.
}
  \item{...}{
further arguments to control display of output. See \link[=options-display]{display options} for details.
}
}
\details{
\code{bal.tab.mimids()} generates a list of balance summaries for the \code{mimids} or \code{wimids} object given, and functions similarly to \code{summary.mimids()} and \code{summary.wimids()} in \pkg{MatchIt.mice}.

The threshold values (\code{m.threshold}, \code{v.threshold}, and \code{ks.threshold}) control whether extra columns should be inserted into the Balance table describing whether the balance statistics in question exceeded or were within the threshold. Including these thresholds also creates summary tables tallying the number of variables that exceeded and were within the threshold and displaying the variables with the greatest imbalance on that balance measure.

}
\value{
If clusters are not specified, an object of class \code{"bal.tab.imp"} containing balance summaries for each imputation and a summary of balance across imputations. See \code{\link{bal.tab.imp}} for details.

If clusters are specified, an object of class \code{"bal.tab.imp.cluster"} containing summaries between and across all clusters and imputations.
}

\author{
Noah Greifer 
}

\seealso{
\code{\link{bal.tab}} for details of calculations.
}
\examples{
library(MatchIt.mice)
data("lalonde_mis", package = "cobalt")

#Imputing the missing data
imp <- mice(lalonde_mis, m = 5)

#Matching using within-imputation propensity scores
mim.out1 <- matchitmice(treat ~ age + educ + race + 
                       married + nodegree + re74 + re75, 
                       data = imp, approach = "within")
bal.tab(mim.out1)

#Matching using across-imputation average propensity scores
mim.out2 <- matchitmice(treat ~ age + educ + race + 
                       married + nodegree + re74 + re75, 
                       data = imp, approach = "across")
bal.tab(mim.out2, data = complete(imp, "long"))

#Weighting using within-imputation propensity scores
## Note: These weights are incorrect as of 
## version 3.0.1 of MatchIt.mice
wim.out <- weightitmice(treat ~ age + educ + race + 
                       married + nodegree + re74 + re75, 
                       data = imp, approach = "within")
bal.tab(wim.out)

}