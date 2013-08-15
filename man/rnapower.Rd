\name{rnapower}
\alias{rnapower}
\title{RNA-seq power computation}
\description{Sample size and power computation for RNA-seq studies}
\usage{
rnapower(depth, n, n2 = n, cv, cv2 = cv, effect, alpha, power)
}
\arguments{
  \item{depth}{average depth of coverage for the transcript or gene of
    interest. Common values are 5-20, any numeric value >0 is valid. }
  \item{n}{sample size in group 1 (or both)}
  \item{n2}{sample size in group 2}
  \item{cv}{biological coefficient of variation in group 1 (or both).}
  \item{cv2}{biological coefficient of variation in group 2}
  \item{effect}{target effect size}
  \item{alpha}{size of the test statistic (two sided), i.e., the false
  positive rate, a number between 0 and 1}
  \item{power}{power of the test, i.e., the fraction of true positives
    that will be detected, a number between 0 and 1}
}
\details{
  The error in an RNA-seq study arises from two causes: the technical
  variability of the sequencing itself, and the within group biological
  variability of the specimens being compared.  For human samples this
  latter will often be .4-1, for inbred animal lines values are commonly
  .1 or less.  That is, we might see a 2-fold difference
  in mean expression between treated and untreated samples, and at the
  same time a 50\% variation in expression between samples within the
  control or treatment group.  This would correspond to an effect size
  of 2 and a CV of 0.5.  As a general rule, sequencing depths of more
  than 5/CV^2 will lead to only minor gains in study efficiency and/or
  power, whereas addition of further samples is always efficatious.

  Depth is a required argument; any one of the others may be left
  missing and the function will solve for it.  An argument may be a
  vector, in which case a vector of values is returned.  If multiple
  arguments are vectors a matrix or array of results is returned.

  Common values for \code{alpha} and \code{power} would be .05 and .9.
  The \code{effect} parameter specifies the biological effect that, if
  it were true, the experimenter would want to be able to detect; values
  of 1.5 to 2 are commonly used.  The statements that A has twice the
  expression of B and that B has half that of A are symmetric,
  likewise values of .5 and 2 for \code{effect} will yeild the same
  result.  The estimated CV of expression within group may be the
  most difficult parameter to choose; see the vignette for an in depth
  discussion of this along with recommended values for diffrent types
  of data.

  By default the samples sizes in the two groups are assumed to be
  equal (n2=n) if a second sample size is not given.  Likewise the
  coefficients of variation in the two groups
  are assumed to be equal if \code{cv2} is not specified.
}
\value{a vector, matrix, or array of values.  If the function was used
  with all the arguments except \code{n} supplied then the result will
  be the needed samples size; if it were run with all arguments except
  \code{power} the result will contain power estimates, etc.}
\author{Terry Therneau and Steven Hart}
\references{
  Steven Hart, Terry Therneau, Yuji Zhang, Gregory Poland, and
  Jean-Pierre Kocher.
  Calculating sample size estimates for RNA sequencing data.
  J. Comp. Biology, 2013, in press.
}  
\examples{
# What would the power be for a study with 12 per group, to detect a
#  2-fold change, given deep (50x) coverage?
rnapower(50, cv=.6, effect=2, n=12, alpha=.05)

# Compute a sample size for several combinations of parameters
temp <- rnapower(10, cv=.5, effect=c(1.5, 2), alpha=c(.05, .01, .001),
                   power=c(.8, .9))
round(temp,1)
# Result is an array with dimensions of effect size, alpha, and power
# which contains the sample size per group for each combination.
}
\keyword{statistics}
