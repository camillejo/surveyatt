\name{glmm.qbyqdrop}
\alias{glmm.qbyqdrop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting a GLMM and searching for significant dropout}
\description{Apply a generalized linear mixed model to the survey data generated from the phase.sim function or other survey data formatted in this fashion. This function also applies a successive differences contrast to test question by question differences in the proportion of dropout.}
\usage{
glmm.qbyqdrop<-function(dta, q, alpha=0.05,...)
}
\arguments{
  \item{dta}{data.frame object to be used for GLMM.}
  \item{q}{number of questions in survey.}
  \item{alpha}{significance level. Default is 0.05.}
}
\details{This code takes a data.frame object and applies a generalized linear mixed model (GLMM) incorporating a subject level random effect to account for the repeated measures nature of survey data. The \code{glmer} function from the \emph{lme4} package is used to apply the GLMM with the following formula:\code{Drop~(Q)+(1|ID)}. \code{Drop} is a binary outcome, therefore the binomial family is used. The \code{contr.sdif} function from the \emph{MASS} package is used to create the successive differences contrast. The resulting p-values for each comparison are reported as well as the first and last comparison at which a significant difference is found. }
\value{
  \item{pval}{vector of corrected p-values from each comparison.}
  \item{none}{indicator for whether or not there are no change-points. If 1, then there are no change-points detected in the data. If 0, then at least one change-point was detected.}
  \item{nchange}{number of change-points detected.}
  \item{q.d1}{vector indicating the first question where a significant result was detected.}
  \item{q.d2}{vector indicating the last question where a significant result was detected.}
%% ...
}
\references{
add reference for FDR
}
\author{
Camille Hochheimer
}
%\note{
%%  ~~further notes~~
%}

\section{Warning }{Due to previous convergence issues, the following settings are invoked within the \code{glmer} function:\code{control=glmerControl(optimizer="bobyqa"), nAGQ=0}. Consequently, the results for the GLMM are not as exact as desired.}

\seealso{
\code{\link{glmer}} for details on GLMM
}

\examples{
## Simulate data and apply GLMM
test.dta<-phase.sim(n=500,q=20,drop.rate=c(0.05,0.1),change.point=10,s=1)
glmm.out<-glmm.qbyqdrop(test.dta$sim.data[[1]],20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
