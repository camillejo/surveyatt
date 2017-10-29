\name{phase.sim}
\alias{phase.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simluating survey dropout data with attrition phases}

\description{Simulate data representing dropout from a survey.  Allows for inclusion of phases of dropout by allowing for different rates of dropout throughout the survey.}

\usage{
phase.sim(n,q,drop.rate,change.point=NULL,s=10000,seed=1234,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{sample size for simulated survey}
  \item{q}{number of questions to be simulated}
  \item{drop.rate}{dropout or attrition rate at each question. Length of this argument should be one more than the number of change-points. If dropout rate is assumed constant throughout this argument should be length one.}
  \item{change.point}{question(s) in survey where the dropout rate changes. If this is left blank, dropout rate assumed constant throughout. Value entered should be an integer between 2 and q. Default is NULL.}
  \item{s}{number of simulated datasets to be created. Default is 10,000.}
  \item{seed}{to create reproducible results. Default is 1234.}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{sim.data}{list containing data frames with dropout data. First column contains IDs, second column contains survey items, and third column contains an indicator for whether or not the participant had dropped out at that question. This frame contains a row for each question for each participant (length=n*q).}
  \item{drops}{a vector of the number of participants who dropout overall within each simulated survey.}
  \item{ndrop}{a matrix holding the locations of dropouts within each simulated survey. Each row corresponds to each simulation and each column represents the survey items where participants dropped out.}

}

\author{Camille Hochheimer}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## Simulate 100 datasets with constant dropout
test1<-phase.sim(n=500,q=20,drop.rate=0.05,s=100)

## Simulate 50 datasets with two phases of dropout
test2<-phase.sim(n=500,q=20,drop.rate=c(0.05,0.1),change.point=12,s=50)

## Simulate 50 datasets with three phases of dropout
test3<-phase.sim(n=500,q=20,drop.rate=c(0.05,0.1,0.05),change.point=c(5,12),s=50)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
