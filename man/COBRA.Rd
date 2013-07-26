\name{COBRA}
\alias{COBRA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  %%  ~~function to do ... ~~
  COBRA
}
\description{
The function COBRA delivers prediction outcomes for a testing sample on
the basis of a training sample and a bunch of basic regression
machines. By default, those machines are wrappers to the R packages
\code{lars}, \code{ridge}, \code{tree} and
\code{randomForest}, covering a somewhat wide spectrum in contemporary
prediction methods for regression. However the most interesting way to use COBRA
is to use any regression method suggested by the context (see argument \code{machines}). COBRA may natively parallelize the computations (use option \code{parallel}).
}
\usage{
COBRA(train.design,
      train.responses,
      split,
      test,
      machines,
      machines.names,
      grid = 200,
      alpha.machines,
      parallel = FALSE,
      nb.cpus = 2,
      plots = FALSE,
      savePlots = FALSE,
      logs = FALSE,
      progress = TRUE,
      path = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train.design}{Mandatory. The design matrix for the training sample.
%%     ~~Describe \code{train.design} here~~
}
  \item{train.responses}{Mandatory. The responses vector for the
    training sample.
%%     ~~Describe \code{train.responses} here~~
}
  \item{split}{Optional. How should COBRA cut the training sample?
%%     ~~Describe \code{split} here~~
}
  \item{test}{Mandatory. The design matrix of the testing sample.
%%     ~~Describe \code{test} here~~
}
  \item{machines}{Optional. Regression basic machines provided by the
    user. This should be a matrix, whose number of rows is the length of
    the training sample (ntrain) plus the length of the testing sample (ntest), and with
    as many columns as machines. Element (i,j) of this matrix is
    assumed to be r_j(X_i), the (scalar) prediction of machine j for
    query point X_i, where i is from 1 to ntrain+ntest.
%%     ~~Describe \code{machines} here~~
}
  \item{machines.names}{Optional. If \code{machines} is provided, a list
    including the names of the machines.
%%     ~~Describe \code{machines.names} here~~
}
  \item{grid}{Optional. How precise is the discretization scheme for
    calibrating the parameters.
%%     ~~Describe \code{grid} here~~
}
\item{alpha.machines}{Optional. Coerce COBRA to use exactly
  \code{alpha.machines}. Obviously this should be a integer between 1 and the
  total number of machines.
}
  \item{parallel}{Optional. If \code{TRUE}, computations will be
    dispatched over available cpus.
%%     ~~Describe \code{parallel} here~~
}
  \item{nb.cpus}{Optional. If \code{parallel}, how many cpus should be
    used. Obviously this should not exceed the number of available cpus!
%%     ~~Describe \code{nb.cpus} here~~
}
  \item{plots}{Optional. If \code{TRUE}, explanatory plots about
    calibrating \code{epsilon} and \code{alpha} (see publication) are generated
    according to the \code{path} variable.
%%     ~~Describe \code{plots} here~~
}
\item{savePlots}{Optional. If \code{TRUE}, plots are saved as .pdf files
  according to \code{path}, otherwise they pop up in the R IDE.
%%     ~~Describe \code{plots} here~~
}
  \item{logs}{Optional. If \code{TRUE}, quadratic risks over the
    training sample for all machines and COBRA are written in the file
    "risks.txt" according to the \code{path} variable.
%%     ~~Describe \code{logs} here~~
}
  \item{progress}{Optional. If \code{TRUE}, a progress bar and final
    quadratic errors are printed.
%%     ~~Describe \code{progress} here~~
}
  \item{path}{Optional. If \code{plots} or \code{logs} are \code{TRUE},
    where should the corresponding files be created?
%%     ~~Describe \code{path} here~~
}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
  For most users, options \code{grid} and \code{split} should be set to
  their default values.
}
\value{
  Returns a list including only
  \item{predict}{The vector of predicted values.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  \url{http://www.lsta.upmc.fr/doct/guedj/index.html}

  G. Biau, A. Fischer, B. Guedj and J. D. Malley (2013), COBRA: A
  Nonlinear Aggregation Strategy. \url{http://arxiv.org/abs/1303.2236} and \url{http://hal.archives-ouvertes.fr/hal-00798579}
%% ~put references to the literature/web site here ~
}
\author{Benjamin Guedj
%%  ~~who you are~~
}
%% \note{
%%   %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  COBRA-package
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n <- 500
d <- 30
ntrain <- 400
X <- replicate(d,2*runif(n = n)-1)
Y <- X[,1]^2 + X[,3]^3 + exp(X[,10]) + rnorm(n = n, sd = .1)
train.design <- as.matrix(X[1:ntrain,])
train.responses <- Y[1:ntrain]
test <- as.matrix(X[-(1:ntrain),])
test.responses <- Y[-(1:ntrain)]

## using the default machines
if(require(lars) && require(tree) && require(ridge) &&
require(randomForest))
{
res <- COBRA(train.design = train.design,
             train.responses = train.responses,
             test = test)

print(cbind(res$predict,test.responses))
plot(test.responses,res$predict,xlab="Responses",ylab="Predictions",pch=3,col=2)
abline(0,1,lty=2)
}

## using own machines
machines.names <- c("Soothsayer","Dummy")
machines <- matrix(nr = n, nc = 2, data = 0)
machines[,1] <- Y+rnorm(n = n, sd=.1)          ## soothsayer
machines[,2] <- mean(train.responses)          ## dummy prediction, averaging train.responses

res2 <- COBRA(train.design = train.design,
              train.responses = train.responses,
              test = test,
              machines = machines,
              machines.names = machines.names)

print(cbind(res2$predict,test.responses))
plot(test.responses,res2$predict,xlab="Responses",ylab="Predictions",pch=3,col=2)
abline(0,1,lty=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ nonlinear }
\keyword{ aggregation }% __ONLY ONE__ keyword per line
\keyword{ prediction }
\keyword{ regression }
\keyword{ machinelearning }
