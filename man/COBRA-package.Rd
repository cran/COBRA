\name{COBRA-package}
\alias{COBRA-package}
%\alias{COBRA}
\docType{package}
\title{
Nonlinear Aggregation of Predictors
}
\description{The function COBRA delivers prediction outcomes for a testing sample on
the basis of a training sample and a bunch of basic regression
machines. By default, those machines are wrappers to the R packages
\code{lars}, \code{ridge}, \code{tree} and
\code{randomForest}, covering a minimal spectrum in contemporary
prediction methods for regression. However the most interesting way to use COBRA
is to use any regression method suggested by the context (see argument \code{machines}). COBRA may natively parallelize the computations (use option \code{parallel}).
}
\details{
\tabular{ll}{
Package: \tab COBRA\cr
Type: \tab Package\cr
Version: \tab 0.99.4\cr
Date: \tab 2013-07-30\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Benjamin Guedj

Maintainer: Benjamin Guedj <benjamin.guedj@upmc.fr>
}
\references{
  \url{http://www.lsta.upmc.fr/doct/guedj/index.html}

  G. Biau, A. Fischer, B. Guedj and J. D. Malley (2013), COBRA: A
  Nonlinear Aggregation Strategy. \url{http://arxiv.org/abs/1303.2236} and \url{http://hal.archives-ouvertes.fr/hal-00798579}
}
\keyword{ package }
\seealso{
\code{COBRA}
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


%% if(require(snowfall)){
%% res3 <- COBRA(train.design = train.design,
%%               train.responses = train.responses,
%%               test = test,
%%               machines = machines,
%%               machines.names = machines.names,
%%               parallel=TRUE,
%%               nb.cpus=2)
%%}
}
