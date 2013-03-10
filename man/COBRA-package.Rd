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
\code{lars}, \code{ridge}, \code{FNN}, \code{tree} and
\code{randomForest}, however any regression predictor may be
used. COBRA may natively parallelize the computations.
}
\details{
\tabular{ll}{
Package: \tab COBRA\cr
Type: \tab Package\cr
Version: \tab 0.99\cr
Date: \tab 2013-03-10\cr
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
  Nonlinear Aggregation Strategy. \url{http://hal.archives-ouvertes.fr/hal-00798579}
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
res <- COBRA(train.design = train.design,
             train.responses = train.responses,
             test = test)
print(cbind(res$predict,test.responses))
}