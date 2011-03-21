\name{arfima.whittle.loglik}
\alias{arfima.whittle.loglik}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Whittle loglikelihood for ARFIMA models
}
\description{
Computes the Whittle loglikelihood of an ARFIMA model for the data.
}
\usage{
arfima.whittle.loglik(theta, series, nar, nma, fixed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
parameter value at which the loglikelihood function should be evaluated
}
  \item{series}{
a numerical vector of time series
}
  \item{nar}{
number of AR parameters \code{p}
}
  \item{nma}{
number of MA parameters \code{q}
}
  \item{fixed}{
option to fix parameters to 0. By default no parameters are fixed. Works as counting the number of parameters in the model and setting the specified ones in \code{fixed} to 0. The order of parameters is
    1) long memory parameter \code{d}; 2) AR parameters; 3) MA parameters. For example, for an \code{ARFIMA(2,d,1)} model there are 4 parameters; to fix the second AR parameter \code{phi_2} the 3rd 
	parameter has to be set equal to zero, i.e. \code{fixed=3}
}
}
\details{
The approximate maximum-likelihood estimates consist in

\code{L(theta) = -(1/2*n)[log(|Gamma_{theta})-t(Y)*Gamma^{-1}_{theta}*Y]}

with \code{(Gamma_{theta})_{ij} = gamma_{theta}(i-j) = int_{-pi}^{pi} f_{theta}(lambda)e^{(i-j)*lambda*k} dlambda}
}
\value{
	\item{L}{
value of the loglikelihood
}
	\item{sigma2}{
innovation variance; estimated by the ratio of theoretical spectrum and periodogram
}
}
\references{
Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.
}
\author{
Georg M. Goerg
}
\seealso{
\code{arfima.whittle}, \code{per}
}
\examples{
data(MammothCreek)
x=MammothCreek-mean(MammothCreek)
param=c(0.11, 0.39, -0.28)
loglik=arfima.whittle.loglik(theta=param, series=x, nar=1, nma=1)
loglik$L
loglik$sigma

####### For other examples, see 'arfima.whittle' function
}
\keyword{arfima}
\keyword{loglikelihood}% __ONLY ONE__ keyword per line
