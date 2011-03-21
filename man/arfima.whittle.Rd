\name{arfima.whittle}
\alias{arfima.whittle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Whittle Estimation for ARFIMA models
}
\description{
Estimates the parameters of an \code{ARFIMA(p,d,q)} model via the Whittle method. Parameters can be fixed to 0 by the \code{fixed} option.
}
\usage{
arfima.whittle(series, nar = 0, nma = 0, fixed = NA)
}
\arguments{
  \item{series}{
a numerical vector with time series
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
	parameter has to be set equal to zero, i.e. \code{fixed=3}.
}
}
\value{
	\item{arfima}{ 
an object of class \code{arfima}
}
	\item{call}{ 
function call
}
	\item{method}{ 
method of Whittle estimation
}
\item{data}{ 
time serie
}
\item{par}{ 
estimated parameter \code{theta}
}
\item{loglik}{ 
value of the loglikelihood function
}
\item{d}{ 
estimated \code{d}
}
\item{ar}{ 
estimated AR parameters
}
\item{ma}{ 
estimated MA parameters
}
\item{nar}{ 
number of AR parameters \code{p}
}
\item{nma}{ 
number of MA parameters \code{q}
}
\item{fixed}{ 
indicator of possible fixed values
}
\item{sd.innov}{ 
innovation standard deviation; estimated by the ratio of theoretical spectrum and periodogram
}
\item{model}{ 
character string of the specified \code{ARFIMA(p,d,q)} model
}
\item{spectrum}{ 
theoretical spectrum implied by the parameter estimates
}
}
\references{
Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.
}
\author{
Wilfredo Palma, Georg M. Goerg
}

\seealso{
\code{\link{arima}}
}
\examples{
data(MammothCreek)
x=MammothCreek-mean(MammothCreek)
plot(x)
x11()
acf(x) ## sample ACF displays long range dependence

mod=arfima.whittle(x) # estimate just the long memory parameter
summary(mod) # highly significant; but AR and/or MA parameters have to be considered too

mod11=arfima.whittle(x, nar=1, nma=1) # AR parameter non-significant
summary(mod11)

mod01=arfima.whittle(x, nar=0, nma=1) # d increased even more
summary(mod11)
plot(mod11)

res=residuals(mod11)
acf(res)
Box.test(res, 10, type="Ljung")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arfima}
\keyword{estimation}% __ONLY ONE__ keyword per line
