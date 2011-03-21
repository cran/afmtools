\name{rho.sowell}
\alias{rho.sowell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exact autocovariance function for ARFIMA models
}
\description{
Computes the exact autocovariance function for the general ARFIMA process by means of the Sowell method
}
\usage{
rho.sowell(object = NULL, lag.max = NULL, ar = object$ar, ma = object$ma, d = object$d, sd.innov = object$sd.innov, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted time-series model of ARFIMA class
}
  \item{lag.max}{
the maximum number of lags for a exact Variance-Covariance matrix (Sowell)
}
  \item{ar}{
a numeric vector containing the autoregressive polynomial parameters. By default, the value of \code{object} is selected
}
  \item{ma}{
a numeric vector containing the moving average polynomial parameters. By default, the value of \code{object} is selected
}
  \item{d}{
a numeric value of fractional parameter. By default, the value of \code{object} is selected
}
  \item{sd.innov}{
a numeric value of innovation standard deviation. By default, the value of \code{object} is selected
}
  \item{plot}{
a \code{TRUE/FALSE} indicator for plotting the exact autocovariance function
}
}
\details{
For the case of \code{lag.max=h<=50}, the asymptotic autocovariance calculus is obtained by Sowell algorithm (See References).
In the case of \code{lag.max=h>50}, the asymptotic autocovariance calculus is obtained using \code{c_{gamma}*|h|^{2*d-1}}, where

\code{c_{gamma} = [sigma^{2}_{varepsilon}/pi]*[|Theta(1)|^2/|Phi(1)|^2]*Gamma(1-2*d)*sin(pi*d)}

and \code{(sigma^{2}_{varepsilon}/pi)*|Theta(1)|^2/|Phi(1)|^2} is obtained by \code{\link{spectrum.arma}} function multiplied by 2.
}
\value{
\item{av}{ 
a numeric vector of autocovariance function values
}
\item{ac}{ 
a numeric vector of autocorrelation function values
}
}
\references{
Sowell F. (1992). Maximum likelihood estimation of stationary univariate fractionally integrated time series models. \emph{Journal of Econometrics} 53, 165-188.

Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.
}
\author{
Javier Contreras-Reyes
}
\seealso{
\code{\link{acf}}, \code{\link{spectrum.arma}}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
mod <- arfima.whittle(y,nar=1,nma=1)
lag=12
rho.sowell(object=mod, lag.max=lag-1)
acf(y,lag=lag)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arfima}
\keyword{autocovariance}
