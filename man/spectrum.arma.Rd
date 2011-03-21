\name{spectrum.arma}
\alias{spectrum.arma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Theoretical spectrum of an ARMA model
}
\description{
Computes the spectrum for an ARMA model. Adapted code by David Stoffer (See references).
}
\usage{
spectrum.arma(ar = 0, ma = 0, sd.innov = 1)
}
\arguments{
  \item{ar}{
AR parameters
}
  \item{ma}{
MA parameters
}
  \item{sd.innov}{
standard deviation of input white noise
}
}
\details{
The spectrum (spectral density) of a covariance stationary \code{Y_t ~ ARMA(p,q)} process equals 
\code{f_Y(lambda)=(sigma^2_{varepsilon}/2*pi)*[|Theta(e^{i*lambda})|^2/|Phi(e^{i*lambda})|^2]}
where \code{Phi(z)} and \code{Theta(z)} are the lag polynomials of the ARMA model.
}
\value{
an \code{R} function \code{f(lambda)}.
}
\references{
Brockwell P. \& Davis R. (1991). \emph{Time Series: Theory and Methods}. Springer.

Shumway R. \& Stoffer D. (2006). \emph{Time Series Analysis and Its Applications: With R Examples}. Springer. \url{http://www.stat.pitt.edu/stoffer/tsa2/index.html}.
}
\author{
Georg M. Goerg
}
\seealso{
\code{spectrum}
}
\examples{
plot(spectrum.arma(ar=c(0.2, 0.7), ma=-0.5), 0, pi)
plot(spectrum.arma(ar=c(0, 0.7), ma=-0.5), 0, pi, add=TRUE, lty=2)
}
\keyword{spectrum}
\keyword{arma}
