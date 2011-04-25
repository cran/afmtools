\name{spectrum.arma}
\alias{spectrum.arma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Theoretical spectrum of an ARMA model
}
\description{
Computes the spectrum for an ARMA model. Adapted code by David Stoffer (Shumway R. & Stoffer D., 2006).
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
The spectrum (spectral density) of a covariance stationary \eqn{Y_t \sim \mbox{ARMA}(p,q)} process equals 
\eqn{f_Y(\lambda)=\frac{\sigma^2_{\varepsilon}}{2\pi}\frac{|\Theta(e^{i\lambda})|^2}{|\Phi(e^{i\lambda})|^2}}
where \eqn{\Phi(z)} and \eqn{\Theta(z)} are the lag polynomials of the ARMA model.
}
\value{
an \code{R} function \eqn{f(\lambda)}.
}
\references{
Brockwell P. & Davis R. (1991). \emph{Time Series: Theory and Methods}. Springer.

Shumway R. & Stoffer D. (2006). \emph{Time Series Analysis and Its Applications: With R Examples}. Springer. \url{http://www.stat.pitt.edu/stoffer/tsa2/index.html}.
}
\author{
Georg M. Goerg
}
\seealso{
\code{spectrum}
}
\examples{
x=spectrum.arma(ar=c(0.2, 0.7), ma=-0.5)
u=spectrum.arma(ar=c(0, 0.7), ma=-0.5)

plot(x, xlim=c(0, pi), ylim=c(0,1))
plot(u, xlim=c(0, pi), add=TRUE, lty=2)
}
\keyword{spectrum}
\keyword{arma}
