\name{spectrum.arfima}
\alias{spectrum.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Theoretical spectrum of an ARFIMA model
}
\description{
Computes the spectral density \eqn{f(\lambda|\Theta)} for an ARFIMA\eqn{(p,d,q)} model with
parameter vector \eqn{\Theta=(d,\phi_1,...,\phi_p,\theta_1,...,\theta_q,\sigma_{\varepsilon})}.
}
\usage{
spectrum.arfima(d = 0, ar = 0, ma = 0, sd.innov = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
long memory parameter between -1 and 0.5
}
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
The spectral density of fractional noise \eqn{\{X_t\}} it is 

\deqn{f_X(\lambda)=\frac{\sigma^2_{\varepsilon}}{2\pi}(2\sin(\lambda/2))^{2d}}

For an \eqn{Y_t ~ \mbox{ARFIMA}(p,d,q)} process it is 

\deqn{f_Y(\lambda)=\frac{\sigma^2_{\varepsilon}}{2\pi}(2\sin(\lambda/2))^{2d}\frac{|\Theta(e^{i\lambda})|^2}{|\Phi(e^{i\lambda})|^2}}

where \eqn{\Phi(z)} and \eqn{\Theta(z)} are the lag polynomials of the ARMA part of the ARFIMA model.
}
\value{
an \code{R} function \eqn{f(\lambda)}.
}
\references{
Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.
}
\author{
Georg M. Goerg
}
\seealso{
\code{\link{spectrum}}
}
\examples{
x=spectrum.arfima(d=0.2)
u=spectrum.arfima(d=-0.2)

plot(x, xlim=c(0.01, pi), ylim=c(0,1))
plot(u, xlim=c(0.01, pi), add=TRUE, lty=2)
}
\keyword{spectrum}

