\name{spectrum.arfima}
\alias{spectrum.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Theoretical spectrum of an ARFIMA model
}
\description{
Computes the spectral density \code{f(lambda|Theta)} for an \code{ARFIMA(p,d,q)} model with parameter vector \code{Theta=(d,phi_1,...,phi_p,theta_1,...,theta_q,sigma_{varepsilon})}.
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
The spectral density of fractional noise \code{{X_t}} it is \code{f_X(lambda) = (sigma^2_{varepsilon}/2*pi)*(2*sin(lambda/2))^{2*d}}.

For an \code{Y_t ~ ARFIMA(p,d,q)} process it is 

\code{f_Y(lambda) = (sigma^2_{varepsilon}/2*pi)*(2*sin(lambda/2))^{2d}*[|Theta(e^{i*lambda})|^2/|Phi(e^{i*lambda})|^2]}

where \code{Phi(z)} and \code{Theta(z)} are the lag polynomials of the \code{ARMA} part of the \code{ARFIMA} model.
}
\value{
an \code{R} function \code{f(lambda)}.
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
plot(spectrum.arfima(d=0.2), 0.01, pi, ylim=c(0,1))
plot(spectrum.arfima(d=-0.2), 0.01, pi, add=TRUE, lty=2)
}
\keyword{spectrum}
\keyword{arfima}% __ONLY ONE__ keyword per line
