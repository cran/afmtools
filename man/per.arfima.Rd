\name{per.arfima}
\alias{per.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Periodogram by Fast Fourier Transform
}
\description{
Computes the periodogram of a signal \code{z} by the Fast Fourier Transform (FFT).
}
\usage{
per.arfima(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
univariate signal / time series z
}
}
\value{
	\item{freq}{
Fourier frequencies \eqn{\omega_j=2\pi j/n}, where \eqn{n} is the length of \eqn{z}, and \eqn{j=1,...,[n/2]}
}
	\item{spec}{
corresponding value of the periodogram \eqn{I(\omega_j)}
}
}
\references{
Brockwell, P. and Davis, R. (1991). \emph{Time Series: Theory and Methods}. Springer. New York.
}
\author{
Wilfredo Palma, Georg M. Goerg
}
\seealso{
\code{\link{spectrum.arfima}}, \code{\link{spectrum}}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
x=per.arfima(y)
plot(x$freq, x$spec, type="l", xlab="frequencies", 
ylab="spectrum", main="Periodogram")
# cbind(x$spec,x$freq)

#### For a ARMA model:
sim <- arima.sim(n=1000, list(ar = c(0.56), 
ma = c(-0.63)), sd = sqrt(1.5))
x=per.arfima(sim)
plot(x$freq, x$spec, type="l", xlab="frequencies", 
ylab="spectrum", main="Periodogram")

#### For USAccDeaths database

par(mfrow=c(1,2))
data(USAccDeaths)
x=per.arfima(USAccDeaths)
plot(x$freq,x$spec, type="l", xlab="frequencies", 
ylab="spectrum", main="Periodogram")

USAstand=(USAccDeaths-mean(USAccDeaths))/sqrt(var(USAccDeaths))
z=per.arfima(USAstand)
plot(z$freq, z$spec, type="l", xlab="frequencies", 
ylab="spectrum", main="Periodogram")
}
\keyword{loglikelihood}
\keyword{arfima}