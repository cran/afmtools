\name{rho.arma}
\alias{rho.arma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Autocovariance function for ARMA models
}
\description{
Computes the autocorrelation function for a given ARMA process. Done via an inverse Fourier transform.
}
\usage{
rho.arma(object = NULL, lag.max = NULL, ar = object$ar, ma = object$ma)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted time-series model of ARFIMA class
}
  \item{lag.max}{
the maximum number of lags
}
  \item{ar}{
a numeric vector containing the autoregressive polynomial parameters. By default, the value of \code{object} is selected
}
  \item{ma}{
a numeric vector containing the mobile average polynomial parameters. By default, the value of \code{object} is selected
}
}
\details{
This function is used for \code{plot} function to plotting the ARFIMA class objects
}
\value{
a numerical value of size \code{lag.max} with the autocorrelations function values. In addition, a graphical output is produced if \code{plot==TRUE} with the autocorrelation function
}
\references{
Brockwell P. \& Davis R. (1991). \emph{Time Series: Theory and Methods}. Springer.
}
\author{
Georg M. Goerg
}
\seealso{
\code{\link{acf}}, \code{\link{rho.sowell}}, \code{\link{spectrum.arma}}
}
\examples{
#### For a ARFIMA model:
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
mod <- arfima.whittle(y,nar=1,nma=1)
lag=12
r.arma=rho.arma(object=mod, lag.max=lag)

n=length(y)
par(mfrow=c(1,2))
plot(r.arma,type="h",main="ACF",ylab=" ",ylim=c(-2/sqrt(n),1))
abline(h=c(0,-2/sqrt(n),2/sqrt(n)),col=c(1,2,2),lty=c(1,2,2))
acf(y,lag=lag)

#### For a ARMA model:

sim <- arima.sim(n=1000, list(ar = c(0.8897), ma = c(0.2488)), sd = 1)
r.arma=rho.arma(lag.max=6, ar=0.8897, ma = 0.2488)
n=length(sim)

par(mfrow=c(1,2))
plot(r.arma,type="h",main="ACF",ylab=" ",ylim=c(-2/sqrt(n),1))
abline(h=c(0,-2/sqrt(n),2/sqrt(n)),col=c(1,2,2),lty=c(1,2,2))
acf(sim,lag=6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arma}
\keyword{autocovariance}
