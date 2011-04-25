\name{rho.sowell}
\alias{rho.sowell}
\title{
Exact autocovariance function for ARFIMA models
}
\description{
Computes the exact autocovariance function for the general ARFIMA process by means of the Sowell method
}
\usage{
rho.sowell(object = NULL, lag.max = NULL, ar = object$ar, 
           ma = object$ma, d = object$d, 
		   sd.innov = object$sd.innov, plot = TRUE)      
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
For the case of \eqn{\code{lag.max}=h\leq50} (See Palma (2007) example), the asymptotic 
autocovariance calculus is obtained by Sowell algorithm (See References).
In the case of \eqn{\code{lag.max}=h>50}, the asymptotic autocovariance calculus is obtained 
using \eqn{c_{\gamma}|h|^{(2d-1)}}, where

\deqn{c_{\gamma}=\frac{\sigma^{2}_{\varepsilon}}{\pi}\frac{|\Theta(1)|^2}{|\Phi(1)|^2}\Gamma(1-2d)\sin(\pi d)}

and \eqn{(\sigma^{2}_{\varepsilon}/\pi)|\Theta(1)|^2/|\Phi(1)|^2} is obtained by 
\code{\link{spectrum.arma}} function multiplied by 2.
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

# Particular case (Palma, 2007; page 69)

av.arfima1d1 <- function(h,d,phi,theta) {
	C <- function(d,h,phi) { # function defined in (3.23), page 48. 
	gamma0=gamma(1-2*d)*gamma(h+d)/(gamma(1-d)*gamma(d)*gamma(1+h-d))
	res=gamma0*(phi^(2)*hypergeo(d+h,1,1-d+h,phi)+ hypergeo(d-h,1,1-d-h,phi)-1)
	return(res)
	}
	av=rep(NA,h+1)
	for(i in 0:h){
		if(i<=50){
			a=theta*C(d,-i,-phi)+(1+theta^2)*C(d,1-i,-phi)+theta*C(d,2-i,-phi)
			b=phi*(phi^2-1)
			av[i+1]=Re(a/b)
		}
		else {
			a=(abs(i)^(2*d-1))*((1+theta)^2)*gamma(1-2*d)*sin(pi*d)
			b=pi*(1+phi)^2
			av[i+1]=a/b
		}
	}
	return(round(av,4))
}

av.arfima1d1(h=998, d=0.4, phi=0.5, theta=0.2)
rho.sowell(lag.max=998, ar=0.5, ma=0.2, d=0.4, sd.innov=1)$av

av.arfima1d1(h=998, d=0.499, phi=-0.9, theta=-0.3)
rho.sowell(lag.max=998, ar=-0.9, ma=-0.3, d=0.499, sd.innov=1)$av
}
\keyword{autocovariance}
