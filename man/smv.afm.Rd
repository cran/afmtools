\name{smv.afm}
\alias{smv.afm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exact and asymptotic sample mean variance
}
\description{
Compute the exact sample mean variance and the asymptotic sample mean variance of a time series.
}
\usage{
smv.afm(object,lag.max=NULL,comp=TRUE)
}
\arguments{
  \item{object}{
a fitted time-series model of ARFIMA class.
}
  \item{lag.max}{
the maximum number of lags for a exact Variance-Covariance matrix (Sowell).
}
  \item{comp}{
indicator for exact/asymptotic (respectively) sample mean variance.
}
}
\details{
The exact mean variance of \eqn{\overline{y}} is given by
 
\deqn{\mbox{Var}(\overline{y})=\frac{1}{n}\left[2\sum_{j=1}^{n-1}(1-\frac{j}{n})\gamma(j)+\gamma(0)\right]}

where the Sowell exact autocovariance function \eqn{\gamma(.)} is used. In addition, for 
large \eqn{n}, we have the asymptotic mean variance

\deqn{\mbox{Var}(\overline{y})\sim \frac{c_{\gamma}}{d(2d+1)}n^{2d-1}} 

where \eqn{c_{\gamma}=\frac{\sigma^{2}_{\varepsilon}}{\pi}\frac{|\Theta(1)|^2}{|\Phi(1)|^2}\Gamma(1-2d)\sin(\pi d)}
and \eqn{(\sigma^{2}_{\varepsilon}/\pi)|\Theta(1)|^2/|\Phi(1)|^2} is obtained by \code{\link{spectrum.arma}} function multiplied by 2.
}
\value{
  \item{v}{
a numeric value of the exact sample mean variance or the asymptotic sample mean variance such as chosen by the user.
}
}
\references{
Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.
}
\author{
Javier E. Contreras-Reyes.
}
\seealso{
\code{\link{spectrum.arma}}, \code{\link{rho.sowell}}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
mod <- arfima.whittle(y,nar=1,nma=1)
smv.afm(object=mod,comp=TRUE)
smv.afm(object=mod,comp=FALSE)
}
\keyword{variance}
\keyword{arfima}