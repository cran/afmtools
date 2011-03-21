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
The exact mean variance of \code{y} is given by 
\code{Var(mean(y))=(1/n)*[2*sum_{j=1}^{n-1}(1-j/n)*gamma(j)+gamma(0)]}
where the Sowell exact autocovariance function \code{gamma(.)} is used. In addition, for large \code{n}, we have the asymptotic mean variance

\code{Var(mean(y)) ~ c_{gamma}*n^{2d-1}/(d*(2d+1))} 

with \code{c_{gamma}=sigma^{2}_{varepsilon}*Gamma(1-2d)*sin(pi d)*|Theta(1)|^2/pi*|Phi(1)|^2}
and \code{(sigma^{2}_{varepsilon}/pi)*|Theta(1)|^2/|Phi(1)|^2} is obtained by \code{\link{spectrum.arma}} function multiplied by 2.
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
Javier Contreras-Reyes.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arfima}
\keyword{variance}% __ONLY ONE__ keyword per line
