\name{var.afm}
\alias{var.afm}
\title{
Parameter Variance-Covariance matrix for ARFIMA models
}
\description{
Method to estimate the Parameter Variance-Covariance matrix associated to \code{ARFIMA(p,d,q)} process.
Correspond to a generalization of the illustrated method for \code{ARFIMA(1,d,1)} model by Palma (2007). 
This method, is based in the explicit formula obtained by the derivatives of the parameters log-likelihood gradients.
}
\usage{
var.afm(phi = 0, theta = 0, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phi}{
a numeric vector containing the autoregressive polynomial parameters. By default, is 0.
}
  \item{theta}{
a numeric vector containing the mobile average polynomial parameters. By default, is 0.
}
  \item{n}{
sample total size.
}
}
\details{
We considering the distribution convergence \code{sqrt{n}*(widehat{alpha} - alpha) => N[0,Gamma^{-1}(alpha)]} so that \code{s.d(widehat{alpha}_i)=sqrt{(1/n)*(Gamma^{-1})_{ii}}} with

\code{Gamma=(1/4*pi)int_{-pi}^{pi}[nabla log f(lambda)]t([nabla log f(lambda)])dlambda}

where the gradients are obtained by the exact derivatives (see references). Finally, the integrals are obtained by \code{integrate} function.
}
\value{
  \item{G}{
the exact Variance-Covariance matrix.
}
 \item{SD}{
a vector with standard deviations associated to \code{d}, \code{phi_i} and \code{theta_j} parameters, with \code{i=1,...,p} and \code{j=1,...,q}.
}
}
\references{
Contreras-Reyes J. \& Palma W. (2011). Forecasting tools package for ARFIMA models. Preprint.

Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.
}
\author{
Javier Contreras-Reyes, Wilfredo Palma
}

\seealso{
\code{\link{arfima.whittle}}, \code{\link{summary.arfima}}, \code{\link{integrate}}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
S=arfima.whittle(y,nar=1,nma=1)
V=var.afm(phi=S$ar,theta=S$ma,n=length(y))
V$SD
summary(S)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{variance}
\keyword{arfima}
