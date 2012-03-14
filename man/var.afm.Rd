\name{var.afm}
\alias{var.afm}
\title{
Parameter Variance-Covariance matrix for ARFIMA models
}
\description{
Method to estimate the Parameter Variance-Covariance matrix associated to ARFIMA\eqn{(p,d,q)} process.
Correspond to a generalization of the illustrated method for ARFIMA\eqn{(1,d,1)} model by Palma (2007). 
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
We considering the distribution convergence 
\eqn{\sqrt{n}(\widehat{\alpha} - \alpha)\rightarrow N[0,\Gamma^{-1}(\alpha)]} so that 
\eqn{s.d(\widehat{\alpha}_i)=\sqrt{\frac{1}{n}\Gamma_{ii}^{-1}}} with

\deqn{\Gamma=\frac{1}{4\pi}\int_{-\pi}^{\pi}[\nabla \log f(\lambda)][\nabla \log f(\lambda)]'d\lambda}

where the gradients are obtained by the exact derivatives (See Contreras-Reyes J. & Palma W., 2011). 
Finally, the integrals are obtained by \code{\link{integrate}} function.
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
Contreras-Reyes, J., Palma, W. (2011). Estimation, Diagnostic and Forecasting Tools for ARFIMA Models. Pre-print.

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

# ARFIMA(1,d,1) particular case: (Palma, 2007; page 106 & 108).

phi1=-0.5
theta1=0.2
n1 = 100

var.arfima1d1 <- 
function(phi,theta,n)
{
V <- matrix(NA,3,3)

V[1,1]=pi^2/6
V[2,2]=1/(1-phi^2)
V[3,3]=1/(1-theta^2)
V[1,2]=V[2,1]=-log(1+phi)/phi
V[1,3]=V[3,1]=log(1+theta)/theta
V[2,3]=V[3,2]=-1/(1-phi*theta)
V1=(solve(V))/n

return(V1)
}

var.arfima1d1(phi = -phi1, theta = theta1, n = n1)
# note the conversion by definition: phi = -phi1

var.afm(phi = phi1, theta = theta1, n = n1)$G
}
\keyword{variance}
\keyword{estimation}
\keyword{arfima}