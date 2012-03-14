\name{psi.j}
\alias{psi.j}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic Infinite MA expansion
}
\description{
Computes the asymptotic MA coefficients related to a ARFIMA process.
}
\usage{
psi.j(h, ar = 0, ma = 0, d = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{h}{
a numeric value related to the selected coefficient
}
  \item{ar}{
a numeric vector containing the autoregressive polynomial parameters
}
  \item{ma}{
a numeric vector containing the mobile average polynomial parameters
}
  \item{d}{
a numeric value of fractional parameter
}
}
\details{
Under the assumption that the roots of the polynomials \eqn{\Phi(B)} and \eqn{\Theta(B)} 
are outside the closed unit disk and \eqn{d \in (-1,0.5)}, the
\eqn{y_t \sim \mbox{ARFIMA}(p,d,q)} process is stationary, causal, and invertible. In this case 
we can write \eqn{y_t = \Psi(B)\varepsilon_t}. Then,
the MA(\eqn{\infty}) coefficients, \eqn{\psi_j}, associated to \eqn{\Psi(B)} polynomial are 
given by the following asymptotic relationship

\deqn{\psi_j\sim \frac{\Theta(1)}{\Phi(1)}\frac{j^{d-1}}{\Gamma(d)}}
}
\value{
a numeric value of the asymptotic MA coefficients
}
\references{
Kokoszka P. S. & Taqqu M. S. (1995). Fractional ARIMA with stable innovations. \emph{Stochastic Processes and Their Applications} 60, 19-47.

Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.
}
\author{
Javier Contreras-Reyes
}
\seealso{
\code{\link{pi.j}}, \code{\link{predict}}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
mod <- arfima.whittle(y,nar=2,nma=2)
h=1:100
pwd=psi.j(h,ar=mod$ar,ma=mod$ma,d=mod$d)
plot(h,pwd,type="l",ylab=" ")

#### Asymptotic convergence
h1=1:12
psi1=psi.j(h1,ar=-0.4,ma=0.7,d=-0.5)
h2=13:50
psi2=psi.j(h2,ar=-0.4,ma=0.7,d=-0.5)
h=1:50
plot(h,c(psi1,rep(NA,50-12)),type="l",ylab=" ")
lines(h,c(rep(NA,12),psi2),lty=2)
}
\keyword{asymptotic}
\keyword{forecasting}
\keyword{arfima}