\name{pi.j}
\alias{pi.j}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Asymptotic Infinite AR expansion
}
\description{
Computes the asymptotic AR coefficients related to a ARFIMA process.
}
\usage{
pi.j(h, ar = 0, ma = 0, d = 0)
}
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
we can write \eqn{\varepsilon_t = \Pi(B)y_t}. Then,
the AR(\eqn{\infty}) coefficients, \eqn{\pi_j}, associated to \eqn{\Pi(B)} polynomial are 
given by the following asymptotic relationship

\deqn{\pi_j\sim \frac{\Phi(1)}{\Theta(1)}\frac{j^{-d-1}}{\Gamma(-d)}}


}
\value{
a numeric value of the asymptotic AR coefficient
}
\references{
Kokoszka P. S. and Taqqu M. S. (1995). Fractional ARIMA with stable innovations. \emph{Stochastic Processes and Their Applications} 60, 19-47.

Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.
}
\author{
Javier E. Contreras-Reyes
}
\seealso{
\code{\link{psi.j}}, \code{\link{predict}}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
mod <- arfima.whittle(y,nar=2,nma=1)
h=1:100
pwd=pi.j(h,ar=mod$ar,ma=mod$ma,d=mod$d)
plot(h,pwd,type="l",ylab=" ")

#### Asymptotic convergence
h1=1:24
psi1=psi.j(h1,ar=0,ma=0.2,d=0.3)
h2=25:70
psi2=psi.j(h2,ar=0,ma=0.2,d=0.3)
h=1:70
plot(h,c(psi1,rep(NA,70-24)),type="l",ylim=c(0,max(psi1)),ylab=" ")
lines(h,c(rep(NA,24),psi2),lty=3)
}
\keyword{asymptotic}
\keyword{forecasting}
\keyword{arfima}