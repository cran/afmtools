\name{pred.arfima}
\alias{pred.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multistep Ahead Predictor for ARFIMA models
}
\description{
h-step ahead forecast linear method for a finite stretch of data. This method consider two ways to compute the predictions. First, consider
a finite past and a asymptotic converging standard deviation associated to a linear predictor founded by Bondon and Palma for stationary time series.
And second, assume a infinite past with infinite ar and ma expansions mentioned by Kokoszka & Taqqu through \code{\link{psi.j}} and \code{\link{pi.j}} functions.
}
\usage{
pred.arfima(object, ahead = 1, crt = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a fitted time-series model of ARFIMA class
}
  \item{ahead}{
a numeric value with time steps ahead to predict. By default, \code{ahead}=1
}
  \item{crt}{
a numeric value of criteria to select the finite or infinite method. If \code{T} is the length of the inserted data from \code{object}, then; if
\code{T<crt} the finite past method is selected and, in the other case, the infinite past method is used. By default, \code{crt}=50
}
}
\details{
For a small size of number of observations \eqn{T<\code{crt}}, \code{\link{pred.arfima}} 
compute the linear predictor

\deqn{\widehat{y}_{t}(h)=\sum_{i=0}^{t} \phi_{t,i+1}(h)y_{t-i+1}}

where \eqn{\phi_t(h)=[\phi_{t1}(h),\ldots,\phi_{tt}(h)]'} satisfies 
\eqn{\Gamma_t\phi_t(h)=\gamma_t(h)} with \eqn{\gamma_t(h)=[\gamma(h),\ldots,\gamma(t+h-1)]'}. 
The \eqn{\gamma_t(h)} terms are computed by \code{\link{rho.sowell}} function. In the case 
\eqn{T\geq\code{crt}}, \code{\link{pred.arfima}} compute the linear predictor

\deqn{\widehat{y}_{t}(h)=\sum_{j=0}^{t} \sum_{i=0}^{h-1}\psi_i\pi_{j+h-i}y_{t-j}}

with \eqn{\psi_i} and \eqn{\pi_j} are computed by \code{\link{psi.j}} and 
\code{\link{pi.j}} functions.
}
\value{
	\item{pred}{
predictions
}
	\item{pred.err}{
estimated standard errors
}
}
\references{
Bondon P. & Palma W. (2007). A class of antipersitent processes. \emph{Journal of Time Series Analysis} 28, 261-273.

Karamata J. (1930). Sur un mode de croissance r\'eguliere des fonctions. Mathematica (Cluj) 4, 38-53.

Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.
}
\author{
Javier Contreras-Reyes, Wilfredo Palma
}
\seealso{
\code{\link{predict}}, \code{\link{rho.sowell}}, \code{\link{psi.j}}, \code{\link{pi.j}}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
mod1 <- arfima.whittle(y,nar=0,nma=0)
pred.arfima(object=mod1,ahead=1)

x=y[(length(y)-300):length(y)]
mod2 <- arfima.whittle(x,nar=1,nma=1)
pred.arfima(object=mod2,ahead=3) # increase the prediction error

m=15
A <- rep(NA,m)
z=y[(length(y)-100):length(y)]
t=length(z)

for(i in 1:m) {
	mod <- arfima.whittle(z,nar=1,nma=1)
	A[i]=pred.arfima(object=mod,ahead=1)$pred
	z=c(z,A[i])
}

plot(c(z[1:t],rep(NA,m)),type="l")
lines(c(rep(NA,t),z[(t+1):length(z)]),type="l",col="red")

####### Prediction Error bands #######

y=MammothCreek-mean(MammothCreek)
r=30 # number of predictions
conf.level = 95
obs=y[1:(length(y)-r)]
nobs=1:r
y.real=y[(length(y)-r+1):length(y)]
mod.arfima <- arfima.whittle(obs,nar=1,nma=1)
P=rep(NA,r) # predictions
R=matrix(NA,r,2) # forecasting bands

for(i in 1:r) {
	p.arfima <- pred.arfima(object=mod.arfima,ahead=1)
	P[i]=p.arfima$pred
	se = qnorm(1-(1-conf.level/100)/2)*sqrt(p.arfima$pred.err)
	R[i,] = p.arfima$pred + c(se, -se)
	obs=c(obs,y.real[i])
	mod.arfima <- arfima.whittle(obs,nar=1,nma=1)
}

plot(nobs,y.real,type="l",main="Tree Ring predictions", 
ylab="width measurements", ylim=c(-1.4,1.2))
lines(nobs,P,type="l",col="blue")
lines(nobs,R[,1],type="l",col="red",lty=2)
lines(nobs,R[,2],type="l",col="red",lty=2)
legend(1,-0.5,c("observed", "predicted","f.b"),
col=c("black","blue","red"),lty=c(1,1,2))

####### For other examples, see 'gw.test' function
}
\keyword{forecasting}

