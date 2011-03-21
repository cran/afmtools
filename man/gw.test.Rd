\name{gw.test}
\alias{gw.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Giacomini & White test of Predictive Ability
}
\description{
Forecasting evaluation method proposed by Giacomini and White (GW) to compare two vectors of
predictions provided by two time series models.
}
\usage{
gw.test(x, y, p, T, tau, method=c("HAC","NeweyWest","Andrews","LumleyHeagerty"), alternative=c("two.sided","less","greater"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric vector of data values corresponding to model 1
}
  \item{y}{
a numeric vector of data values corresponding to model 2
}
  \item{p}{
a numeric vector of data values corresponding to observations
}
  \item{T}{
sample total size
}
  \item{tau}{
forecasting horizon parameter. If \code{tau=1}, the \emph{Standard Statistic Simple Regression Estimator} method is used. If \code{tau>1}, the chosen \code{method} by the user is used
}
  \item{method}{
only if \code{tau>1} is chosen, \code{method} able to select between a set of Matrix Covariance Estimation methods, such as \code{HAC}, \code{NeweyWest}, \code{Andrews} and \code{LumleyHeagerty}
}
  \item{alternative}{
a character string specifying the alternative hypothesis, must be one of \code{two.sided} (default), \code{greater} or \code{less}
}
}
\value{
	\item{statistic}{ 
the value of the GW statistic
}
	\item{alternative}{ 
a character string describing the alternative hypothesis
}
	\item{p.value}{ 
the p-value for the test
}
	\item{method}{ 
a character string indicating what type of Matrix Covariance Estimation method was performed
}
	\item{data.name}{ 
a character string giving the name(s) of the data
}
}
\references{
Giacomini R. \& White H. (2006). Tests of Conditional Predictive Ability. \emph{Econometrica} \bold{74}, 6.
}
\author{
Javier Contreras-Reyes
}
\seealso{
\code{pred.arfima}, \code{predict}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
r=30
obs=y[1:(length(y)-r)]
nobs=1:r
y.real=y[(length(y)-r+1):length(y)]
mod.arima <- arima(obs, order = c(1,1,1))
mod.arfima <- arfima.whittle(obs,nar=1,nma=1)
P=matrix(NA,r,2) # matrix of predictions

for(i in 1:r) {
	p.arma <- predict(mod.arima,ahead=1)
	p.arfima <- pred.arfima(object=mod.arfima,ahead=1)
	P[i,]=c(p.arma$pred, p.arfima$pred)
	obs=c(obs,y.real[i])
	mod.arima <- arima(obs, order = c(1,1,1))
	mod.arfima <- arfima.whittle(obs,nar=1,nma=1)
}

plot(nobs,abs(P[,1]-y.real),type="l",col="red",main="Absolute Forecasting Errors",
ylab=expression(abs(y-hat(y))))
lines(nobs,abs(P[,2]-y.real),type="l",col="blue")
legend(0,0.9,c("ARIMA(1,1,1)","ARFIMA(1,d,1)"),col=c("red","blue"),lty=c(1,1))

### Giacomini-White test
gw.test(x=P[,1],y=P[,2],p=y.real,T=length(y),tau=1,method="HAC",alternative="less")
gw.test(x=P[,2],y=P[,1],p=y.real,T=length(y),tau=1,method="Andrews",alternative="greater")
}
\keyword{test}

