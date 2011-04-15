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
gw.test(x, y, p, T, tau, 
        method=c("HAC","NeweyWest","Andrews","LumleyHeagerty"), 
        alternative=c("two.sided","less","greater"))
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
Giacomini R. & White H. (2006). Tests of Conditional Predictive Ability. \emph{Econometrica} \bold{74}, 6.

Contreras-Reyes, J., Idrovo, B. (2011). En busca de un modelo Benchmark Univariado para
predecir la tasa de desempleo de Chile. \emph{Cuadernos de Economia}, 30, 55, 105-125.
}
\author{
Javier E. Contreras-Reyes
}
\seealso{
\code{predict}
}
\examples{
r = 100
s = 30
y = arima.sim(n = r, list(ar = c(0.8, -0.4), 
ma = c(-0.2, 0.3)),  sd = sqrt(0.18))
y.real = y[(length(y)-s+1):length(y)]
obs = y[1:(r-s)]
mod.arma <- arima(obs, order = c(1,0,1))
mod.arima <- arima(obs, order = c(1,1,1))
P = matrix(NA, s, 2) # matrix of predictions

for(i in 1:s) {
	p.arma <- predict(mod.arma, ahead = 1)
	p.arima <- predict(mod.arima, ahead = 1)
	P[i,] = c(p.arma$pred, p.arima$pred)
	obs = c(obs, y.real[i])
	mod.arma <- arima(obs, order = c(1,0,1))
	mod.arima <- arima(obs, order = c(1,1,1))
}

gw.test(x = P[,1], y = P[,2], p = y.real, T = length(y), 
tau = 1, method = "HAC", alternative = "less")
}
\keyword{test}
\keyword{forecasting}
