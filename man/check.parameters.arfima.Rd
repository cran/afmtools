\name{check.parameters.arfima}
\alias{check.parameters.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter check for ARFIMA models
}
\description{
Given the autoregressive and moving average and fractional coefficients, this function checks the
stationarity of the process.
}
\usage{
check.parameters.arfima(d = 0, ar = 0, ma = 0, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
long memory parameter
}
  \item{ar}{
autoregressive parameters
}
  \item{ma}{
moving average parameters
}
  \item{plot}{
should the output be plotted
}
}
\value{
	\item{d}{ 
\code{TRUE/FALSE} indicators for the fractional parameter
}
	\item{ar}{ 
\code{TRUE/FALSE} indicators for the autoregressive polynomial
}
	\item{ma}{ 
\code{TRUE/FALSE} indicators for the mobile average polynomial
}
	\item{model.OK}{ 
\code{TRUE/FALSE} indicators for the whole model
}
}
\references{
Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics, New Jersey.

Shumway, R., Stoffer, D. (2006). \emph{Time Series Analysis and Its Applications: With R Examples}, Springer. \url{http://www.stat.pitt.edu/stoffer/tsa2/index.html}.
}
\author{
Georg M. Goerg
}
\section{Side Effects}{
A graphical output is produced if \code{plot==TRUE}. This plot consists of the unit circle and the inverse AR and MA Roots.
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
mod <- arfima.whittle(y,nar=1,nma=1)
check.parameters.arfima(d=mod$d,ar=mod$ar,ma=mod$ma,plot=TRUE)

###### Some particular cases ######
# Invertible process
check.parameters.arfima(d=0.2,ar=-0.4,ma=0.7,plot=TRUE)
# No invertible process 
check.parameters.arfima(d=0.2,ar=-0.4,ma=c(0.7,-0.9),plot=TRUE)
# Causal process 
check.parameters.arfima(d=-0.5,ar=0.4,ma=-0.7,plot=TRUE) 
# No causal process
check.parameters.arfima(d=-0.5,ar=-1.4,ma=-0.7,plot=TRUE) 
# Fractional parameter in the(-1,0.5) interval
check.parameters.arfima(d=-0.7,ar=0.4,ma=0,plot=TRUE)
# Fractional parameter outside the (-1,0.5) interval 
check.parameters.arfima(d=0.5,ar=0.4,ma=0,plot=TRUE) 
}
\keyword{diagnostic}
\keyword{arfima}