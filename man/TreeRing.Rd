\name{TreeRing}
\alias{TreeRing}
\docType{data}
\title{
Tree Rings Time Series Data
}
\description{
Width measurements Time Series Data corresponding to Campito Tree Rings.
}
\usage{data(TreeRing)}
\format{
\code{TreeRing} is a data frame with 1164 observations.
}
\source{
Available at URL \url{http://www.mat.puc.cl/~wilfredo/english/data/whitemtn.htm}. Developed by Hipel K.W. & McLeod A.I. (1994)
}
\references{
Contreras-Reyes, J.E., Palma, W. (2012). Statistical Analysis of Autoregressive Fractionally Integrated Moving Average Models. arXiv:1208.1728v1.

Palma, W., Olea, R. (2010). An efficient estimator for Gaussian locally stationary processes. \emph{The Annals of Statistics} 38, 2958-2997.
}
\examples{
data(TreeRing)
str(TreeRing)
plot(TreeRing)

y=TreeRing-mean(TreeRing)
plot(y,type="l",xlab="Time (annual)",ylab="width measurements")
acf(y,lag=100,main="")
}
\keyword{datasets}
