\name{MammothCreek}
\alias{MammothCreek}
\docType{data}
\title{
Tree Ring Time Series Data
}
\description{
Centered annual pinus longaeva tree ring time series width measurements at Mammoth Creek, Utah, from 0 A.D to 1989 A.D.
}
\usage{data(MammothCreek)}
\format{
\code{MammothCreek} is a data frame with 1990 observations.
}
\source{
Available at the National Climatic Data Center, URL \url{http://www.ncdc.noaa.gov/paleo/metadata/noaa-tree-3376.html}
}
\references{
Graybill D. A. (1990). Pinus longaeva tree ring data. \emph{Mammoth Creek}, Utah, National Climatic Data Center.

Palma W. & Olea R. (2010). An efficient estimator for Gaussian locally stationary processes. \emph{The Annals of Statistics} \bold{38}, 2958--2997.
}
\examples{
data(MammothCreek)
str(MammothCreek)
plot(MammothCreek)

y=MammothCreek-mean(MammothCreek)
plot(y,type="l",xlab="Time (annual)",ylab="width measurements")
acf(y,lag=100,main="")
}
\keyword{datasets}
