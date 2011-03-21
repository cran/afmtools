\name{afmtools-package}
\alias{afmtools-package}
\alias{afmtools}
\docType{package}
\title{
Estimation, diagnostic and Forecasting functions for ARFIMA models}
\description{
A collection of estimation, forecasting and diagnostic tools for autoregressive fractionally
integrated moving-average process (ARFIMA).
}
\details{
\tabular{ll}{
Package: \tab afmtools\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2011-03-21\cr
License: \tab GPL >=2\cr
LazyLoad: \tab yes\cr
}
}
\section{Functions}{
The package includes several functions. The following ones are those more relevant for practical use:

summary, plot, print, residuals and tsdiag options associated to \code{arfima} class object and residuals diagnostic. 
\code{\link{arfima.whittle}} and \code{\link{arfima.whittle.loglik}} for Whittle estimation (produce an \code{arfima} class object) 
and the log-likelihood function. \code{\link{gw.test}} and \code{\link{pred.arfima}} are the forecasting tools. 
\code{\link{spectrum.arfima}}, \code{\link{rho.sowell}} and \code{\link{var.afm}} produce the spectrum density, autocovariance function 
and parameter variance for ARFIMA models, respectively.

It is suggested that the user starts by reading the documentation of (some of) these functions.
}
\section{Requirements}{
  R 2.12.2. Packages \code{fracdiff}, \code{polynom}, \code{longmemo}, \code{sandwich} and \code{hypergeo}.
}
\author{
Javier Contreras-Reyes, Seismological Service, Department of Geophysics, Universidad de Chile. <Email:  jcontr@dgf.uchile.cl>. 

Georg M. Goerg, Department of Statistics, Carnegie Mellon University. <Email: gmg@stat.cmu.edu>

Wilfredo Palma, Department of Statistics, Faculty of Mathematics, Pontificia Universidad Cat\'olica de Chile. <Email: wilfredo@mat.puc.cl>

Please send comments, error reports, etc. to the maintainer Javier Contreras-Reyes.
}

\section{Licence}{
This package and its documentation are usable under the terms of the "GNU General Public License", a copy of which is distributed with the
package. While the software is freely usable, it would be appreciated if a reference is inserted in publications or other work which makes use of it.
}

\references{
Bondon P. and Palma W. (2007). A class of antipersitent processes. \emph{Journal of Time Series Analysis} 28, 261-273.

Brockwell, P. and Davis, R. (1991). \emph{Time Series: Theory and Methods}. Springer. New York.

Contreras-Reyes J. and Palma W. (2011). Forecasting tools package for ARFIMA models. Preprint.

Giacomini R. and White H. (2006). Tests of Conditional Predictive Ability. \emph{Econometrica} 74, 6.

Graybill D. A. (1990). Pinus longaeva tree ring data. \emph{Mammoth Creek}, Utah, National Climatic Data Center.

Kokoszka P. S. and Taqqu M. S. (1995). Fractional ARIMA with stable innovations. \emph{Stochastic Processes and Their Applications} 60, 19-47.

Ljung G. M. and Box G. E. P. (1978). On a measure of lack of fit in time series models. \emph{Biometrika} 65, 297-303.

Palma W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.

Shumway, R. and Stoffer, D. (2006). \emph{Time Series Analysis and Its Applications: With R Examples}, Springer. \url{http://www.stat.pitt.edu/stoffer/tsa2/index.html}

Sowell F. (1992). Maximum likelihood estimation of stationary univariate fractionally integrated time series models. \emph{Journal of Econometrics} 53, 165-188.
}
\keyword{arfima}
\keyword{longmemory}
\keyword{package}
