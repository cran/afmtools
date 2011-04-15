\name{afmtools-package}
\alias{afmtools-package}
\alias{afmtools}
\docType{package}
\title{
Estimation, Diagnostic and Forecasting functions for ARFIMA models}
\description{
A collection of estimation, forecasting and diagnostic tools for autoregressive fractionally
integrated moving-average process (ARFIMA).
}
\details{
\tabular{ll}{
Package: \tab afmtools\cr
Type: \tab Package\cr
Version: \tab 0.1.8\cr
Date: \tab 2011-03-22\cr
License: \tab GPL >= 2\cr
LazyLoad: \tab yes\cr
}
}
\section{Functions}{
The package includes several functions. The following ones are those more relevant for practical use:

summary, plot, print, residuals and tsdiag options associated to \code{arfima} class object and residuals diagnostic. 
\code{\link{arfima.whittle}} and \code{\link{arfima.whittle.loglik}} for Whittle estimation (produce an \code{arfima} class object) 
and the log-likelihood function. \code{\link{gw.test}} is the forecasting tool. 
\code{\link{spectrum.arfima}}, \code{\link{rho.sowell}} and \code{\link{var.afm}} produce the spectrum density, autocovariance function 
and parameter variance for ARFIMA models, respectively.

It is suggested that the user starts by reading the documentation of (some of) these functions.
}
\section{Requirements}{
  R >= 2.6.0
  
  Packages \code{fracdiff}, \code{polynom}, \code{longmemo}, \code{sandwich} and \code{hypergeo}.
}
\author{
Javier E. Contreras-Reyes, Instituto de Fomento Pesquero, Valparaiso, Chile.
<Email: jecontrr@mat.puc.cl>.

Georg M. Goerg, Department of Statistics, Carnegie Mellon University. <Email: gmg@stat.cmu.edu>

Wilfredo Palma, Departamento de Estadistica, Facultad de Matematicas, Pontificia Universidad Catolica de Chile. <Email: wilfredo@mat.puc.cl>

Maintainer: Javier E. Contreras-Reyes <jecontrr@mat.puc.cl>
}

\section{Licence}{
This package and its documentation are usable under the terms of the \emph{GNU General Public License}, a copy of which is distributed with the
package. While the software is freely usable, it would be appreciated if a reference is inserted in publications or other work which makes use of it.
}

\references{
Brockwell, P., Davis, R. (1991). \emph{Time Series: Theory and Methods}. Springer. New York.

Contreras-Reyes, J., Idrovo, B. (2011). En busca de un modelo Benchmark Univariado para
predecir la tasa de desempleo de Chile. \emph{Cuadernos de Economia}, 30, 55, 105-125.

Contreras-Reyes, J.E., Palma, W. (2012). Statistical Analysis of Autoregressive Fractionally Integrated Moving Average Models. arXiv:1208.1728v1.

Giacomini, R., White, H. (2006). Tests of Conditional Predictive Ability. \emph{Econometrica} 74, 6.

Graybill, D.A. (1990). Pinus longaeva tree ring data. \emph{Mammoth Creek}, Utah, National Climatic Data Center.

Hassler, U., Kokoszka, P. (2010). Impulse Responses of Fractionally Integrated Processes with Long Memory. \emph{Econometric Theory}, 26, 6, 1885-1861.

Karamata, J. (1930). Sur un mode de croissance r\'eguliere des fonctions. Mathematica (Cluj) 4, 38-53.

Kokoszka, P.S., Taqqu M.S. (1995). Fractional ARIMA with stable innovations. \emph{Stochastic Processes and Their Applications} 60, 19-47.

Ljung, G.M., Box G.E.P. (1978). On a measure of lack of fit in time series models. \emph{Biometrika} 65, 297-303.

Palma, W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.

Palma, W., Olea, R. (2010). An efficient estimator for Gaussian locally stationary processes. \emph{The Annals of Statistics} 38, 2958-2997.

Shumway, R., Stoffer, D. (2006). \emph{Time Series Analysis and Its Applications: With R Examples}, Springer. \url{http://www.stat.pitt.edu/stoffer/tsa2/index.html}

Sowell, F. (1992). Maximum likelihood estimation of stationary univariate fractionally integrated time series models. \emph{Journal of Econometrics} 53, 165-188.

Whittle, P. (1951). \emph{Hypothesis Testing in Time Series Analysis}. Hafner, New York.
}
\keyword{arfima}
\keyword{package}
