\name{arfima-methods}
\alias{plot.arfima}
\alias{summary.arfima}
\alias{residuals.arfima}
\alias{print.arfima}
\alias{tsdiag.arfima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for fitted ARFIMA models
}
\description{
\code{summary}, \code{print}, \code{residuals}, \code{tsdiag} and \code{plot} methods for class \code{arfima} model. 
A equivalent function of \code{summary} is provided in \code{afmtools} package called \code{summary.arfima}.
\code{tsdiag} its a generic diagnostic function which produces several plots of the residual from a fitted ARFIMA model.
}
\usage{
\method{plot}{arfima}(x, ...)
\method{residuals}{arfima}(object, ...)
\method{summary}{arfima}(object, ...)
\method{print}{arfima}(x, ...)
\method{tsdiag}{arfima}(object, gof.lag = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object, x}{ object of class \code{arfima}; usually a result of a call to \code{\link{arfima.whittle}}}
    \item{gof.lag}{the maximum number of lags for a Portmanteau goodness-of-fit test}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
\code{plot}{ 
produces 4 figure: 1) AR and MA roots of the model; 2) Sample ACF and theoretical ACF implied by the 
estimates; 3) Periodogram and theoretical spectrum implied by the estimates; 4) Sample ACF of Residuals
}
\code{summary}{ 
(and basically the same for \code{print}) gives a summary output in the \code{summary.lm} style - 
i.e. parameter estimates, standard errors, significance, etc.
}
\code{residuals}{ 
gives the residuals from the estimated ARFIMA model. This is not implemented directly via the 
\eqn{AR(\infty)} represenation of an ARFIMA(p,d,q) process, but using a \dQuote{trick}: first the 
original series is differenced with \eqn{\widehat{d}} using \code{diffseries} in the package 
\code{fracdiff}. Consequently an ARMA(p,q) should remain. Now instead of estimating the parameters 
again, an ARMA(p,q) model where ALL parameters are fixed to the Whittle estimates is \dQuote{estimated} 
with \code{\link[stats]{arima}} and then the residuals are obtained.
}
}
\value{
for \code{residuals} a vector of class \code{ts}; for \code{summary} and \code{print} the fitted 
ARFIMA model object. For \code{tsdiag} , produce plots of standardized residuals, autocorrelation 
function of the residuals, and the p-values of a Portmanteau test for all lags up to \code{gof.lag}.
}
\references{
Contreras-Reyes, J.E., Palma, W. (2012). Statistical Analysis of Autoregressive Fractionally Integrated Moving Average Models. arXiv:1208.1728v1.

Palma, W. (2007). \emph{Long Memory Time Series: Theory and Methods}. Wiley Series in Probability and Statistics. New Jersey.

Ljung, G.M., Box, G.E.P. (1978). On a measure of lack of fit in time series models. \emph{Biometrika} \bold{65}, 297--303.
}
\author{
Georg M. Goerg, Javier E. Contreras-Reyes
}
\seealso{
\code{\link{Box.test}}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
mod <- arfima.whittle(y,nar=1,nma=1)
res=residuals(mod)
summary(mod)
print(mod)
plot(res)
acf(res)
tsdiag(mod, gof.lag=12)
}
\keyword{arfima}

