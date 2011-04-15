\name{ir.arfima}
\alias{ir.arfima}

\title{
Impulse Response Function for ARFIMA models
}
\description{
Compute the Impulse Response (IR) associated to ARFIMA(p,d,q) model. Consider the \code{Normal} and 
\code{Asymptotic} (Hassler and Kokoszka, 2010) method(s). See below for details
}
\usage{
ir.arfima(h, model, plot.it = TRUE)
}

\arguments{
  \item{h}{a numeric value for number of lags}
  \item{model}{a fitted time-series model of ARFIMA class}
  \item{plot.it}{a \code{TRUE/FALSE} indicator for plotting the impulse response function}
}
\details{
For a fractional noise ARFIMA(0,d,0) process, we have the MA coefficients

\deqn{\psi_j=\frac{\Gamma(j+d)}{\Gamma(j+1)\Gamma(d)}.}

For a general ARFIMA(p,d,q) process we have the asymptotic approximation of MA coefficients

\deqn{\psi_j\sim\left(\frac{1+\sum_{i=1}^{q}\theta_i}{1-\sum_{i=1}^{p}\phi_i}\right)\frac{j^{d-1}}{\Gamma(d)}}

for a large \eqn{j}. The IR function is computed first via the \code{Normal} method

\deqn{R_j=\sum_{i=0}^{j}\psi_i\eta_{j-i},}

where \eqn{\eta_t=\Gamma(t+d)/(\Gamma(t+1)\Gamma(d))}. Then, the IR function is computed via the \code{Asymptotic} method

\deqn{R_j\sim\frac{j^{d-1}}{\Gamma(d)}\sum_{i=0}^{\infty}\psi_i,}

for a large \eqn{j}.
}
\value{
\item{h}{a numeric vector with values \code{1} to \code{h}}
\item{RE}{a numeric vector of \code{Normal} impulse response function}
\item{RA}{a numeric vector of \code{asymptotic} impulse response function}
\item{psi.j}{a numeric value of the asymptotic MA coefficients}
\item{int}{a numeric vector of coordinates of the intersection between \code{Normal} IR and \code{Asymptotic} IR 
functions}
}
\references{
Hassler, U., Kokoszka, P. (2010). Impulse Responses of Fractionally Integrated Processes with Long Memory. \emph{Econometric Theory}, 26, 6, 1885-1861.

Contreras-Reyes, J., Palma, W. (2011). Estimation, Diagnostic and Forecasting Tools for ARFIMA Models. Pre-print.
}
\author{
Javier E. Contreras-Reyes
}

\seealso{
\code{\link{psi.j}}, \code{\link{arfima.whittle}}
}
\examples{
data(MammothCreek)
y=MammothCreek-mean(MammothCreek)
mod <- arfima.whittle(y,nar=0,nma=0)
mod1 <- arfima.whittle(y,nar=1,nma=1)

ir.arfima(50, mod)
ir.arfima(50, mod1)
ir.arfima(100, mod)
ir.arfima(100, mod1)
ir.arfima(150, mod)
ir.arfima(150, mod1)
}
\keyword{asymptotic}
\keyword{arfima}
