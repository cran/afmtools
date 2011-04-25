arfima.whittle <-
function(series, nar=0, nma=0, fixed=NA)
{     	
  # Computes Whittle MLE for ARFIMA(p,d,q) models     	

if (any(is.na(series))) 
stop("missing values not allowed in time series")
if (is.matrix(series) && ncol(series) > 2) 
stop("multivariate time series not allowed")
	
d.ml=fracdiff(series)$d
res.arma=diffseries(series, d.ml)

mod.arma=arima(res.arma, order=c(nar, 0, nma))

start.p=c(d.ml, mod.arma$coef[-length(mod.arma$coef)])
names(start.p)[1] = "d"
start.p[fixed]=0

aux=function(theta, series, nar, nma, fixed) {
arfima.whittle.loglik(theta, series, nar, nma, fixed)$L
}

fit=suppressWarnings(nlm(f = aux, p = start.p, series, nar, nma, fixed))
param=fit$estimate

OUT=NULL
OUT$call = match.call()
OUT$method = c("Whittle")
OUT$data = series
OUT$par = param
OUT$loglik=fit$minimum
OUT$d=param[1]
OUT$ar=param[2:(nar+1)]
OUT$ma=param[(nar+2):length(param)]
if (nar ==0) OUT$ar=0
if (nma ==0) OUT$ma=0
OUT$nma = nma
OUT$nar = nar
OUT$fixed= fixed
OUT$sd.innov = arfima.whittle.loglik(param, series, nar, nma)$sigma
OUT$model=paste("ARFIMA(", nar, ",d,", nma, ")", sep="")

OUT$spectrum = spectrum.arfima(ar=OUT$ar, ma=OUT$ma, d=OUT$d, sd.innov=OUT$sd.innov)
OUT$method=c("Whittle")
class(OUT) = c("arfima")
OUT$call=match.call()
invisible(OUT)
}

