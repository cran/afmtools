residuals.arfima <-
function(object, ...)
{

x=object$data
x=x-mean(x)
res=diffseries(x, object$d)

nar=object$nar
nma=object$nma

if (nar !=0 || nma !=0) {
if (nar !=0 && nma!=0) FIXED = c(object$ar, object$ma)
if (nar ==0 ) FIXED = object$ma
if (nma == 0) FIXED = object$ar

mod.arma=suppressWarnings(arima(res, order=c(object$nar, 0, object$nma), fixed = FIXED, include.mean=FALSE))
res=mod.arma$res
}
res[is.na(res)] = 0
ts(res)
}

