tsdiag.arfima <-
function(object, gof.lag = 10, ...)
{
	res=residuals(object)
	pv=rep(NA,gof.lag)
	res.est=(res-mean(res))/sqrt(var(res))
	for(i in 1:gof.lag) pv[i]=Box.test(res,lag=i)$p.value

	par(mfrow=c(3,1))
	plot(res.est,type="h",main="Standardized Residuals",xlab="Time",ylab="")
	abline(h=0)
	acf(res,lag=gof.lag,main="ACF of Residuals",xlab="Lag",ylab="ACF")
	plot(pv,type="p",main="p values for Ljung-Box statistic",xlab="lag",ylab="p value",ylim=c(0,1))
	abline(h=0.05,lty=2,col="blue")
}

