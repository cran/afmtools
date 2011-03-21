plot.arfima <-
function(x, ...) 
{

n=length(x$data)
op=par(no.readonly=TRUE)

B=per.arfima(x$data)
acov.fn=ckARMA0(100, H = x$d+0.5)
corr.fn=acov.fn/acov.fn[1]
corr.fn.mirror=c(rev(corr.fn[-1]), corr.fn)

corr.arma=rho.arma(x)
#corr.arma.mirror=c(rev(corr.arma[-1]), corr.arma)

#corr=convolve(corr.fn.mirror, corr.arma.mirror)
#we=convolve(corr.fn.d, rev(corr.arma), type = "o")

par(mfrow=c(2,2), mar=c(2,2,3,1))
check.parameters.arfima(d=x$d, ar=x$ar, ma=x$ma, plot=TRUE)

acf(x$data, main="ACF")
if (x$d == 0) lines(0:(length(corr.arma)-1), corr.arma, lwd=2, col=2)
if (x$nar ==0 && x$nma==0) lines(0:(length(corr.fn)-1), corr.fn, lwd=2, col=2)

plot(B$freq, B$spec, log="y", main="Theoretical vs. Empirical Spectrum \n Logarithmic scale")
plot(x$spectrum, 0, pi, add=TRUE, col=2, lwd=2)

res=residuals(x)
acf(res, main="ACF of Residuals")

par(op)
}
