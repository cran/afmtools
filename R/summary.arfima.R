summary.arfima <-
function(object, ...) 
{

if (object$nma == 0 && object$nar > 0) {
Gamma_0=CetaARIMA(eta = c(H = object$d, phi=object$ar), m = 10^5,  p = object$nar, q = 0)
R.names=c("d", paste("phi", 1:object$nar))
}
if (object$nar ==0 && object$nma > 0) {
Gamma_0=CetaARIMA(eta = c(H = object$d, psi=object$ma), m = 10^5, p=0, q = object$nma)
R.names=c("d", paste("theta", 1:object$nma))
}
if (object$nar==0 && object$nma==0) {
Gamma_0=CetaARIMA(eta = c(H = object$d), m = 10^5,  p = 0, q = 0)
R.names=c("d")
}
if (object$nar >0 && object$nma >0) {
Gamma_0=CetaARIMA(eta = c(H = object$d, phi=object$ar, psi=object$ma), m = 10^5,  p = object$nar, 
q = object$nma)
R.names=c("d", paste("phi", 1:object$nar), paste("theta", 1:object$nma))
}
	theta=object$par
	Hess = solve(Gamma_0)
	se = diag(Hess)/sqrt(length(object$data))
    	tval = theta/se
    	TAB = cbind(Esimate = theta, StdErr = se, t.value = tval, 
        	p.value = 2 * (1 - pnorm(abs(tval))))
    	dimnames(TAB) = list(names(tval), c(" Estimate", " Std. Error", 
        			" t value", "Pr(>|t|)"))
    	rownames(TAB)=R.names
    	OUT = NULL
    	OUT$call = object$call
    	OUT$coefmat = TAB
    	OUT$sd.innov = object$sd.innov
    	OUT$method=object$method
    	class(OUT) = c("summary.arfima")
 	OUT
}

