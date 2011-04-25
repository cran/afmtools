spectrum.arma <-
function(ar=0,ma=0,sd.innov=1)
{
    	if (is.na(ar[1])) ar=0
	if (is.na(ma[1])) ma=0
	nar=length(ar)
	nma=length(ma)
	if (nar == 1 && ar ==0) nar=0
	if (nma == 1 && ma ==0) nma=0
	M=check.parameters.arfima(d=0, ar=ar, ma=ma)
	if (!M$Total.OK) cat("WARNING: Model is not OK.")
 	ar.poly <- c(1, -ar)
 	z.ar <- polyroot(ar.poly)
 	ma.poly <- c(1, ma)
 	z.ma <- polyroot(ma.poly)
	Phi.z=as.function(polynomial(coef = ar.poly)) 
    	Theta.z=as.function(polynomial(coef = ma.poly))
    	k=function(lambda) Theta.z(exp(-1i*lambda))/Phi.z(exp(-1i*lambda))
    	if (nma == 0) k=function(lambda) 1/Phi.z(exp(-1i*lambda))
    	spec=function(lambda) (Mod(k(lambda))^2)*sd.innov^2/(2*pi) 
	invisible(spec) 
}

