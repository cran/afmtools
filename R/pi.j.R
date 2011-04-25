pi.j <-
function(h,ar=0,ma=0,d=0) 
{
	options(warn=-1)
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
    	ele=function(lambda) Phi.z(exp(-1i*lambda))/Theta.z(exp(-1i*lambda))
    	if (nma == 0) ele=function(lambda) Phi.z(exp(-1i*lambda))
	r=(h^(-d-1))*Re(ele(0))/gamma(-d) 
	if(h==0) r=1
	return(r) 
}

