smv.afm <-
function(object,lag.max=NULL,comp=TRUE)
{
	n=length(object$data)

	if(comp==TRUE){
		saux=1-seq(1,n-1,1)/n
		gam=rho.sowell(object=object,lag.max=(n-1),plot=FALSE)$av[-n]
		v=(2*saux%*%gam + gam[1])/n
	}
	else {
		cg=spectrum.arma(ar=-object$ar,ma=object$ma,sd.innov=object$sd.innov)
		cgamma=function(lambda) 2*cg(lambda)*gamma(1-2*object$d)*sin(pi*object$d)
		v=cgamma(0)*(n^(2*object$d-1))/((1+2*object$d)*object$d)	
	}
	return(as.numeric(round(v,4)))
}

