
rho.sowell <-
function(object=NULL, lag.max=NULL, ar=object$ar, ma=object$ma, 
d=object$d, sd.innov=object$sd.innov, plot=TRUE) {

	if(length(ar)==1 & ar[1]==0) warning("No Parameters for AR", "\n")
 	if(length(ma)==1 & ma[1]==0) warning("No Parameters for MA", "\n")
	check.parameters.arfima(0,ar,ma)
	n=length(object$data)
	if (is.null(lag.max)) lag.max=floor(10*(log10(n)))

	p=length(ar)
	q=length(ma)	
	theta=c(1,ma)
	rho=1/polyroot(c(1,-ar))
	
	C = function(h,r){
		gamma0=gamma(1-2*d)*gamma(h+d)/(gamma(1-d)*gamma(d)*gamma(1+h-d))
		gamma0*(r^(2*p)*hypergeo(d+h,1,1-d+h,r)+ hypergeo(d-h,1,1-d-h,r)-1) 
	}

	chi = function(j) {
		aux1=1	
		for(i in 1:p) aux1=aux1*(1-rho[i]*rho[j])
		if(p>1){
			ar.aux=rho[-j]
			aux2=1
			for(i in 1:(p-1)) aux2=aux2*(rho[j]-ar.aux[i]) 		
		}
		if(p==1) aux2=1
		(rho[j]*aux1*aux2)^1 
	}
	
	psi.ma = function(l) {
		k1=max(0,l)
		k2=min(q,q+l)
		sum=0
		for(i in k1:k2) sum=sum+theta[i+1]*theta[i+1-l]
		sum 
	}

	b=rep(NA,(lag.max+1))
	
	for (h in 1:length(b)) {
		if(h<=50){
		q.aux=seq(-q,q,1)
		b[h]=0
		for (i in 1:length(q.aux)){
			for (j in 1:p){ 
				if(length(ar)==1 & ar[1]==0){
					chi.aux=1
					C.aux=gamma(1-2*d)*gamma(h+d-q.aux[i])/(gamma(1-d)*gamma(d)*gamma(1+q.aux[i]-d-h))
				}
				if(ar[1]!=0){ 
					chi.aux=chi(j)
					C.aux=C(p+q.aux[i]-(h-1),rho[j])
				} 
				if(length(ma)==1 & ma[1]==0) psi.aux=1
				if(ma[1]!=0) psi.aux=psi.ma(q.aux[i])
				b[h] = b[h] + sd.innov*psi.aux*chi.aux*C.aux
			}
		}
		}
		if(h>50) b[h]=(abs(h)^(2*d-1))*sd.innov*((1+ma)^2)*gamma(1-2*d)*sin(pi*d)/(pi*(1+ar)^2)
	}
	av=round(Re(b),4)
	ac=round(av/av[1],4)
	if (plot) {
	  plot(c(1:lag.max-1), ac[1:12],type="h",main="Sowell ACF",ylab="ACF",xlab="Lag", ylim=c(-2/sqrt(n),1))
        #abline(h=c(0,-2/sqrt(n),2/sqrt(n)),col=c(1,4,4),lty=c(1,2,2))
		abline(h=c(0),col=c(1),lty=c(1))

    	}
	list(av=av, ac=ac)
}