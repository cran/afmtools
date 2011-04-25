pred.arfima <-
function(object, ahead=1, crt=50)
{
	T=length(object$data)

	if(T < crt){
		G <- matrix(NA,T,T)
		pred <- rep(0,ahead)
		pred.err <- rep(NA,ahead)
		b.f <- matrix(NA,2,ahead)
		gam.aux=rho.sowell(object,lag.max=T+ahead-1,plot=FALSE)$av
	
		# Exact Autocovariance Gamma Matrix:
		for(i in 1:T) for(j in 1:T) G[i,j]=gam.aux[abs(i-j)+1] 

		for(h in 1:ahead){
			gam <- gam.aux[(h+1):(T+h)]
			phi=solve(G)%*%gam
			for (i in 1:T) pred[h] = pred[h] + phi[i]*object$data[T+1-i]
			sum=0
			for(j in 0:(h-1)) sum = sum + psi.j(j)
			pred.err[h]= object$sd.innov + sum*object$d*tan(pi*object$d)/(pi*T) 	
		}
	}

	else {
		pred <- rep(0,ahead)
		pred.err <- rep(NA,ahead)
		b.f <- matrix(NA,2,ahead)
		sum=0

		for(h in 1:ahead){ 
			for (j in 0:(T-1)) for(i in 0:(h-1)){
				pred[h] = pred[h] + psi.j(h=i,ar=object$ar,ma=object$ma,d=object$d)*pi.j(h=j+h-i,ar=object$ar,ma=object$ma,d=object$d)*object$data[T-j]
			}
			for(k in 0:(h-1)) sum = sum + psi.j(h=k,ar=object$ar,ma=object$ma,d=object$d)^2
			pred.err[h]= object$sd.innov*sum
		}	
	}

	list(pred=pred,pred.err=pred.err)
}

