arfima.whittle.loglik <-
function(theta, series, nar, nma, fixed)   	
{ 

theta[fixed]=0

d=theta[1]
phi = NA
if (nar >0) phi=theta[2:(nar+1)]
theta_ma=NA
if (nma >0) theta_ma=theta[(nar+2):length(theta)]

# Parameter (zero) restrictions

#if (!is.na(phi[1])) phi[fixed$ar]=0
#if (!is.na(theta_ma[1])) theta_ma[fixed$ma]=0
#if (!is.na(fixed$d)) d=fixed$d

M=check.parameters.arfima(d=d, ar=phi, ma=theta_ma)
if (M$Total.OK) {
# 	
#   Whittle Loglikelihood 	
# 	

if (any(is.na(series))) 
stop("missing values not allowed in time series")
if (is.matrix(series) && ncol(series) > 2) 
stop("multivariate time series not allowed")

 series <- series - mean(series)  
 n = length(series) 	
 # Periodogram starting with frequency omega_1 ( > 0 )
 P = per.arfima(series)
 a = P$spec
 w = P$freq	
 #   	
 #Spectral Density:   	
 #   	
 S.arfima=spectrum.arfima(d=d, ar=phi, ma=theta_ma)
 b <- S.arfima(w)  	
 #   	
 #  Calculate sigma^2   	
 #   	
 sigma2 <- (2 * sum(a/b))/n   	
 #   	
 #  Whittle Log-likelihood   	
 #   	
 loglik <- 2 * pi * (sum(log(b)) + sum(a/b)/sigma2)   
 L=loglik/n + pi * log(sigma2)
}
else {L = 10^10; sigma2=0}

return(list(L=L, sigma=sqrt(sigma2)))   	
}

