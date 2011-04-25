rho.arma <-
function(object=NULL, lag.max=NULL, ar=object$ar, ma=object$ma) 
{

n=length(object$data)
if (is.null(lag.max)) lag.max=floor(10 * (log10(n)))

a = spectrum.arfima(ar=ar, ma=ma, d=0, sd.innov=1)
b=NULL

for (j in 0:lag.max) {
aux.r=function(lambda) {
Re(a(lambda)*exp(-1i*lambda*j))
}
aux.i=function(lambda) {
Im(a(lambda)*exp(-1i*lambda*j))
}
b[j]=integrate(aux.r, -pi,pi)$value+1i*integrate(aux.i, -pi,pi)$value
}
b=Re(b)
ac=b/b[1]

return(ac)
}

