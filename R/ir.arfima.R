ir.arfima <-
function(h, model, plot.it=TRUE) {

ar=round(model$ar,3)
ma=round(model$ma,3)
d=round(model$d,3)
Ej <- Aj <- rep(NA, h)
eta <- function(j, d) gamma(j+d)/(gamma(j+1)*gamma(d))

if(ar!=0 | ma!=0){
psjE = c(1, psi.j(c(1:h), ar, ma, d))
for(i in 1:h) Ej[i] = sum(psjE[1:(i+1)] * rev(eta(c(0:i), d)))
}

else {
psjE = eta(0:h,d)
for(i in 1:h) Ej[i] = sum(eta(0:i,d) * rev(eta(c(0:i),d)))
}

S = sum(psjE)
for(i in 1:h) Aj[i] = S*i^(d-1)/gamma(d)

if(plot.it) {
k=max(Ej, Aj, na.rm = FALSE)
p=length(ar); if(ar==0 & p==1) {p=0}
q=length(ma); if(ma==0 & q==1) {q=0}
plot(c(1:h), Ej, type="l", xlab="j", ylab=expression(R[j]), 
main=paste("Impulse Response: h=", h, 
", ARFIMA(", p, ",", d, ",", q, ")", sep = ""), ylim=c(0,k))
lines(c(1:h), Aj, lty=2)
legend(round(h/2,0), k, c("Normal","Asymptotic"), lty=c(1,2))
error = abs(Ej-Aj)
m = min(error, na.rm=TRUE)
pos=order(error)[1]
#abline(h=Aj[error==m],col="blue",lty=1)
#abline(v=pos,col="blue",lty=1)
g = round(Aj[error==m],3)
#text(pos+12, Aj[error==m]+0.02, paste("(", pos, ",", g, ")"))
return(list(h=c(1:h), RE=Ej, RA=Aj, psi.j=psjE, int=c(pos, g)))
}

return(list(h=c(1:h), RE=Ej, RA=Aj, psi.j=psjE))
}

