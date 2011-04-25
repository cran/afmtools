check.parameters.arfima <-
function(d=0, ar=0, ma=0, plot=FALSE)
{

if (is.na(ma[1])) ma=0
if (is.na(ar[1])) ar=0

nar=length(ar)
nma=length(ma)

if (nar == 1 && ar ==0) nar=0
if (nma == 1 && ma ==0) nma=0

model=NULL
model$d.OK=TRUE
model$ar.OK=TRUE
model$ma.OK=TRUE
model$Total.OK=TRUE

if (d <=-1 || d >= 0.5) {
warning("WARNING: Parameter d not adequate.", "\n")  
model$d.OK=FALSE
model$Total.OK=FALSE
}
    # check causality
     ar.poly <- c(1, -ar)
     z.ar <- polyroot(ar.poly)
    
	if(any(abs(z.ar) <= 1)) {
warning("WARNING: Model Not Causal", "\n")
model$ar.OK=FALSE
model$Total.OK=FALSE
}    # check invertibility

     ma.poly <- c(1, ma)
     z.ma <- polyroot(ma.poly)
    
	if(any(abs(z.ma) <= 1)) {
warning("WARNING: Model Not Invertible", "\n")
model$ma.OK=FALSE
model$Total.OK=FALSE
}
#     if(any(abs(z.ma) <= 1) || any(abs(z.ar) <= 1) ) stop("Try Again")
    #

    # check (near) parameter redundancy [i.e. are any roots (approximately) equal]  
      
   if ( ((ar == 0 & nar == 1) || (ma == 0 & nma ==1))) {
	 for (i in 1:nar) {
       if ( (ar == 0 & nar == 1) || (ma == 0 & nma ==1) ) break
       if(any(abs(z.ar[i]-z.ma[1:nma]) < 1e-03)) {warning("WARNING: Parameter Redundancy", "\n"); break}
       }
}
if (plot) {
plot(1/z.ar, ylim=range(c(-1,1, 1/abs(z.ar), 1/abs(z.ma))), xlim=range(c(-1,1, 1/abs(z.ar), 1/abs(z.ma))), lwd=2, main="Inverse AR and MA Roots")
points(1/z.ma, pch=19, lwd=2)
abline(h=0)
abline(v=0)

legend("topleft", c(paste("AR Roots:", nar) , paste("MA Roots:", nma)), pch=c(1,19), col=1)

aux=function(lambda) {exp(1i*lambda)}
lines(aux((0:100)/100*2*pi), lwd=2)
}
invisible(model)
}

