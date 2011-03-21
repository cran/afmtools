var.afm <-
function(phi=0,theta=0,n)
{
	if(phi[1]==0) warning("No Parameters for AR", "\n")
 	if(theta[1]==0) warning("No Parameters for MA", "\n")
	check.parameters.arfima(0,phi,theta)
	
	p<-length(phi)
	q<-length(theta)
	np=1+p+q	
	par1=c(1,-phi)
	par2=c(1,theta)

	grad.int <- function(lambda,x,pos) {
		if(pos==1){ g=-log(2*(1-cos(lambda)))
				 h=1}
		if(pos!=1){ if(x==0) { p=par1
				     	     sig=-2}
				if(x==1) { p=par2
				     	     sig=2}
			g=0 # numerador
			for (i in 1:length(p)){
				m=sig*p[i]*cos((pos-i)*lambda)
				g=m+g}
			h=sum(p^2) # denominador
			for (i in 2:length(p)){
				for (j in 1:(i-1)){
					r=2*p[i]*p[j]*cos((i-j)*lambda)
					h=r+h}}}
			g/h }

	prod.grad <- function(lambda,posA,posB,x.A,x.B){
		grad.int(lambda,x.A,posA)*grad.int(lambda,x.B,posB) }
	
	mp <- function(p,q){
		r=rep(1,np)
		r[2:(1+p)]=2:(1+p)
		r[(2+p):np]=2:(1+q)
		U=V=matrix(NA,np,np) 
		for(i in 1:np) for(j in i:np) {
			U[i,j]=U[j,i]=r[i]
			V[i,j]=V[j,i]=r[j] }
		list(U=U,V=V) }
 
	U=mp(p,q)$U
	V=mp(p,q)$V
	G<-matrix(0,np,np)

	for (i in 1:np) for (j in i:np){
		x.A=0
		if(i>(p+1) & j>(p+1)) x.A=1		     
		x.B=0
		if(i>(p+1) | j>(p+1)) x.B=1				
		G[i,j]=G[j,i]=integrate(prod.grad,lower=0,upper=pi,posA=U[i,j],posB=V[i,j],x.A=x.A,x.B=x.B)$value/(2*pi) }
	
	names=c(paste("d",sep=""),paste("ar",1:p,sep=""),paste("ma",1:q,sep=""))
	rownames(G)=names
	colnames(G)=names
	if(theta[1]==0) G=G[-(p+2),-(p+2)]
	if(phi[1]==0) G=G[-2,-2]	
	G=round(solve(G)/n,5)	
	SD=round(sqrt(diag(G)),5)
	list(G=G,SD=SD)
}

