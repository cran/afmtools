gw.test <-
function(x,y,p,T,tau,method=c("HAC","NeweyWest","Andrews","LumleyHeagerty"), alternative=c("two.sided","less","greater"))
{

if (is.matrix(x) && ncol(x) > 2) 
stop("multivariate time series not allowed")
if (is.matrix(y) && ncol(y) > 2) 
stop("multivariate time series not allowed")
if (is.matrix(p) && ncol(p) > 2) 
stop("multivariate time series not allowed")

	# x: predicciones modelo 1
	# y: predicciones modelo 2
	# p: observaciones
	# T: sample total size
	# tau: horizonte de prediccion
	# method: if tau=1, method=NA. if tau>1, methods
	# alternative: "two.sided","less","greater"
	
	if(NCOL(x) > 1) stop("x is not a vector or univariate time series")     
	if(tau < 1) stop("Predictive Horizon must to be a positive integer")     
	if(length(x) != length(y)) stop("size of x and y difier")

	alternative <- match.arg(alternative)     
	DNAME <- deparse(substitute(x)) 

	l1=abs(x-p)
	l2=abs(y-p)
	dif=l1-l2
	q=length(dif)
	m=T-q	
	n=T-tau-m+1	
	delta=mean(dif)
	mod <- lm(dif~0+rep(1,q))

	if(tau==1){
		re=summary(mod)
            STATISTIC = re$coefficients[1,3]
		if (alternative == "two.sided") PVAL <- 2 * pnorm(-abs(STATISTIC))
		else if (alternative == "less") PVAL <- round(pnorm(STATISTIC),4)
     		else if (alternative == "greater") PVAL <- round(pnorm(STATISTIC, lower.tail = FALSE),4)     
		names(STATISTIC) <- "Normal Standad"
		METHOD <- "Standard Statistic Simple Regression Estimator" 
	}
	
	if(tau>1){
		
		if(method=="HAC"){ 
			METHOD <- "HAC Covariance matrix Estimation"
			ds=sqrt(vcovHAC(mod)[1,1])
		}
		if(method=="NeweyWest"){ 
			METHOD <- "Newey-West HAC Covariance matrix Estimation"
			ds=sqrt(NeweyWest(mod,tau)[1,1])
		}
		if(method=="LumleyHeagerty"){ 
			METHOD <- "Lumley HAC Covariance matrix Estimation"
			ds=sqrt(weave(mod)[1,1])
		}
		if(method=="Andrews"){ 
			METHOD <- "kernel-based HAC Covariance matrix Estimator"
			ds=sqrt(kernHAC(mod)[1,1])
		}
		#STATISTIC = sqrt(n)*delta/ds
		STATISTIC = delta/ds
		if (alternative == "two.sided") PVAL <- 2 * pnorm(-abs(STATISTIC))     
		else if (alternative == "less") PVAL <- pnorm(STATISTIC)
     		else if (alternative == "greater") PVAL <- pnorm(STATISTIC, lower.tail = FALSE)     
		names(STATISTIC) <- "Normal Standar"
	}
	structure(list(statistic = STATISTIC, alternative = alternative,p.value = PVAL, method = METHOD, data.name = DNAME))
}

