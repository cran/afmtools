per.arfima <-
function(z) 
{
	if (any(is.na(z))) 
	stop("missing values not allowed in time series")
	if (is.matrix(z) && ncol(z) > 2) 
	stop("multivariate time series not allowed")

	n = length(z)
 	OUT = NULL
 	OUT$spec=(Mod(fft(z))^2/(2 * pi * n))[2:(n%/%2 + 1)]
 	m = length(OUT$spec)
 	OUT$freq = (2 * pi * (1:m))/n 
 	invisible(OUT)
}

