spectrum.arfima <-
function(d=0, ar=0,ma=0, sd.innov=1)
{
	if (d >= 0.5) warning("WARNING: Process is not stationary.", "\n")
	if (d < -1) warning("WARNING: Process is not invertible.", "\n")

	S.arma=spectrum.arma(ar=ar,ma=ma,sd.innov=sd.innov)
	S.arfima=function(lambda) S.arma(lambda)*(2*sin(lambda/2))^(-2*d)
	invisible(S.arfima)
}

