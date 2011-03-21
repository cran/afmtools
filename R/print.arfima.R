print.arfima <-
function(x, ...) 
{
    cat("Call: ")
    print(x$call)
    cat("Estimation method: ")
    cat(x$method)
    cat("\n")
    cat("\n Parameter estimates:\n")
    printCoefmat(summary(x)$coefmat, signif.stars = TRUE)
    cat("\n Estimated standard deviation of residuals (ratio of periodogram and spectrum):\n")
    cat(x$sd.innov)
    cat("\n")
}
