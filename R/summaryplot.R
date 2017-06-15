# Project: Copula Package
# 
# Author: Gebruiker
###############################################################################


## ADD SUMMARY, PRINT AND PLOT

#' @export
summary.Sunclarco <- function(object,...){
	cat("Execution Time:",object$info$runtime,"mins\n")
	cat("Copula:",object$info$parameters$copula,"\n")
	cat("Marginal Survival Distribution:",object$info$parameters$marginal,"\n")
	cat("Loglikelihood:",object$logllh,"\n\n")
	print(object$Parameters)
	cat("\n")
	temp <- matrix(object$Kendall_Tau,nrow=1,ncol=2)
	rownames(temp) <- "Kendall's Tau"
	colnames(temp) <- c("Estimate","StandardError")
	print(temp)
}


#' @export
print.Sunclarco <- function(x,...){
	summary(x)
}

