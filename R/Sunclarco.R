# Project: Copula Package
# 
# Author: lucp8394
###############################################################################


#' @title Survival analysis for unbalanced clusters using archimedes copula's.
#'
#' @description Multivariate survival data occur in many different disciplines. 
#' If the correlation in the data is of interest by itself, two main modeling tools exist, the frailty model (Duchateau and Janssen, 2008) and the copula model.
#' The use of copula modeling has been restricted due to the lack of software on the one hand.     
#' On the other hand, copula models mainly deal theoretically with small clusters, all of equal size, i.e., containing the same number of subjects.    
#' The Sunclarco package can handle large unbalanced clusters, i.e., of varying size.
#' This allows the use of copula models to data sets that could previously not be handled, e.g., multicentre cancer clinical trials.
#' Furthermore, the Sunclarco package is flexible in terms of the baseline hazard (Weibull, piecewise exponential, unspecified (using partial likelihood) and in terms of the copula function (Clayton and Gumbel-Hougaard).
#' 
#' 
#' @author Leen Prenen
#' @author Ewoud De Troyer
#' @author Roel Braekers
#' @author Luc Duchateau
#' 
#' @references  Prenen L, Braekers R, Duchateau L (2017). Extending the Archimedean copula methodology to model multivariate survival data grouped in clusters of variable size. Journal of the Royal Statistical Society, 6, 1-24.
#' @references  Duchateau L, Janssen P. (2008). The frailty model. Spinger Verlag.
#'
#' @docType package
#' @name Sunclarco
NULL
