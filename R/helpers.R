# Project: Copula Package
# 
# Author: Gebruiker
###############################################################################


varphi         <- function(theta,t){(1+theta*t)^(-1/theta)}
varphi.inverse <- function(theta,t){(t^(-theta)-1)/theta}
varphi.prime   <- function(theta,t){-(1+theta*t)^(-1/theta-1)}



varphiGH         <- function(theta,t){exp(-t^theta)}
varphiGH.inverse <- function(theta,t){(-log(t))^(1/theta)}
varphiGH.prime   <- function(theta,t){-theta*exp(-t^theta)*t^(theta-1)}



coeffs <-function(alpha,p){
	coeff <- vector("list",length=p)
	coeff[[1]][1]<-1
	coeff[[2]][1]<-1-alpha
	coeff[[2]][2]<-1
	for (i in 3:p){
		coeff[[i]][1]<-prod(seq(1,i-1)-alpha)
		coeff[[i]][i]<-1
		for (j in 2:(i-1)){
			coeff[[i]][j]<-coeff[[i-1]][j-1]+coeff[[i-1]][j]*((i-1)-j*alpha)
		}
	}
	return(coeff)}

logdth.deriv_GumbHoug <- function(d,t,theta){
	K <- length(d)
	alpha <- theta
	delta <- theta
	gamma <- 0
	#varphi <- c()
	logderiv <- c()
	coeff <- coeffs(alpha,max(d)+3)
	for (i in 1:K){
		if (d[i]>=1){
			l <- 1:d[i]
			logterms <- log(coeff[[d[i]]])+l*log(delta)+(l*alpha-d[i])*log(gamma+t[i])
			logderiv[i] <- log(varphiGH(theta,t[i]))+max(logterms)+log(sum(exp(logterms-max(logterms))))} #op (-1)^d[i] na
		else #0'th derivative
		{logderiv[i]<-log(varphiGH(theta,t[i]))}
	}
	return(list(coeff=coeff,logderiv=logderiv))}



if(getRversion() >= "2.15.1"){
  globalVariables(c("temp_formula"))
  
}