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





factor2cont <- function(data,factorcovariates,baselevels){
  
  if(!is.null(baselevels)){
    if(!(names(baselevels)%in%colnames(data))){stop("baselevels: one of the variables not in data.")}
    if(!(names(baselevels)%in%factorcovariates)){stop("baselevels: one of the variables not a factor.")}
    incorrect_base <- which(sapply(1:length(baselevels),FUN=function(i){  !(baselevels[i]%in%levels(data[,names(baselevels)[i]]))  }))
    if(length(incorrect_base)>0){stop(paste0("Error in baselevels parameter: ",paste0(paste0("\"",baselevels[incorrect_base],"\""),collapse=", ")," not an available level."))}
  }
  
  
  extra_columns <- lapply(as.list(factorcovariates),FUN=function(cov){
    fvec <- data[,cov]
    temp <- do.call(cbind,lapply(as.list(levels(fvec)),FUN=function(label){return((fvec==label)+0)}))
    colnames(temp) <- paste0(cov,"_",levels(fvec))
    return(temp)
  })
  
  # NEED TO CHOOSE BASE
  if(!is.null(baselevels)){
    basenames <- sapply(1:length(baselevels),FUN=function(x){paste0(names(baselevels)[x],"_",baselevels[x])})
  }else{
    basenames <- unlist(lapply(extra_columns,FUN=function(x){return(colnames(x)[1])}))
  }
  
  # cbind all new factor columns
  extra_columns <- do.call(cbind,extra_columns)
  # delete base levels
  extra_columns <- extra_columns[,!(colnames(extra_columns)%in%basenames)]
  
  return(
    list(
      extra_columns=extra_columns,
      factorbasenames=basenames,
      newcovariates=colnames(extra_columns)
      
    )
  )
}



insert_factorbase <- function(df,factorbasenames){
  
  insertRow <- function(existingDF, newrow, r) {
    existingDF <- rbind(existingDF,newrow)
    existingDF <- existingDF[order(c(1:(nrow(existingDF)-1),r-0.5)),]
    row.names(existingDF) <- 1:nrow(existingDF)
    return(existingDF)  
  }
  
  variablenames <- sapply(factorbasenames,FUN=function(x){
    strsplit(x,split="_")[[1]][1]
  })
  
  
  # Add base levels to parameter df (in the right spot)
  for(i in 1:length(variablenames)){
    r <- which(grepl(pattern=paste0("beta\\_",variablenames[i]),x=rownames(df)))
    temp <- df[r,]
    df <- df[-r,]
    temp <- rbind(c(0,0),temp)
    rownames(temp)[1] <- paste0("beta_",factorbasenames[i]," (base)")
    df <- rbind(df,temp)
  }
  return(df)
}


if(getRversion() >= "2.15.1"){
  globalVariables(c("temp_formula"))
  
}