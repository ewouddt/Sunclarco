# Project: Copula Package
# 
# Author: lucp8394
###############################################################################


########################################################
## 1-stage: Clayton copula model with Weibull margins ##
########################################################


# p : lambda, rho, theta, beta's 
# NOTE THAT ORDER WAS CHANGED!!!!

loglik.1stage_weibCL <- function(p,data,covariates,status,time,clusters,ClusterData,ClusterDataList){
	
	
	lambda <- exp(p[1])
	rho    <- exp(p[2])
	theta  <- exp(p[3])
		
	betas   <- p[4:length(p)]
	names(betas) <- covariates
	

	if(length(covariates)==1){
		cov.lincomb <- (betas * data[,covariates])
	}else{
		
		# ALTERNATIVE 1 (use vectorised solution)
#				cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas * row)})
		
		# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
		data.list <- split(as.matrix(data[,covariates,drop=FALSE]), seq(nrow(data)))
		cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas*row)}))
		names(cov.lincomb) <- NULL
		rm(list="data.list")
	}
	
	s <- exp(-lambda*exp(cov.lincomb)*data[,time]^rho) #u_ij
	f <- lambda*rho*data[,time]^(rho-1)*exp(cov.lincomb)*s #-du_ij/dy_ij
	
	
	sumG <- 1:length(ClusterDataList)
	sumH <- 1:length(ClusterDataList)
	
	for(i in 1:length(ClusterDataList)){
		ClusterDataList[[i]]$S <- s[data[,clusters]==i] #contains estimated survival probabilities
		ClusterDataList[[i]]$F <- f[data[,clusters]==i] #contains estimated distribution
		
		# C[[i]]* correct here?  (cause for GH they are the same)
		ClusterDataList[[i]]$G <- ClusterDataList[[i]][,status]*log(-ClusterDataList[[i]]$F/varphi.prime(theta,varphi.inverse(theta,ClusterDataList[[i]]$S)))
		
		
#		C[[i]]*log(-F[[i]]/varphi.prime(theta,varphi.inverse(theta,S[[i]])))
		
		ClusterDataList[[i]]$H <- varphi.inverse(theta,ClusterDataList[[i]]$S)
		
		sumG[i] <- sum(ClusterDataList[[i]]$G[ClusterDataList[[i]][,status]==1])
		sumH[i] <- sum(ClusterDataList[[i]]$H )
	}
	
	loglik <- sapply(ClusterData$ClusterEvents,function(x) ifelse(x==0,0,sum(log(1/theta+seq(0,x-1)))))-
			(1/theta)*log(theta)+ sumG + (-ClusterData$ClusterEvents-1/theta)*log(sumH+1/theta)
	
	return(-sum(loglik))
}



loglik.1stage_pweCL <- function(p,cutpoints,num_pieces,data,time,status,covariates,clusters,ClusterData,ClusterDataList,stage2part=FALSE){
	
	if(stage2part){
		lambdas <- exp(p[1:num_pieces])
		betas   <- p[(num_pieces+1):length(p)] 
		names(betas) <- covariates
	}else{
		lambdas <- exp(p[1:num_pieces])
		theta  <- exp(p[num_pieces+1])
		betas   <- p[(num_pieces+2):length(p)] 
		names(betas) <- covariates
	}

	
	
	# Remember that theta and beta are switched  + beta became betas
	
	haz <- approx(cutpoints[1:num_pieces],lambdas,xout=data[,time],method='constant',rule=2)$y
	DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
	Inter <- c(0,cumsum(lambdas*DiffInter))
	cumhaz <- approx(cutpoints,Inter,xout=data[,time],method='linear',rule=2)$y

	if(length(covariates)==1){
		cov.lincomb <- (betas * data[,covariates])
	}else{
		
		# ALTERNATIVE 1 (use vectorised solution)
#				cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas * row)})
		
		# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
		data.list <- split(as.matrix(data[,covariates,drop=FALSE]), seq(nrow(data)))
		cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas*row)}))
		names(cov.lincomb) <- NULL
		rm(list="data.list")
	}
	
	s <- exp(-cumhaz*exp(cov.lincomb)); #u_ij 
	f <- haz*exp(cov.lincomb)*s; #-du_ij/dy_ij
	
	if(stage2part){
		loglik <- data[,status]*log(f)+(1-data[,status])*log(s)
		return(-sum(loglik))
	}
	else{	
 
		sumG <- 1:length(ClusterDataList)
		sumH <- 1:length(ClusterDataList)
	
		for (i in 1:length(ClusterDataList)){
		
		
			ClusterDataList[[i]]$S <- s[data[,clusters]==i]
			ClusterDataList[[i]]$F <- f[data[,clusters]==i]
		
			# correct that "C[[i]] * " is not in front of it anymore for G?
			ClusterDataList[[i]]$G <- log(-ClusterDataList[[i]]$F / varphi.prime(theta,varphi.inverse(theta,ClusterDataList[[i]]$S)))
			ClusterDataList[[i]]$H <- varphi.inverse(theta,ClusterDataList[[i]]$S)
		
			sumG[i] <- sum(ClusterDataList[[i]]$G[ClusterDataList[[i]][,status]==1])
			sumH[i] <- sum(ClusterDataList[[i]]$H )

		}
	
		# Same as previous logllh
		loglik <- sapply(ClusterData$ClusterEvents,function(x) ifelse(x==0,0,sum(log(1/theta+seq(0,x-1)))))-
				(1/theta)*log(theta)+ sumG + (-ClusterData$ClusterEvents-1/theta)*log(sumH+1/theta)
	
		return(-sum(loglik))
	}
}




loglik.1stage_GH <- function(p,cutpoints,num_pieces,data,time,status,covariates,clusters,ClusterData,ClusterDataList,marginal,stage2part=FALSE){
	
	if(marginal=="PiecewiseExp"){
		if(stage2part){
			lambdas <- exp(p[1:num_pieces])
			betas   <- p[(num_pieces+1):length(p)]
			
		}else{
			lambdas <- exp(p[1:num_pieces])
			theta  <- exp(p[num_pieces+1])/(1+exp(p[num_pieces+1]))
			betas   <- p[(num_pieces+2):length(p)]
		}
		
		
		haz <- approx(cutpoints[1:num_pieces],lambdas,xout=data[,time],method='constant',rule=2)$y
		DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
		Inter <- c(0,cumsum(lambdas*DiffInter))
		cumhaz <- approx(cutpoints,Inter,xout=data[,time],method='linear',rule=2)$y
		

		if(length(covariates)==1){
			cov.lincomb <-  (betas * data[,covariates])
		}else{
			
			# ALTERNATIVE 1 (use vectorised solution)
#				cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas * row)})
			
			# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
			data.list <- split(as.matrix(data[,covariates,drop=FALSE]), seq(nrow(data)))
			cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas*row)}))
			names(cov.lincomb) <- NULL
			rm(list="data.list")
		}
		
		s <- exp(-cumhaz*exp(cov.lincomb)) #u_ij 
		f <- haz*exp(cov.lincomb)*s #-du_ij/dy_ij
	}
	else if(marginal=="Weibull"){
		lambda <- exp(p[1])
		rho    <- exp(p[2])
		theta  <- exp(p[3])/(1+exp(p[3]))#between 0 and 1
		betas   <- p[4:length(p)]
		names(betas) <- covariates
		
		# beta and theta switched around + beta is betas
		

		if(length(covariates)==1){
			cov.lincomb <-  (betas * data[,covariates])
		}else{
		
			# ALTERNATIVE 1 (use vectorised solution)
#					cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas * row)})
		
			# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
			data.list <- split(as.matrix(data[,covariates,drop=FALSE]), seq(nrow(data)))
			cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas*row)}))
			names(cov.lincomb) <- NULL
			rm(list="data.list")
		}
	
		s <- exp(-lambda*exp(cov.lincomb)*data[,time]^rho) #u_ij
		f <- lambda*rho*data[,time]^(rho-1)*exp(cov.lincomb)*s #-du_ij/dy_ij
	}
	
	if(stage2part){
		loglik <- data[,status]*log(f)+(1-data[,status])*log(s)
		return(-sum(loglik))
	}else{
		sumG <- 1:length(ClusterDataList)
		sumH <- 1:length(ClusterDataList) 
		
		for (i in 1:length(ClusterDataList)){
			ClusterDataList[[i]]$S <- s[data[,clusters]==i]
			ClusterDataList[[i]]$F <- f[data[,clusters]==i]
			
			ClusterDataList[[i]]$G <- log(-ClusterDataList[[i]]$F/varphiGH.prime(theta,varphiGH.inverse(theta,ClusterDataList[[i]]$S)))
			ClusterDataList[[i]]$H <- varphiGH.inverse(theta,ClusterDataList[[i]]$S)
			
			sumG[i] <- sum(ClusterDataList[[i]]$G[ClusterDataList[[i]][,status]==1])
			sumH[i] <- sum(ClusterDataList[[i]]$H )
			
			
		}
		
		loglik <- sumG+logdth.deriv_GumbHoug(ClusterData$ClusterEvents,sumH,theta)$logderiv;
		return(-sum(loglik))
	}
	

}



loglik.2stage_CL <- function(p,status,ClusterData,ClusterDataList){
	
	theta  <- exp(p)
			
	sumG <- 1:length(ClusterDataList)
	sumH <- 1:length(ClusterDataList)
	
	for (i in 1:length(ClusterDataList)){
		
		ClusterDataList[[i]]$G <- ClusterDataList[[i]][,status] * log(-1/varphi.prime(theta,varphi.inverse(theta,ClusterDataList[[i]]$S)))
		ClusterDataList[[i]]$H <- varphi.inverse(theta,ClusterDataList[[i]]$S)
		 
		sumG[i] <- sum(ClusterDataList[[i]]$G[ClusterDataList[[i]][,status]==1])
		sumH[i] <- sum(ClusterDataList[[i]]$H )
				
	}
	
	loglik <- sapply(ClusterData$ClusterEvents,function(x) ifelse(x==0,0,sum(log(1/theta+seq(0,x-1)))))-
			(1/theta)*log(theta)+ sumG + (-ClusterData$ClusterEvents-1/theta)*log(sumH+1/theta)
	return(-sum(loglik))
		

}




loglik.2stage_GH <- function(p,status,ClusterData,ClusterDataList,marginal){
	
	theta  <- exp(p)/(1+exp(p))

	sumG <- 1:length(ClusterDataList)
	sumH <- 1:length(ClusterDataList)
	
	for (i in 1:length(ClusterDataList)){
		
		if(marginal=="Weibull"){
#			ClusterDataList[[i]]$G <- log(-ClusterDataList[[i]]$F/varphiGH.prime(theta,varphiGH.inverse(theta,ClusterDataList[[i]]$S)))
			ClusterDataList[[i]]$G <- log(-1/varphiGH.prime(theta,varphiGH.inverse(theta,ClusterDataList[[i]]$S)))
			
			
		}else if(marginal=="PiecewiseExp" | marginal=="Cox"){
			ClusterDataList[[i]]$G <- ClusterDataList[[i]][,status] * log(-1/varphiGH.prime(theta,varphiGH.inverse(theta,ClusterDataList[[i]]$S)))
		}
				
		ClusterDataList[[i]]$H <- varphiGH.inverse(theta,ClusterDataList[[i]]$S)
		
		sumG[i] <- sum(ClusterDataList[[i]]$G[ClusterDataList[[i]][,status]==1])
		sumH[i] <- sum(ClusterDataList[[i]]$H )
		
	}
	
	loglik <- sumG+logdth.deriv_GumbHoug(ClusterData$ClusterEvents,sumH,theta)$logderiv
	return(-sum(loglik))
	
}


