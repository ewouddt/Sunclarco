# Project: Copula Package
# 
# Author: Gebruiker
###############################################################################


# some logllh
loglik2.pwe_stage1jack <- function(p,num_pieces,status,time.k,Z.k,cutpoints,data.k){
	
	lambdas.k <- exp(p[1:num_pieces])
	betas.k <- p[(num_pieces+1):length(p)]
	
	
	haz.k <- approx(cutpoints[1:num_pieces],lambdas.k,xout=time.k,method='constant',rule=2)$y
	DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
	Inter.k <- c(0,cumsum(lambdas.k*DiffInter))
	cumhaz.k <- approx(cutpoints,Inter.k,xout=time.k,method='linear',rule=2)$y
	

	if(dim(Z.k)[2]==1){
		cov.lincomb <- (betas.k * Z.k[,1])
	}else{
		
		# ALTERNATIVE 1 (use vectorised solution)
#				cov.lincomb <- apply(Z.k,MARGIN=1,FUN=function(row){sum(betas.k * row)})
		
		# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
		data.list <- split(as.matrix(Z.k), seq(nrow(Z.k)))
		cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas.k*row)}))
		names(cov.lincomb) <- NULL
		rm(list="data.list")
	}
	
	s.k <- exp(-cumhaz.k*exp(cov.lincomb)); #u_ij 
	f.k <- haz.k*exp(cov.lincomb)*s.k; #-du_ij/dy_ij
	
	loglik <- data.k[,status]*log(f.k)+(1-data.k[,status])*log(s.k)
	return(-sum(loglik))
}


loglik.2stagejack_CL <- function(p,ClusterData,ClusterDataList,status,k){
	
	theta  <- exp(p)
	
	
	sumG <- 1:length(ClusterDataList)
	sumH <- 1:length(ClusterDataList)
	
	for (i in 1:length(ClusterDataList)){
		ClusterDataList[[i]]$G <- ClusterDataList[[i]][,status]*log(-1/varphi.prime(theta,varphi.inverse(theta,ClusterDataList[[i]]$Sk)))
		ClusterDataList[[i]]$H <- varphi.inverse(theta,ClusterDataList[[i]]$Sk)
		sumG[i] <- sum(ClusterDataList[[i]]$G[ClusterDataList[[i]][,status]==1]);
		sumH[i] <- sum(ClusterDataList[[i]]$H)
	}
	
	loglik <- sapply(ClusterData$ClusterEvents,function(x) ifelse(x==0,0,sum(log(1/theta+seq(0,x-1)))))-
			(1/theta)*log(theta)+ sumG + (-ClusterData$ClusterEvents-1/theta)*log(sumH+1/theta)
	return(-sum(loglik[-k]))
	
}

loglik.2stagejack_GH <- function(p,ClusterData,ClusterDataList,status,k){
	
	theta  <- exp(p)/(1+exp(p))

	sumG <- 1:length(ClusterDataList)
	sumH <- 1:length(ClusterDataList)
	
	for (i in 1:length(ClusterDataList)){
		ClusterDataList[[i]]$G <- ClusterDataList[[i]][,status]*log(-1/varphiGH.prime(theta,varphiGH.inverse(theta,ClusterDataList[[i]]$Sk)))
		ClusterDataList[[i]]$H <- varphiGH.inverse(theta,ClusterDataList[[i]]$Sk)
		sumG[i] <- sum(ClusterDataList[[i]]$G[ClusterDataList[[i]][,status]==1]);
		sumH[i] <- sum(ClusterDataList[[i]]$H)
	}
	
	loglik <- sumG+logdth.deriv_GumbHoug(ClusterData$ClusterEvents,sumH,theta)$logderiv
	
	return(-sum(loglik[-k]))
		
}

#Robust variance via grouped jackknife estimator

jack_2stage_pweCL <- function(data,covariates,status,time,clusters,ClusterData,ClusterDataList,num_pieces,init.values,cutpoints,verbose){
	
	betas2jack_pwe   <-  matrix(NA,nrow=length(ClusterDataList),ncol=length(covariates))
	lambdas2jack_pwe <- matrix(NA,nrow=length(ClusterDataList),ncol=num_pieces)
	theta2jack_pweCL <- 1:length(ClusterDataList)
	
	

	# applying logllh 
	if(verbose){
	  cat("Jackknife started (",length(ClusterDataList)," clusters ): \n")
	  pb <- txtProgressBar(min=1,max=length(ClusterDataList),style=3)
	}
	
	for(k in 1:length(ClusterDataList)){
		
		data.k <- data[data[,clusters]!=k,]
		Z.k <- data.k[,covariates,drop=FALSE]
		time.k <- data.k[,time]
		
		
		res2.pwe_stage1jack <- optim(init.values[-(num_pieces+1)],loglik2.pwe_stage1jack,
				num_pieces=num_pieces,status=status,time.k=time.k,Z.k=Z.k,cutpoints=cutpoints,data.k=data.k,
				method="BFGS")
		
		lambdas2jack_pwe[k,] <- exp(res2.pwe_stage1jack$par[1:num_pieces])
		betas2jack_pwe[k,] <- res2.pwe_stage1jack$par[(num_pieces+1):length(res2.pwe_stage1jack$par)]
		
		
		#inpluggen in Clayton loglikelihood
		haz.kk <- approx(cutpoints[1:num_pieces],lambdas2jack_pwe[k,],xout=time.k,method='constant',rule=2)$y
		DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
		Inter.kk <- c(0,cumsum(lambdas2jack_pwe[k,]*DiffInter))
		cumhaz.kk <- approx(cutpoints,Inter.kk,xout=time.k,method='linear',rule=2)$y
		

		if(dim(Z.k)[2]==1){
			cov.lincomb <- (betas2jack_pwe[k,] * Z.k[,1])
		}else{
			
			# ALTERNATIVE 1 (use vectorised solution)
#				cov.lincomb <- apply(Z.k,MARGIN=1,FUN=function(row){sum(betas2jack_pwe[k,] * row)})
			
			# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
			data.list <- split(as.matrix(Z.k), seq(nrow(Z.k)))
			cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas2jack_pwe[k,]*row)}))
			names(cov.lincomb) <- NULL
			rm(list="data.list")
		}
		
		s.kk <- exp(-cumhaz.kk*exp(cov.lincomb))
		
		for(i in 1:length(ClusterDataList)){
			if(i!=k){
				ClusterDataList[[i]]$Sk <- s.kk[data.k[,clusters]==i]
			}
			else {
				ClusterDataList[[i]]$Sk <- rep(0.123,length(ClusterDataList[[i]][,time]))
			}
		}
		

		res2jack_pweCL <- optim(init.values[num_pieces+1],loglik.2stagejack_CL,
				ClusterData=ClusterData,ClusterDataList=ClusterDataList,status=status,k=k,
				method="BFGS")


		theta2jack_pweCL[k] <- exp(res2jack_pweCL$par) 
		
		if(verbose){
		  setTxtProgressBar(pb,value=k)
		}
	}
	if(verbose){
	  cat("\n Jackknife ended \n")
	  close(pb)
	}
	return(list(lambdas=lambdas2jack_pwe,theta=theta2jack_pweCL,betas=betas2jack_pwe))

}


jack_2stage_coxCL <- function(data,covariates,status,time,clusters,ClusterData,ClusterDataList,init.values,verbose){
	
	betas_jack   <-  matrix(NA,nrow=length(ClusterDataList),ncol=length(covariates))
	theta_jack <- 1:length(ClusterDataList)
	
	if(verbose){
	  cat("Jackknife started (",length(ClusterDataList)," clusters ): \n")
	  pb <- txtProgressBar(min=1,max=length(ClusterDataList),style=3)
	}
		

	for(k in 1:length(ClusterDataList)){
		
		data.k <- data[data[,clusters]!=k,]
		
		Z.k <- data.k[,covariates,drop=FALSE]
		
		
		eval(parse(text=paste0("temp_formula <- Surv(",time,",",status,") ~ ",covariates[1])))
		
		
		if(length(covariates)>1){
			for(i.covariate in covariates[-1]){
				eval(parse(text=paste0("temp_formula <- update(temp_formula, ~ . + ",i.covariate,")")))
			}
		}		
		
		PH_Z.k <- coxph(temp_formula,data=data.k)
		betas.k <- as.numeric(PH_Z.k$coef) 
				
		newdata.zero <- as.data.frame(matrix(0,nrow=1,ncol=length(covariates)))
		colnames(newdata.zero) <- names(PH_Z.k$coeff)
				
		fit0.k <- survfit(PH_Z.k,newdata=newdata.zero)
		S0.k <- fit0.k$surv
		
		#estimated survival probabilities without k'th cluster
		
		for(i in 1:length(ClusterDataList)){
			
			Sk.temp <- c(1:length(ClusterDataList[[i]][,time]))
									
			for(j in 1:length(ClusterDataList[[i]][,time])){
				if(i!=k){
					
					if(length(covariates)==1){
						cov.lincomb <- betas.k*ClusterDataList[[i]][,covariates][j]
					}else{
						
						cov.lincomb <- sum(betas.k*as.numeric(ClusterDataList[[i]][j,covariates]))
						
					}
								
#					beta.k*HEIF[[i]][j] # ORIGINAL
					Sk.temp[j] <- (S0.k[fit0.k$time==ClusterDataList[[i]][,time][j]])^exp(cov.lincomb)
					
				}
				else {
					Sk.temp[j] <-0.123
				}
			}
			ClusterDataList[[i]]$Sk <- Sk.temp
			
		}
		
		res_jack <- optim(init.values,loglik.2stagejack_CL,ClusterData=ClusterData,ClusterDataList=ClusterDataList,status=status,k=k,control=list(maxit=3000))

		betas_jack[k,] <- betas.k
		theta_jack[k] <- exp(res_jack$par)
		
		if(verbose){
		  setTxtProgressBar(pb,value=k)
		}
	}
	if(verbose){
	  cat("\n Jackknife ended \n")
	  close(pb)
	}
	
	return(list(betas=betas_jack,theta=theta_jack))
}


jack_2stage_pweGH <- function(data,covariates,status,time,clusters,ClusterData,ClusterDataList,num_pieces,init.values,cutpoints,verbose){
	
	betas2jack_pwe   <-  matrix(NA,nrow=length(ClusterDataList),ncol=length(covariates))
	lambdas2jack_pwe <- matrix(NA,nrow=length(ClusterDataList),ncol=num_pieces)
	theta2jack_pweGH <- 1:length(ClusterDataList)
	
	
	
	# applying logllh 
	if(verbose){
	  cat("Jackknife started (",length(ClusterDataList)," clusters ): \n")
	  pb <- txtProgressBar(min=1,max=length(ClusterDataList),style=3)
	}
	
	for(k in 1:length(ClusterDataList)){
		
		data.k <- data[data[,clusters]!=k,]
		Z.k <- data.k[,covariates,drop=FALSE]
		time.k <- data.k[,time]
		
		
		res2.pwe_stage1jack <- optim(init.values[-(num_pieces+1)],loglik2.pwe_stage1jack,
				num_pieces=num_pieces,status=status,time.k=time.k,Z.k=Z.k,cutpoints=cutpoints,data.k=data.k,
				method="BFGS")
		
		
		lambdas2jack_pwe[k,] <- exp(res2.pwe_stage1jack$par[1:num_pieces])
		betas2jack_pwe[k,] <- res2.pwe_stage1jack$par[(num_pieces+1):length(res2.pwe_stage1jack$par)]
		
		#inpluggen in Clayton loglikelihood
		haz.kk <- approx(cutpoints[1:num_pieces],lambdas2jack_pwe[k,],xout=time.k,method='constant',rule=2)$y
		DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
		Inter.kk <- c(0,cumsum(lambdas2jack_pwe[k,]*DiffInter))
		cumhaz.kk <- approx(cutpoints,Inter.kk,xout=time.k,method='linear',rule=2)$y
		
		
		if(dim(Z.k)[2]==1){
			cov.lincomb <- (betas2jack_pwe[k,] * Z.k[,1])
		}else{
			
			# ALTERNATIVE 1 (use vectorised solution)
#				cov.lincomb <- apply(Z.k,MARGIN=1,FUN=function(row){sum(betas2jack_pwe[k,] * row)})
			
			# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
			data.list <- split(as.matrix(Z.k), seq(nrow(Z.k)))
			cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas2jack_pwe[k,]*row)}))
			names(cov.lincomb) <- NULL
			rm(list="data.list")
		}
		
		s.kk <- exp(-cumhaz.kk*exp(cov.lincomb))
		
		for(i in 1:length(ClusterDataList)){
			if(i!=k){
				ClusterDataList[[i]]$Sk <- s.kk[data.k[,clusters]==i]
			}
			else {
				ClusterDataList[[i]]$Sk <- rep(0.123,length(ClusterDataList[[i]][,time]))
			}
		}
		
		
		res2jack_pweGH <- optim(init.values[num_pieces+1],loglik.2stagejack_GH,
				ClusterData=ClusterData,ClusterDataList=ClusterDataList,status=status,k=k,
				method="BFGS")

		
		theta2jack_pweGH[k] <- exp(res2jack_pweGH$par)/(1+exp(res2jack_pweGH$par))
				
		if(verbose){
      setTxtProgressBar(pb,value=k)
		}
		
	}
	
	if(verbose){
	  cat("\n Jackknife ended \n")
	  close(pb)
	}
	return(list(lambdas=lambdas2jack_pwe,theta=theta2jack_pweGH,betas=betas2jack_pwe))
	
	
}

jack_2stage_coxGH <- function(data,covariates,status,time,clusters,ClusterData,ClusterDataList,init.values,verbose){
	
	betas_jackGH <-  matrix(NA,nrow=length(ClusterDataList),ncol=length(covariates))
	theta_jackGH <- 1:length(ClusterDataList)
	
	if(verbose){
	  cat("Jackknife started (",length(ClusterDataList)," clusters ): \n")
	  pb <- txtProgressBar(min=1,max=length(ClusterDataList),style=3)
	}
	
	for(k in 1:length(ClusterDataList)){
		
		data.k <- data[data[,clusters]!=k,]
		Z.k <- data.k[,covariates,drop=FALSE]
	
		eval(parse(text=paste0("temp_formula <- Surv(",time,",",status,") ~ ",covariates[1])))
		
		
		if(length(covariates)>1){
			for(i.covariate in covariates[-1]){
				eval(parse(text=paste0("temp_formula <- update(temp_formula, ~ . + ",i.covariate,")")))
				
			}
		}		
		
		PH_Z.k <- coxph(temp_formula,data=data.k)
		betas.k <- as.numeric(PH_Z.k$coef) 
		
		newdata.zero <- as.data.frame(matrix(0,nrow=1,ncol=length(covariates)))
		colnames(newdata.zero) <- names(PH_Z.k$coeff)
		
		fit0.k <- survfit(PH_Z.k,newdata=newdata.zero)
		S0.k <- fit0.k$surv
		
		#estimated survival probabilities without k'th cluster
		
		for(i in 1:length(ClusterDataList)){
			
			Sk.temp <- c(1:length(ClusterDataList[[i]][,time]))
			
			for(j in 1:length(ClusterDataList[[i]][,time])){
				if(i!=k){
					
					if(length(covariates)==1){
						cov.lincomb <- betas.k*ClusterDataList[[i]][,covariates][j]
					}else{
						
						cov.lincomb <- sum(betas.k*as.numeric(ClusterDataList[[i]][j,covariates]))
						
					}
					
#					beta.k*HEIF[[i]][j] # ORIGINAL
					Sk.temp[j] <- (S0.k[fit0.k$time==ClusterDataList[[i]][,time][j]])^exp(cov.lincomb)
					
				}
				else {
					Sk.temp[j] <-0.123
				}
			}
			ClusterDataList[[i]]$Sk <- Sk.temp
			
		}
		
	
		
		res_jackGH <- optim(init.values,loglik.2stagejack_GH,ClusterData=ClusterData,ClusterDataList=ClusterDataList,status=status,k=k,control=list(maxit=3000))

		
		betas_jackGH[k,] <- betas.k
		theta_jackGH[k] <- exp(res_jackGH$par)/(1+exp(res_jackGH$par))
		
		if(verbose){
      setTxtProgressBar(pb,value=k)
		}
	}
	if(verbose){
	  cat("\n Jackknife ended \n")
	  close(pb)
	}
	
	
	return(list(betas=betas_jackGH,theta=theta_jackGH))
		
}


#######################################################################################
#######################################################################################



