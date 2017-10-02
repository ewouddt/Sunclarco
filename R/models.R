# Project: Copula Package
# 
# Author: lucp8394
###############################################################################



## EXAMPLE DATA ##
#' @title Insemination Data
#'
#' @description In dairy cattle, the calving interval (the time between two calvings) should be optimally 
#' between 12 and 13 months.
#' One of the main factors determining the length of the calving interval is the time from parturition
#' to the time of first insemination (Time variable in the data set).
#' The data set includes 181 clusters (Herd variable in the data set) of different sizes.
#' The parity of the cow (0 if multiparous; 1 if primiparous) is added as a covariate (Heifer in the data set). 
#'
#' @format A dataframe with 10513 rows and 5 columns.
#' @name insem
NULL


## IMPORTS ##

#' @importFrom survival survreg Surv cluster survfit coxph
#' @importFrom stats aggregate approx optim quantile update
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL



CopulaModel_1stage <- function(data,time,status,clusters,covariates,init.values=NULL,marginal="Weibull",copula="Clayton",n.piecewise=20,verbose=FALSE,factorbasenames=NULL,optim.method=NULL,lower=-Inf,upper=Inf){
	
	
	###################
	### PREPARATION ###
	###################
	

	### Check parameters  (check if combination possible, length of init.values (lambda's and betas!) )
	if(!(marginal %in% c("Weibull","PiecewiseExp"))){stop(paste0("Parameter 'marginal' can not be ",marginal,". It should be either \"Weibull\" or \"PiecewiseExp\" for 1-stage approach"),call.=FALSE)}
	if(!(copula %in% c("Clayton","GH"))){stop(paste0("Parameter 'copula' can not be ",copula,". It should be either \"Clayton\" or \"GH\" for 1-stage approach"),call.=FALSE)}
	
	## Automatic initial parameters + naming parameters correctly
	if(is.null(init.values)){
		init.values <- estimate_parameters(data=data,time=time,status=status,clusters=clusters,covariates=covariates,n.piecewise=n.piecewise,marginal=marginal,copula=copula,stage=1)
	}else{
	  
	  # Checking correct length of given parameters
	  if(marginal=="Weibull"){
	    correct_length <- 3+length(covariates)
	  }
	  else if(marginal=="PiecewiseExp"){
	    correct_length <- n.piecewise+1+length(covariates)
	  }
	  
	 if(verbose){cat("Initial Parameters:\n")}
	  
	  # All initial values are provided: convert to vector and transform
	  if(length(unlist(init.values))==correct_length){
	    theta_temp <- ifelse(copula=="GH",log(init.values$theta[1]/(1-init.values$theta[1])),log(init.values$theta[1]))
	    
	    if(marginal=="Weibull"){
	      init.values.temp <- c(log(init.values$lambda[1]),log(init.values$rho[1]),theta_temp,init.values$beta)
	      if(verbose){cat("lambda = ",init.values$lambda[1],"\n")}
	      
	   }else if(marginal=="PiecewiseExp"){
	      init.values.temp <- c(log(init.values$lambda),theta_temp,init.values$beta)
	      if(verbose){cat("lambda's = ",paste0(init.values$lambda,collapse="; "),"\n")}
	   }
	   if(verbose){ 
	   cat("theta = ",init.values$theta[1],"\n")
	   cat("beta's = ",paste0(init.values$beta,collapse="; "),"\n")
	   cat("\n")
	   }
	   
	   init.values <- init.values.temp
	    
	  }else{
	    # Not all initial values are provided:
	    # generate all, overwrite with transformed provided values
	    
	    init.values.temp <- estimate_parameters(data=data,time=time,status=status,clusters=clusters,covariates=covariates,n.piecewise=n.piecewise,marginal=marginal,copula=copula,stage=1,verbose=FALSE)
	    
	    if(marginal=="Weibull"){
	      if("lambda"%in%names(init.values)){
	        init.values.temp[1] <- log(init.values$lambda[1])
	      }
	      if("rho"%in%names(init.values)){
	        init.values.temp[2] <- log(init.values$rho[1])
	      }
	      if("theta"%in%names(init.values)){
	        init.values.temp[3] <- ifelse(copula=="GH",log(init.values$theta[1]/(1-init.values$theta[1])),log(init.values$theta[1]))
	      }
	      if("beta"%in%names(init.values)){
	        init.values.temp[(4:(4+length(init.values$beta)-1))] <- init.values$beta
	      }
	      if(verbose){
	      cat("lambda = ",exp(init.values.temp[1]),"\n")
	      cat("rho = ",exp(init.values.temp[2]),"\n")
	      cat("theta = ",ifelse(copula=="GH",exp(init.values.temp[3])/(1+exp(init.values.temp[3])),exp(init.values.temp[3])),"\n")
	      cat("beta's = ",paste0(init.values.temp[4:length(init.values.temp)],collapse="; "),"\n")
	      cat("\n")
	      }
	    
	    }else if(marginal=="PiecewiseExp"){
	      
	      if("lambda"%in%names(init.values)){
	        init.values.temp[1:(length(init.values$lambda))] <- log(init.values$lambda)
	      }
	      if("theta"%in%names(init.values)){
	        init.values.temp[n.piecewise+1] <- ifelse(copula=="GH",log(init.values$theta[1]/(1-init.values$theta[1])),log(init.values$theta[1]))
	      }
	      if("beta"%in%names(init.values)){
	        init.values.temp[((n.piecewise+2):(n.piecewise+2+length(init.values$beta)-1))] <- init.values$beta
	      }
	      if(verbose){
	      cat("lambda's = ",paste0(exp(init.values.temp[1:n.piecewise]),collapse="; "),"\n")
	      cat("theta = ",ifelse(copula=="GH",exp(init.values.temp[n.piecewise+1])/(1+exp(init.values.temp[n.piecewise+1])),exp(init.values.temp[n.piecewise+1])),"\n")
	      cat("beta's = ",paste0(init.values.temp[(n.piecewise+2):length(init.values.temp)],collapse="; "),"\n")
	      cat("\n")
	      }
	        
	    }
	    
	    init.values <- init.values.temp
	    
	  }
	  
	  
	  
    ## OLD INITIAL VALUE CODE (WHEN INPUT WAS STILL VECTOR!)
		# if(length(init.values)!=correct_length){stop(paste0("Parameter 'init.values' has an incorrect length. With given parameters it should be of length ",correct_length),call.=FALSE)}
		# 
		# # Checking if parameter is within bounds + transforming with log
		# if(marginal=="Weibull"){
		# 	if(init.values[1]<=0){stop("Lambda parameter should be strictly larger than 0.",call.=FALSE)}
		# 	init.values[1] <- log(init.values[1]) # lambda
		# 	
		# 	if(init.values[2]<=0){stop("Rho parameter should be strictly larger than 0.",call.=FALSE)}
		# 	init.values[2] <- log(init.values[2]) # rho
		# 	
		# 	if(copula=="GH"){
		# 		if(init.values[3]<=0 | init.values[3]>=1){stop("Theta parameter should be between 0 and 1 for copula=\"GH\".",call.=FALSE)}
		# 		init.values[3] <- log(init.values[3]/(1-init.values[3])) # theta
		# 	}else{
		# 		if(init.values[3]<=0){stop("Theta parameter should be strictly larger than 0.",call.=FALSE)}
		# 		init.values[3] <- log(init.values[3]) # theta
		# 	}
		# 
		# }else if(marginal=="PiecewiseExp"){
		# 	if(sum(init.values[1:n.piecewise]<=0)>0){stop("Lambda parameter should be strictly larger than 0.",call.=FALSE)}
		# 	init.values[1:n.piecewise] <- log(init.values[1:n.piecewise]) # lambdas
		# 	
		# 	if(copula=="GH"){
		# 		if(init.values[n.piecewise+1]<=0 | init.values[n.piecewise+1]>=1){stop("Theta parameter should be between 0 and 1.",call.=FALSE)}
		# 		init.values[n.piecewise+1] <- log(init.values[n.piecewise+1]/(1-init.values[n.piecewise+1])) # theta
		# 	}else{
		# 		if(init.values[n.piecewise+1]<=0){stop("Theta parameter should be strictly larger than 0.",call.=FALSE)}
		# 		init.values[n.piecewise+1] <- log(init.values[n.piecewise+1]) # theta
		# 	}
		# }
	}
		
		
	### CLUSTER VARIABLE NEEDS TO BE CHECKED! IT SHOULD BE A NUMBER FROM 1 TO N. IF NEEDED IT HAS TO BE TRANSFORMED!
	clusters_labels <- sort(unique(data[,clusters]))
	if(!all((clusters_labels-c(1:length(clusters_labels)))==0)){
	  
		newlabels <- c(1:length(clusters_labels))
		
		temp <- sapply(data[,clusters],FUN=function(x){
					index <- which(x == clusters_labels)
					return(newlabels[index])
				})
		data[,clusters] <- temp
	}
	
	### Sort Data: Cluster (needed for rle)
	data <- data[order(data[,clusters]),]
				
	
	### Get Information
	Cluster <- data[,clusters]
	nclusterss <- length(levels(as.factor(Cluster)))
	ClusterData <- data.frame(Cluster=c(1:nclusterss),ClusterSubjects=rle(Cluster)$lengths,ClusterEvents=aggregate(data[,status],by=list(Cluster),FUN=sum)[,2])
	ClusterDataList <- vector("list",length=nclusterss)
	
	for(i.list in 1:nclusterss){
		ClusterDataList[[i.list]] <- data[Cluster==i.list,c(time,status,covariates)]
	}
	
	if(marginal=="PiecewiseExp"){
		
		#define cutpoints for piecewise exponential
		
		# Sort data on time
		data_time <- data[order(data[,time]),]
		#define cutpoints such that all intervals contain same number of events
		cutpoints <- quantile(data_time[,time][data_time[,status]==1], probs = seq(0, 1, length.out=(n.piecewise+1)),names=FALSE)
		cutpoints[1] <- 0
		cutpoints[length(cutpoints)] <- max(cutpoints)+1000  
		num_pieces <- length(cutpoints)-1 #number of intervals between cutpoints
		
		}
	
	#########################
	### CLAYTON & WEIBULL ###
	#########################
	
	if(marginal=="Weibull" & copula=="Clayton"){ 
	  
	  if(is.null(optim.method)){
	    optim.method <- "Nelder-Mead"
	  }
			
		if(verbose){cat("Stage 1 initiated. \n")}
		res1_weibCL <- optim(init.values,
				loglik.1stage_weibCL,
				data=data,time=time,status=status,clusters=clusters,covariates=covariates,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				hessian=TRUE,control=list(maxit=3000),method=optim.method,lower=lower,upper=upper) 
		if(verbose){cat("Stage 1 finalized. \n \n")}
		
		# Get Estimates and put in dataframe
		lambda1_weibCL <- exp(res1_weibCL$par[1]) 
		rho1_weibCL    <- exp(res1_weibCL$par[2]) 
		theta1_weibCL  <- exp(res1_weibCL$par[3]) 
		
		betas_weibCL   <- res1_weibCL$par[4:length(res1_weibCL$par)] 
		
		
		# Preparation for standard errors
		VarCov1_weibCL <- diag(c(lambda1_weibCL,rho1_weibCL,theta1_weibCL,rep(1,length(betas_weibCL))))%*%solve(res1_weibCL$hessian)%*%diag(c(lambda1_weibCL,rho1_weibCL,theta1_weibCL,rep(1,length(betas_weibCL))))
		
		stderrtheta1_weibCL <- sqrt(VarCov1_weibCL[3,3])
		
		#standard errors for lambda, rho,theta, betas
		stderr_all <- sqrt(diag(VarCov1_weibCL))
		
	
		#kendall's tau
		tau1_weibCL <- theta1_weibCL/(2+theta1_weibCL) 
		stderrtau1_weibCL <- (2/(theta1_weibCL+2)^2)*stderrtheta1_weibCL 
		
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambda1_weibCL,rho1_weibCL,theta1_weibCL,betas_weibCL),StandardErrors=c(stderr_all))
		rownames(out1) <- c("lambda","rho","theta",paste0("beta_",covariates))
		
		out2 <- c(tau1_weibCL,stderrtau1_weibCL)
		names(out2) <- c("Estimates","StandardErrors")
		
		out3 <- VarCov1_weibCL
		rownames(out3) <- colnames(out3) <- rownames(out1)
		
		out4 <- res1_weibCL$value
	}
	
	###########################
	### CLAYTON & PIECEWISE ###
	###########################
	
	if(marginal=="PiecewiseExp" & copula=="Clayton"){
		
		if(is.null(optim.method)){
			optim.method <- "BFGS"
		}
		
		if(verbose){cat("Stage 1 initiated. \n")}	
		res1_pweCL <- optim(init.values,
				loglik.1stage_pweCL,cutpoints=cutpoints,num_pieces=num_pieces,data=data,status=status,time=time,covariates=covariates,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				hessian=TRUE,method=optim.method,lower=lower,upper=upper)
		
		if(verbose){cat("Stage 1 finalized. \n \n")}
		
		# Get Estimates and put in dataframe
		lambdas1_pweCL <- exp(res1_pweCL$par[1:num_pieces]) 
		theta1_pweCL <- exp(res1_pweCL$par[num_pieces+1]) 
		betas_pweCL <- res1_pweCL$par[(num_pieces+2):length(res1_pweCL$par)] 
		
		# Preparation for standard errors
		VarCov1_pweCL <- diag(c(lambdas1_pweCL,theta1_pweCL,rep(1,length(betas_pweCL))))%*%solve(res1_pweCL$hessian)%*%diag(c(lambdas1_pweCL,theta1_pweCL,rep(1,length(betas_pweCL))))
		stderrtheta1_pweCL <- sqrt(VarCov1_pweCL[num_pieces+1,num_pieces+1])

		#standard errors for lambdas, theta, betas
		stderr_all <- sqrt(diag(VarCov1_pweCL))
	
		#kendall's tau
		tau1_pweCL <- theta1_pweCL/(2+theta1_pweCL)
		stderrtau1_pweCL <- (2/(theta1_pweCL+2)^2)*stderrtheta1_pweCL
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambdas1_pweCL,theta1_pweCL,betas_pweCL),StandardErrors=c(stderr_all))
		rownames(out1) <- c(paste0("lambda",c(1:length(lambdas1_pweCL))),"theta",paste0("beta_",covariates))
		
		out2 <- c(tau1_pweCL,stderrtau1_pweCL)
		names(out2) <- c("Estimates","StandardErrors")
				
		out3 <- VarCov1_pweCL
		rownames(out3) <- colnames(out3) <- rownames(out1)
		
		out4 <- res1_pweCL$value
	}
	
	#################################
	### GUMBEL-HOUGAARD & WEIBULL ###
	#################################
	
	if(marginal=="Weibull" & copula=="GH"){
		
	  if(is.null(optim.method)){
	    optim.method <- "Nelder-Mead"
	  }
		
		if(verbose){cat("Stage 1 initiated. \n")}
		res1_weibGH <- optim(init.values,
				loglik.1stage_GH,data=data,cutpoints=NULL,num_pieces=NULL,time=time,status=status,clusters=clusters,covariates=covariates,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				marginal=marginal,
				hessian=TRUE,control=list(maxit=3000),method=optim.method,lower=lower,upper=upper) 
		if(verbose){cat("Stage 1 finalized. \n \n")}
		
		# Get Estimates and put in dataframe
		lambda1_weibGH <- exp(res1_weibGH$par[1]) 
		rho1_weibGH    <- exp(res1_weibGH$par[2]) 
		theta1_weibGH  <- exp(res1_weibGH$par[3])/(1+exp(res1_weibGH$par[3]))
		betas1_weibGH   <- res1_weibGH$par[4:length(res1_weibGH$par)] 
		
		# Preparation for standard errors
		VarCov1_weibGH <- diag(c(lambda1_weibGH,rho1_weibGH,exp(res1_weibGH$par[3])/(1+exp(res1_weibGH$par[3]))^2,rep(1,length(betas1_weibGH))))%*%solve(res1_weibGH$hessian)%*%diag(c(lambda1_weibGH,rho1_weibGH,exp(res1_weibGH$par[3])/(1+exp(res1_weibGH$par[3]))^2,rep(1,length(betas1_weibGH))))
		stderrtheta1_weibGH <- sqrt(VarCov1_weibGH[3,3]) 
		
		#standard errors for lambda, rho, theta, betas
		stderr_all <- sqrt(diag(VarCov1_weibGH))
		
		#kendall's tau
		tau1_weibGH <- 1-theta1_weibGH
		stderrtau1_weibGH <- stderrtheta1_weibGH
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambda1_weibGH,rho1_weibGH,theta1_weibGH,betas1_weibGH),StandardErrors=c(stderr_all))
		rownames(out1) <- c("lambda","rho","theta",paste0("beta_",covariates))
	
		out2 <- c(tau1_weibGH,stderrtau1_weibGH)
		names(out2) <- c("Estimates","StandardErrors")
		
		out3 <- VarCov1_weibGH
		rownames(out3) <- colnames(out3) <- rownames(out1)
		
		out4 <- res1_weibGH$value
		
	}
	###################################
	### GUMBEL-HOUGAARD & PIECEWISE ###
	###################################
	
	if(marginal=="PiecewiseExp" & copula=="GH"){
		
		if(is.null(optim.method)){
			optim.method <- "BFGS"
		}

		if(verbose){cat("Stage 1 initiated. \n")}	
		res1_pweGH <- optim(init.values,
				loglik.1stage_GH,data=data,cutpoints=cutpoints,num_pieces=num_pieces,status=status,time=time,covariates=covariates,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				marginal=marginal,
				hessian=TRUE,method=optim.method,lower=lower,upper=upper)
		
		if(verbose){cat("Stage 1 finalized. \n \n")}
		
		# Get Estimates and put in dataframe
		lambdas1_pweGH <- exp(res1_pweGH$par[1:num_pieces])
		theta1_pweGH <- exp(res1_pweGH$par[num_pieces+1])/(1+exp(res1_pweGH$par[num_pieces+1]))  
		betas1_pweGH <- res1_pweGH$par[(num_pieces+2):length(res1_pweGH$par)] 

		# Preparation for standard errors
				
		VarCov1_pweGH <- diag(c(lambdas1_pweGH,exp(res1_pweGH$par[num_pieces+1])/(1+exp(res1_pweGH$par[num_pieces+1]))^2,rep(1,length(betas1_pweGH))))%*%solve(res1_pweGH$hessian)%*%diag(c(lambdas1_pweGH,exp(res1_pweGH$par[num_pieces+1])/(1+exp(res1_pweGH$par[num_pieces+1]))^2,rep(1,length(betas1_pweGH))))
		stderrtheta1_pweGH <- sqrt(VarCov1_pweGH[num_pieces+1,num_pieces+1]) 
		
		#standard errors for lambdas, theta, betas
		stderr_all <- sqrt(diag(VarCov1_pweGH))
				
		#kendall's tau
		tau1_pweGH <- 1-theta1_pweGH
		stderrtau1_pweGH <- stderrtheta1_pweGH
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambdas1_pweGH,theta1_pweGH,betas1_pweGH),StandardErrors=c(stderr_all))
		rownames(out1) <- c(paste0("lambda",c(1:length(lambdas1_pweGH))),"theta",paste0("beta_",covariates))
		
		out2 <- c(tau1_pweGH,stderrtau1_pweGH)
		names(out2) <- c("Estimates","StandardErrors")
		
		out3 <- VarCov1_pweGH
		rownames(out3) <- colnames(out3) <- rownames(out1)
		
		out4 <- res1_pweGH$value
	}
		
		# Give init.values correct naming + save all parameters
		names(init.values) <- rownames(out1)
		parameters <- list(time=time,status=status,clusters=clusters,covariates=covariates,init.values=init.values,
			marginal=marginal,copula=copula,n.piecewise=n.piecewise)
		
		# Add factor base levels (if necessary)
		if(!is.null(factorbasenames)){
		  out1 <- insert_factorbase(df=out1,factorbasenames=factorbasenames)
		}
	
		out <- list(Parameters=out1,Kendall_Tau=out2,ParametersCov=out3,logllh=out4,parameter.call=parameters)
		class(out) <- "Sunclarco"
		return(out)

}




CopulaModel_2stage <- function(data,time,status,clusters,covariates,init.values=NULL,marginal="Weibull",copula="Clayton",n.piecewise=20,verbose=FALSE,factorbasenames=NULL,optim.method=NULL,lower=-Inf,upper=Inf){
	
	###################
	### PREPARATION ###
	###################
	

  ## Check if combination possible
	if(!(marginal %in% c("Weibull","PiecewiseExp","Cox"))){stop(paste0("Parameter 'marginal' can not be ",marginal,". It should be either \"Weibull\" or \"PiecewiseExp\" for 2-stage approach"),call.=FALSE)}
	if(!(copula %in% c("Clayton","GH"))){stop(paste0("Parameter 'copula' can not be ",copula,". It should be either \"Clayton\" or \"GH\" for 2-stage approach"),call.=FALSE)}
	
	
	## Automatic initial parameters + naming parameters correctly
	if(is.null(init.values)){
		
		init.values <- estimate_parameters(data=data,time=time,status=status,clusters=clusters,covariates=covariates,n.piecewise=n.piecewise,marginal=marginal,copula=copula,stage=2)
	}else{
	  
		### Check parameters  ( length of init.values (lambda's and betas!) )
		if(marginal=="Weibull" | marginal=="Cox"){
			correct_length <- 1
		}
		else if(marginal=="PiecewiseExp"){
			correct_length <- n.piecewise+1+length(covariates)
		}
	  
	  
	  if(verbose){cat("Initial Parameters:\n")}
	  
	  
	  # All initial values are provided: convert to vector and transform
	  if(length(unlist(init.values))==correct_length){
	    
	    theta_temp <- ifelse(copula=="GH",log(init.values$theta[1]/(1-init.values$theta[1])),log(init.values$theta[1]))
	    
	    if(marginal=="Weibull" | marginal=="Cox"){
	      init.values.temp <- c(theta_temp)
	      if(verbose){
	      cat("theta = ",init.values$theta[1],"\n")
	      cat("\n")
	      }
	    }else if(marginal=="PiecewiseExp"){
	      init.values.temp <- c(log(init.values$lambda),theta_temp,init.values$beta)
	      if(verbose){
	      cat("lambda's = ",paste0(init.values$lambda,collapse = "; "),"\n")
	      cat("theta = ",init.values$theta[1],"\n")
	      cat("beta's = ",paste0(init.values$beta,collapse="; "),"\n")
	      cat("\n")
	      }
	    }
	    
	    init.values <- init.values.temp
	    
	  }else{
	    
	    # Not all initial values are provided:
	    # generate all, overwrite with transformed provided values
	    
	    init.values.temp <- estimate_parameters(data=data,time=time,status=status,clusters=clusters,covariates=covariates,n.piecewise=n.piecewise,marginal=marginal,copula=copula,stage=2)
	    
	    if(marginal=="Weibull" | marginal=="Cox"){
	      
	      # Since there is only 1 parameter here, this will never happen, no need to code this tbh
	      if("theta"%in%names(init.values)){
	        init.values.temp[1] <- ifelse(copula=="GH",log(init.values$theta[1]/(1-init.values$theta[1])),log(init.values$theta[1]))
	      }
	      if(verbose){
	      cat("theta = ",ifelse(copula=="GH",exp(init.values.temp[1])/(1+exp(init.values.temp[1])),exp(init.values.temp[1])),"\n")
	      cat("\n")
	      }
	    }else if(marginal=="PiecewiseExp"){
	      
	      if("lambda"%in%names(init.values)){
	        init.values.temp[1:(length(init.values$lambda))] <- log(init.values$lambda)
	      }
	      if("theta"%in%names(init.values)){
	        init.values.temp[n.piecewise+1] <- ifelse(copula=="GH",log(init.values$theta[1]/(1-init.values$theta[1])),log(init.values$theta[1]))
	      }
	      if("beta"%in%names(init.values)){
	        init.values.temp[((n.piecewise+2):(n.piecewise+2+length(init.values$beta)-1))] <- init.values$beta
	      }
	      if(verbose){
	      cat("lambda's = ",paste0(exp(init.values.temp[1:n.piecewise]),collapse="; "),"\n")
	      cat("theta = ",ifelse(copula=="GH",exp(init.values.temp[n.piecewise+1])/(1+exp(init.values.temp[n.piecewise+1])),exp(init.values.temp[n.piecewise+1])),"\n")
	      cat("beta's = ",paste0(init.values.temp[(n.piecewise+2):length(init.values.temp)],collapse="; "),"\n")
	      cat("\n")
	      }
	      
	    }
	    
	    init.values <- init.values.temp
	    
	  }
	  
	  
	  ## OLD INITIAL PARAMETERS CODE
		# if(length(init.values)!=correct_length){stop(paste0("Parameter 'init.values' has an incorrect length. With given parameters it should be of length ",correct_length),call.=FALSE)}
		# 
		# # Transform paramteres + check bounds
		# 
		# if(marginal=="Weibull" | marginal=="Cox"){
		# 	if(copula=="GH"){
		# 		if(init.values[1]<=0 | init.values[1]>=1){stop("Theta parameter should be between 0 and 1.",call.=FALSE)}
		# 		init.values[1] <- log(init.values[1]/(1-init.values[1])) # theta
		# 	}else{
		# 		if(init.values[1]<=0){stop("Theta parameter should be strictly larger than 0.",call.=FALSE)}
		# 		init.values[1] <- log(init.values[1]) # theta
		# 	}
		# 
		# }else if(marginal=="PiecewiseExp"){
		# 	if(sum(init.values[1:n.piecewise]<=0)>0){stop("Lambda parameter should be strictly larger than 0.",call.=FALSE)}
		# 	init.values[1:n.piecewise] <- log(init.values[1:n.piecewise]) # lambdas
		# 
		# 	if(copula=="GH"){
		# 		if(init.values[n.piecewise+1]<=0 | init.values[n.piecewise+1]>=1){stop("Theta parameter should be between 0 and 1.",call.=FALSE)}
		# 		init.values[n.piecewise+1] <- log(init.values[n.piecewise+1]/(1-init.values[n.piecewise+1])) # theta
		# 	}else{
		# 		if(init.values[n.piecewise+1]<=0){stop("Theta parameter should be strictly larger than 0.",call.=FALSE)}
		# 		init.values[n.piecewise+1] <- log(init.values[n.piecewise+1]) # theta
		# 	}
		# }
	  
	  
	}
		

	### CLUSTER VARIABLE NEEDS TO BE CHECKED! IT SHOULD BE A NUMBER FROM 1 TO N. IF NEEDED IT HAS TO BE TRANSFORMED!
	clusters_labels <- sort(unique(data[,clusters]))
	
	
	if(!all((clusters_labels-c(1:length(clusters_labels)))==0)){
		
		newlabels <- c(1:length(clusters_labels))
		
		temp <- sapply(data[,clusters],FUN=function(x){
					index <- which(x == clusters_labels)
					return(newlabels[index])
				})
		data[,clusters] <- temp
	}
	
	### Sort Data: Cluster (needed for rle)
	data <- data[order(data[,clusters]),]
	
	
	### Get Information
	Cluster <- data[,clusters]
	
	nclusterss <- length(levels(as.factor(Cluster)))
	
	ClusterData <- data.frame(Cluster=c(1:nclusterss),ClusterSubjects=rle(Cluster)$lengths,ClusterEvents=aggregate(data[,status],by=list(Cluster),FUN=sum)[,2])
	
	ClusterDataList <- vector("list",length=nclusterss)
	
	for(i.list in 1:nclusterss){
		ClusterDataList[[i.list]] <- data[Cluster==i.list,c(time,status,covariates)]
	}
	
	if(marginal=="PiecewiseExp"){
		
		#define cutpoints for piecewise exponential
		
		# Sort data on time
		data_time <- data[order(data[,time]),]
		#define cutpoints such that all intervals contain same number of events
		cutpoints <- quantile(data_time[,time][data_time[,status]==1], probs = seq(0, 1, length.out=(n.piecewise+1)),names=FALSE)
		cutpoints[1] <- 0
		cutpoints[length(cutpoints)] <- max(cutpoints)+1000  
		num_pieces <- length(cutpoints)-1 #number of intervals between cutpoints
		
	}
	
	### Check JackKnife Threshold ###
	if(marginal=="PiecewiseExp"){
		if(length(ClusterDataList)<=(num_pieces + length(covariates))){stop("Not enough clusters for Jackknife.",call.=FALSE)}
	}
	
	#########################
	### CLAYTON & WEIBULL ###
	#########################
	
	if(marginal=="Weibull" & copula=="Clayton"){
		
		if(verbose){cat("Stage 1 initiated. \n")}

		eval(parse(text=paste0("temp_formula <- Surv(",time,",",status,") ~ cluster(",clusters,")")))
		
		for(i.covariate in covariates){
			eval(parse(text=paste0("temp_formula <- update(temp_formula, ~ . + ",i.covariate,")")))
		}
		
		Surv_result   <- survreg(temp_formula,dist="weibull",robust=TRUE,data=data)

		mu <- Surv_result$coeff[1]
		gammas <- Surv_result$coeff[2:length(Surv_result$coeff)]
		sigma <- Surv_result$scale
		
		lambda2_weib <- as.numeric(exp(-mu/sigma)) 
		rho2_weib    <- as.numeric(1/sigma) 
		betas2_weib   <- as.numeric(-gammas/sigma) 
		
		
		#Surv_result$var gives ROBUST variance-covariance matrix for coef[1],coef[2] and log(scale)
		#use delta method to transform log(scale) to scale
		varmu <-Surv_result$var[1,1]
		
		vargammas <- 	matrix(Surv_result$var[c(2:length(Surv_result$coeff)),c(2:length(Surv_result$coeff))],nrow=length(gammas),ncol=length(gammas))
		
		covmugammas <- sapply(c(2:length(Surv_result$coeff)),FUN=function(x){Surv_result$var[1,x]})
		
		covmusigma <- Surv_result$var[1,dim(Surv_result$var)[2]]*Surv_result$scale
	
		covgammassigma <- sapply(c(2:length(Surv_result$coeff)),FUN=function(x){Surv_result$var[x,dim(Surv_result$var)[2]]*Surv_result$scale})
		
		varsigma <- Surv_result$var[dim(Surv_result$var)[2],dim(Surv_result$var)[2]]*Surv_result$scale^2
		
		#Variance-covariance matrix for marginal parameters lambda, rho, beta
		#Klein & Moeschberger (1997) p.378-379

		varbetas <- matrix(NA,nrow=length(gammas),ncol=length(gammas))
		for (j in 1:length(gammas)){
				for (k in 1:length(gammas)){
					varbetas[j,k] <- vargammas[j,k]/sigma^2-gammas[j]*covgammassigma[j]/sigma^3-gammas[k]*covgammassigma[k]/sigma^3+gammas[j]*gammas[k]*varsigma/sigma^4
				}
		}


		varlam=exp(-2*mu/sigma)*(varmu/sigma^2-2*mu*covmusigma/sigma^3+mu^2*varsigma/sigma^4)
		varrho=varsigma/sigma^4
		covbetaslam=exp(-mu/sigma)*(covmugammas/sigma^2-mu*covgammassigma/sigma^3-gammas*covmusigma/sigma^3+gammas*mu*varsigma/sigma^4)
		covbetasrho=covgammassigma/sigma^3-gammas*varsigma/sigma^4
		covlamrho=exp(-mu/sigma)*(covmusigma/sigma^3-mu*varsigma/sigma^4)
		
	
		# Covariance Matrix: Lambda, Rho, Betas

		VarCov2_lambda_rho <- matrix(c(
						varlam,covlamrho,covbetaslam,
						covlamrho,varrho,covbetasrho
									),nrow=2,ncol=(2+length(covbetaslam)),byrow=TRUE)
		
		VarCov2_betas_1 <- matrix(c(covbetaslam,covbetasrho),nrow=length(covbetaslam),ncol=2,byrow=FALSE)
		VarCov2_betas_2 <- varbetas
		
		VarCov2 <- rbind(VarCov2_lambda_rho,cbind(VarCov2_betas_1,VarCov2_betas_2))
		rownames(VarCov2) <- colnames(VarCov2) <- c("lambda","rho",paste0("beta_",covariates))
		
		IVI <- VarCov2
		
		stderrbetas2_weib <- sapply(c(3:dim(VarCov2)[2]),FUN=function(x){sqrt(VarCov2[x,x])})
		
	
		if(length(covariates)==1){
			cov.lincomb <- (betas2_weib * data[,covariates])
		}else{
			
			# ALTERNATIVE 1 (use vectorised solution)
#			cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas2_weib * row)})
		
			# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
			data.list <- split(as.matrix(data[,covariates,drop=FALSE]), seq(nrow(data)))
			cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas2_weib*row)}))
			names(cov.lincomb) <- NULL
			rm(list="data.list")
		}



		s <- exp(-lambda2_weib*exp(cov.lincomb)*data[,time]^rho2_weib) #u_ij
		f <- lambda2_weib*rho2_weib*data[,time]^(rho2_weib-1)*exp(cov.lincomb)*s #-du_ij/dy_ij
		
		
		for (i in 1:length(ClusterDataList)){
			ClusterDataList[[i]]$S <- s[data[,clusters]==i]
			ClusterDataList[[i]]$F <- f[data[,clusters]==i]
		}
		if(verbose){cat("Stage 1 finalized. \n \n")}
		if(verbose){cat("Stage 2 initiated. \n")}
		
		if(is.null(optim.method)){
			optim.method <- "Brent"
			lower <- min(c(0,init.values-1))
			upper <- init.values+1
		}
		
		res2_weibCL <- optim(init.values,
				loglik.2stage_CL,
				ClusterData=ClusterData,ClusterDataList=ClusterDataList,status=status,
				hessian=TRUE,control=list(maxit=3000),method=optim.method,lower=lower,upper=upper)
		

		theta2_weibCL <- exp(res2_weibCL$par)

		
		if(verbose){
			cat("Estimates: \n")
			cat(c(lambda2_weib,rho2_weib,betas2_weib,theta2_weibCL),fill=1,labels=c("Lambda =","Rho =",paste0("beta_",covariates," ="),"Theta ="))
				}
		
		
		
		#get estimates for I12 en I22 (not using complete one-stage results)
		#estimate for I22
		vartheta2stage_weibCL <- theta2_weibCL*solve(res2_weibCL$hessian)*theta2_weibCL #estimate for 1/I22
		stderrtheta2stage_weibCL <- sqrt(vartheta2stage_weibCL) #estimate for sqrt(1/I22)
		
		#estimate for I12
		#use hessian from one-stage loglikelihood, evaluated at 2-stage solution
		

		res12_weibCL <- optim(c(log(lambda2_weib),log(rho2_weib),log(theta2_weibCL),betas2_weib),
						loglik.1stage_weibCL,data=data,covariates=covariates,status=status,time=time,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
						hessian=TRUE,control=list(maxit=1))
				
		if(verbose){cat("Stage 2 finalized. \n \n")}
		
		VarCov12_weibCL <- diag(c(lambda2_weib,rho2_weib,theta2_weibCL,rep(1,length(betas2_weib))))%*%solve(res12_weibCL$hessian)%*%diag(c(lambda2_weib,rho2_weib,theta2_weibCL,rep(1,length(betas2_weib))))
		
		
		II_weibCL.inv <- VarCov12_weibCL
		II_weibCL <- solve(VarCov12_weibCL)
		II_weibCL11 <- II_weibCL[-3,-3]
		II_weibCL21 <- II_weibCL[3,-3]
		II_weibCL12 <- II_weibCL[-3,3]
		II_weibCL22 <- II_weibCL[3,3]
		
		vartheta2_weibCL <- vartheta2stage_weibCL+II_weibCL21%*%IVI%*%II_weibCL12*vartheta2stage_weibCL^2 
		stderrtheta2_weibCL <- sqrt(vartheta2_weibCL)  
		
				
		#kendall's tau
		tau2_weibCL <- theta2_weibCL/(2+theta2_weibCL) 
		stderrtau2_weibCL <- (2/(theta2_weibCL+2)^2)*stderrtheta2_weibCL # OR stderrtheta2_weibCLa FROM COMPARING?!
		
		# VarCov2 holds the SE of the others
		stderr_all <- c(sqrt(diag(VarCov2)),stderrtheta2_weibCL)
		

		# Correct estimates and standard errors? No transformation required?
		# Which estimates do you take? From the first or from res12_weibCL ?
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambda2_weib,rho2_weib,betas2_weib,theta2_weibCL),StandardErrors=c(stderr_all))
		rownames(out1) <- c(rownames(VarCov2),"theta")

		out2 <- c(tau2_weibCL,stderrtau2_weibCL)
		names(out2) <- c("Estimates","StandardErrors")

		# DOES THETA VARIANCE HAVE TO BE REPLACED?
		out3 <- VarCov2
#		rownames(out3) <- colnames(out3) <- rownames(out1)
		out4 <- res2_weibCL$value
				
#		Give init.values correct naming 
		names(init.values) <- "theta"
	}
	
	###########################
	### CLAYTON & PIECEWISE ###
	###########################
	
	if(marginal=="PiecewiseExp" & copula=="Clayton"){
		

		if(verbose){cat("Stage 1 initiated. \n")}
		
		
		# input lambda's, beta's	
		res2.pwe_stage1 <- optim(init.values[-(num_pieces+1)],loglik.1stage_pweCL,
				cutpoints=cutpoints,num_pieces=num_pieces,data=data,time=time,status=status,covariates=covariates,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,stage2part=TRUE,
				hessian=TRUE,method="BFGS")
	
		# (estimates)
		lambdas2_pwe <- exp(res2.pwe_stage1$par[1:num_pieces])
		betas2_pwe <- res2.pwe_stage1$par[(num_pieces+1):(length(res2.pwe_stage1$par))]

	
		#plug in into Clayton loglikelihood
		haz <- approx(cutpoints[1:num_pieces],lambdas2_pwe,xout=data[,time],method='constant',rule=2)$y
		DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
		Inter <- c(0,cumsum(lambdas2_pwe*DiffInter))
		cumhaz <- approx(cutpoints,Inter,xout=data[,time],method='linear',rule=2)$y
	
		
		if(length(covariates)==1){
			cov.lincomb <- (betas2_pwe * data[,covariates])
		}else{
			
			# ALTERNATIVE 1 (use vectorised solution)
#			cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas2_pwe * row)})
			
			# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
			data.list <- split(as.matrix(data[,covariates,drop=FALSE]), seq(nrow(data)))
			cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas2_pwe*row)}))
			names(cov.lincomb) <- NULL
			rm(list="data.list")
		}
		
		
		s <- exp(-cumhaz*exp(cov.lincomb)) #u_ij 
		f <- haz*exp(cov.lincomb)*s #-du_ij/dy_ij
	

		for(i in 1:length(ClusterDataList)){
			ClusterDataList[[i]]$S <- s[data[,clusters]==i] #contains estimated survival probabilities
			ClusterDataList[[i]]$F <- f[data[,clusters]==i] #contains estimated distribution
		}
		
	
		#plug marginal estimates into loglikelihood
		if(verbose){cat("Stage 1 finalized. \n \n")}
		if(verbose){cat("Stage 2 initialized. \n")}
		
		if(is.null(optim.method)){
			optim.method <- "Brent"
			lower <- min(c(0,init.values[num_pieces+1]-1))
			upper <- init.values[num_pieces+1]+1
		}
		
		# theta input
		res2_pweCL <- optim(init.values[num_pieces+1],
				loglik.2stage_CL,status=status,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				hessian=TRUE,method=optim.method,lower=lower,upper=upper)

		theta2_pweCL <- exp(res2_pweCL$par) 
		
		if(verbose){
			cat("Estimates: \n")
			cat(c(lambdas2_pwe,betas2_pwe,theta2_pweCL),fill=1,labels=c(paste0("Lambda",c(1:n.piecewise)," ="),paste0("beta_",covariates," ="),"Theta ="))
		}
				
				
		
		init.values.jack <-  c(log(lambdas2_pwe),log(theta2_pweCL),betas2_pwe)# lambda's,theta, beta's
		jack_out <- jack_2stage_pweCL(data=data,covariates=covariates,status=status,time=time,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,num_pieces=num_pieces,init.values=init.values.jack,cutpoints=cutpoints,verbose=verbose,optim.method=optim.method,lower=lower,upper=upper)
	
		
		if(verbose){cat("Stage 2 finalized. \n \n")}
		betas2jack_pwe <- jack_out$betas
		theta2jack_pweCL <- jack_out$theta
		lambdas2jack_pwe <- jack_out$lambdas
		
		stderrbetas2_pwe <- c()
		for(i in 1:dim(betas2jack_pwe)[2]){
			stderrbetas2_pwe <- c(stderrbetas2_pwe,sqrt((length(ClusterDataList)-(num_pieces+length(covariates)))/length(ClusterDataList)*sum((betas2jack_pwe[,i]-betas2_pwe[i])^2)))
		}
		## More than 1 beta, still num_pieces+1
		
		stderrtheta2_pweCL <- sqrt((length(ClusterDataList)-(num_pieces+length(covariates)+1))/length(ClusterDataList)*sum((theta2jack_pweCL-theta2_pweCL)^2)) 
		# still num_pieces+2 (see above)


		# How about SE for lambdas?
		stderrlambdas2_pwe <- c()
		for(i in 1:dim(lambdas2jack_pwe)[2]){
			stderrlambdas2_pwe <- c(stderrlambdas2_pwe,sqrt((length(ClusterDataList)-(num_pieces+length(covariates)))/length(ClusterDataList)*sum((lambdas2jack_pwe[,i]-lambdas2_pwe[i])^2)))
		}
				
	
		stderr_all <- c(stderrlambdas2_pwe,stderrtheta2_pweCL,stderrbetas2_pwe)
		
		#kendall's tau
		tau2_pweCL <- theta2_pweCL/(2+theta2_pweCL) 
		stderrtau2_pweCL <- (2/(theta2_pweCL+2)^2)*stderrtheta2_pweCL 
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambdas2_pwe,theta2_pweCL,betas2_pwe),StandardErrors=c(stderr_all))
		rownames(out1) <- c(paste0("lambda",c(1:n.piecewise)),"theta",paste0("beta_",covariates))
		
		out2 <- c(tau2_pweCL,stderrtau2_pweCL)
		names(out2) <- c("Estimates","StandardErrors")
		
		# SHOULD COVARIANCE BE HERE?
		out3 <- NA
#		rownames(out3) <- colnames(out3) <- rownames(out1)
		

		# CORRECT LIKELIHOOD???
		out4 <- res2_pweCL$value
		
#		Give init.values correct naming 
		names(init.values) <- rownames(out1)

	}
	
	#####################
	### CLAYTON & COX ###
	#####################
	
	if(marginal=="Cox" & copula=="Clayton"){
		if(verbose){cat("Stage 1 initiated. \n")}
		
		eval(parse(text=paste0("temp_formula <- Surv(",time,",",status,") ~ cluster(",clusters,")")))
		
		for(i.covariate in covariates){
			eval(parse(text=paste0("temp_formula <- update(temp_formula, ~ . + ",i.covariate,")")))
  	}
		
		# note: as.factor(heifer)
		
		PH_COX <- coxph(temp_formula,data=data)
		betas2_semipar <- as.numeric(PH_COX$coeff) 
		stderrbetas2_semipar <- if(length(PH_COX$coeff)>1) sqrt(diag(PH_COX$var)) else sqrt(PH_COX$var)
		  
		
		
		newdata.zero <- as.data.frame(matrix(0,nrow=1,ncol=length(covariates)))
		colnames(newdata.zero) <- names(PH_COX$coeff)
		
		
		fit0 <- survfit(PH_COX,newdata=newdata.zero) # put all covariates at 0
		S0 <- fit0$surv
		
		for(i in 1:length(ClusterDataList)){
			
			if(length(covariates)==1){
				cov.lincomb <- (betas2_semipar * ClusterDataList[[i]][,covariates])
			}else{
			
				# ALTERNATIVE 1 (use vectorised solution)
#				cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas2_semipar * row)})
			
				# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
				data.list <- split(as.matrix(ClusterDataList[[i]][,covariates,drop=FALSE]), seq(nrow(ClusterDataList[[i]])))
				cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas2_semipar*row)}))
				names(cov.lincomb) <- NULL
				rm(list="data.list")
			}
		
		
			temp.S <- c(1:length(ClusterDataList[[i]][,time]))
			for(j in 1:length(ClusterDataList[[i]][,time])){
				temp.S[j] <- (S0[fit0$time==ClusterDataList[[i]][,time][j]])^exp(cov.lincomb[j])
			}
			
			ClusterDataList[[i]]$S <- temp.S
		}
		
		
		if(verbose){cat("Stage 1 finalized. \n \n")}
		if(verbose){cat("Stage 2 initiated. \n")}
		
		if(is.null(optim.method)){
			optim.method <- "Brent"
			lower <- min(c(0,init.values-1))
			upper <- init.values+1
		}
		
		#plug marginal estimates into loglikelihood 
		
		
		# init.values = theta
		res2_semiparCL <- optim(init.values,loglik.2stage_CL,status=status,ClusterData=ClusterData,ClusterDataList=ClusterDataList,control=list(maxit=3000),method=optim.method,lower=lower,upper=upper)
		theta2_semiparCL <- exp(res2_semiparCL$par)
		
		
		if(verbose){
			cat("Estimates: \n")
			cat(c(betas2_semipar,theta2_semiparCL),fill=1,labels=c(paste0("beta_",covariates," ="),"Theta ="))
		}
		
		#jackknife procedure to get standard error

		init.values.jack <-  c(log(theta2_semiparCL))# theta
		jack_out <- jack_2stage_coxCL(data=data,covariates=covariates,status=status,time=time,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,init.values=init.values.jack,verbose=verbose,optim.method=optim.method,lower=lower,upper=upper)

		stderrtheta2_semiparCL <- sqrt((length(ClusterDataList)-2)/length(ClusterDataList)*sum((jack_out$theta-theta2_semiparCL)^2))

		if(verbose){cat("Stage 2 finalized. \n \n")}
		
				
		#kendall's tau
		tau2_semiparCL <- theta2_semiparCL/(2+theta2_semiparCL)
		stderrtau2_semiparCL <- (2/(theta2_semiparCL+2)^2)*stderrtheta2_semiparCL 
		
		# Standard Errors
	
		stderr_all <- c(stderrtheta2_semiparCL,stderrbetas2_semipar)

		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(theta2_semiparCL,betas2_semipar),StandardErrors=c(stderr_all))
		rownames(out1) <- c("theta",paste0("beta_",covariates))
		
		out2 <- c(tau2_semiparCL,stderrtau2_semiparCL)
		names(out2) <- c("Estimates","StandardErrors")
#		
		# SHOULD COVARIANCE BE HERE?
		out3 <- NA
		##		rownames(out3) <- colnames(out3) <- rownames(out1)
		
		
		# CORRECT LIKELIHOOD???
		out4 <- res2_semiparCL$value
		
#		Give init.values correct naming 
		names(init.values) <- "theta"
		
	}
	
	
	#################################
	### GUMBEL-HOUGAARD & WEIBULL ###
	#################################
	
	if(marginal=="Weibull" & copula=="GH"){
		
		if(verbose){cat("Stage 1 initiated. \n")}
		
		eval(parse(text=paste0("temp_formula <- Surv(",time,",",status,") ~ cluster(",clusters,")")))
		
		for(i.covariate in covariates){
			eval(parse(text=paste0("temp_formula <- update(temp_formula, ~ . + ",i.covariate,")")))
		}
		
		Surv_result   <- survreg(temp_formula,dist="weibull",robust=TRUE,data=data)
		
		mu <- Surv_result$coeff[1]
		gammas <- Surv_result$coeff[2:length(Surv_result$coeff)]
		sigma <- Surv_result$scale
		
		lambda2_weib <- as.numeric(exp(-mu/sigma)) 
		rho2_weib    <- as.numeric(1/sigma) 
		betas2_weib   <- as.numeric(-gammas/sigma) 
		
		#Surv_result$var gives ROBUST variance-covariance matrix for coef[1],coef[2] and log(scale)
		#use delta method to transform log(scale) to scale
		varmu <-Surv_result$var[1,1]
		
		vargammas <- 	matrix(Surv_result$var[c(2:length(Surv_result$coeff)),c(2:length(Surv_result$coeff))],nrow=length(gammas),ncol=length(gammas))
		
		covmugammas <- sapply(c(2:length(Surv_result$coeff)),FUN=function(x){Surv_result$var[1,x]})
		
		covmusigma <- Surv_result$var[1,dim(Surv_result$var)[2]]*Surv_result$scale
		
		covgammassigma <- sapply(c(2:length(Surv_result$coeff)),FUN=function(x){Surv_result$var[x,dim(Surv_result$var)[2]]*Surv_result$scale})
		
		varsigma <- Surv_result$var[dim(Surv_result$var)[2],dim(Surv_result$var)[2]]*Surv_result$scale^2
		
		
		#Variance-covariance matrix for marginal parameters lambda, rho, beta
		#Klein & Moeschberger (1997) p.378-379
		
		# IS THE FOLLOWING STILL OKAY FOR GUMBEL????
		varbetas <- matrix(NA,nrow=length(gammas),ncol=length(gammas))
		for (j in 1:length(gammas)){
			for (k in 1:length(gammas)){
				varbetas[j,k] <- vargammas[j,k]/sigma^2-gammas[j]*covgammassigma[j]/sigma^3-gammas[k]*covgammassigma[k]/sigma^3+gammas[j]*gammas[k]*varsigma/sigma^4
			}
		}
		
		varlam=exp(-2*mu/sigma)*(varmu/sigma^2-2*mu*covmusigma/sigma^3+mu^2*varsigma/sigma^4)
		varrho=varsigma/sigma^4
		covbetaslam=exp(-mu/sigma)*(covmugammas/sigma^2-mu*covgammassigma/sigma^3-gammas*covmusigma/sigma^3+gammas*mu*varsigma/sigma^4)
		covbetasrho=covgammassigma/sigma^3-gammas*varsigma/sigma^4
		covlamrho=exp(-mu/sigma)*(covmusigma/sigma^3-mu*varsigma/sigma^4)
		
		
		VarCov2_lambda_rho <- matrix(c(
						varlam,covlamrho,covbetaslam,
						covlamrho,varrho,covbetasrho
				),nrow=2,ncol=(2+length(covbetaslam)),byrow=TRUE)
		
		VarCov2_betas_1 <- matrix(c(covbetaslam,covbetasrho),nrow=length(covbetaslam),ncol=2,byrow=FALSE)
		VarCov2_betas_2 <- varbetas
		
		VarCov2 <- rbind(VarCov2_lambda_rho,cbind(VarCov2_betas_1,VarCov2_betas_2))
		rownames(VarCov2) <- colnames(VarCov2) <- c("lambda","rho",paste0("beta_",covariates))
		
		IVI <- VarCov2
		
		stderrbetas2_weib <- sapply(c(3:dim(VarCov2)[2]),FUN=function(x){sqrt(VarCov2[x,x])})
		
		
		if(length(covariates)==1){
			cov.lincomb <- (betas2_weib * data[,covariates])
		}else{
			
			# ALTERNATIVE 1 (use vectorised solution)
#			cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas2_weib * row)})
			
			# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
			data.list <- split(as.matrix(data[,covariates,drop=FALSE]), seq(nrow(data)))
			cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas2_weib*row)}))
			names(cov.lincomb) <- NULL
			rm(list="data.list")
		}
		
		
		
		s <- exp(-lambda2_weib*exp(cov.lincomb)*data[,time]^rho2_weib) #u_ij
		f <- lambda2_weib*rho2_weib*data[,time]^(rho2_weib-1)*exp(cov.lincomb)*s #-du_ij/dy_ij
		
		for (i in 1:length(ClusterDataList)){
			ClusterDataList[[i]]$S <- s[data[,clusters]==i]
			ClusterDataList[[i]]$F <- f[data[,clusters]==i]
		}
		if(verbose){cat("Stage 1 finalized. \n \n")}
		if(verbose){cat("Stage 2 initiated. \n")}
		
		if(is.null(optim.method)){
			optim.method <- "BFGS"
		}
		
		res2_weibGH <- optim(init.values,
				loglik.2stage_GH,marginal=marginal,
				ClusterData=ClusterData,ClusterDataList=ClusterDataList,status=status,
				hessian=TRUE,method=optim.method,control=list(maxit=3000),lower=lower,upper=upper)
		
		theta2_weibGH <- exp(res2_weibGH$par)/(1+exp(res2_weibGH$par)) 
		
		if(verbose){
			cat("Estimates: \n")
			cat(c(lambda2_weib,rho2_weib,betas2_weib,theta2_weibGH),fill=1,labels=c("Lambda =","Rho =",paste0("beta_",covariates," ="),"Theta ="))
		}
		
		#get estimates for I12 en I22 (not using complete one-stage results)
		#estimate for I22
		deltameth_der <- exp(res2_weibGH$par)/(1+exp(res2_weibGH$par))^2
		vartheta2stage_weibGH <- deltameth_der^2*solve(res2_weibGH$hessian) #estimate for 1/I22
		stderrtheta2stage_weibGH <- sqrt(vartheta2stage_weibGH) #estimate for sqrt(1/I22)
		
		
		
		#estimate for I12
		#use hessian from one-stage loglikelihood, evaluated at 2-stage solution
		res12_weibGHa <- optim(c(log(lambda2_weib),log(rho2_weib),log(theta2_weibGH/(1-theta2_weibGH)),betas2_weib),loglik.1stage_GH,marginal=marginal,cutpoints=NULL,num_pieces=NULL,data=data,time=time,status=status,covariates=covariates,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,hessian=TRUE,method="BFGS",control=list(maxit=1))
		
		if(verbose){cat("Stage 2 finalized. \n \n")}
		
		VarCov12_weibGHa <- diag(c(lambda2_weib,rho2_weib,deltameth_der,rep(1,length(betas2_weib))))%*%solve(res12_weibGHa$hessian)%*%diag(c(lambda2_weib,rho2_weib,deltameth_der,rep(1,length(betas2_weib))))
		
		II_weibGH.inv <- VarCov12_weibGHa
		II_weibGH <- solve(II_weibGH.inv)
		II_weibGH11 <- II_weibGH[-3,-3]
		II_weibGH21 <- II_weibGH[3,-3]
		II_weibGH12 <- II_weibGH[-3,3]
		II_weibGH22 <- II_weibGH[3,3]
		
		vartheta2_weibGH <- vartheta2stage_weibGH+II_weibGH21%*%IVI%*%II_weibGH12*vartheta2stage_weibGH^2 
		stderrtheta2_weibGH <- sqrt(vartheta2_weibGH) 
		
		#kendall's tau
		tau2_weibGH <- 1-theta2_weibGH
		stderrtau2_weibGH <- stderrtheta2_weibGH  # OR stderrtheta2_weibGHa FROM COMPARING?!
		
		# VarCov2 holds the SE of the others
		stderr_all <- c(sqrt(diag(VarCov2)),stderrtheta2_weibGH)
		
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambda2_weib,rho2_weib,betas2_weib,theta2_weibGH),StandardErrors=c(stderr_all))
		rownames(out1) <- c(rownames(VarCov2),"theta")
		
		out2 <- c(tau2_weibGH,stderrtau2_weibGH)
		names(out2) <- c("Estimates","StandardErrors")
		
		out3 <- VarCov2
#		rownames(out3) <- colnames(out3) <- rownames(out1)
		out4 <- res2_weibGH$value
		
#		Give init.values correct naming 
		names(init.values) <- "theta"
		
	}
	
	###################################
	### GUMBEL-HOUGAARD & PIECEWISE ###
	###################################
	
	if(marginal=="PiecewiseExp" & copula=="GH"){
		
		if(verbose){cat("Stage 1 initiated. \n")}
		
		
		# input lambda's, beta's	
		res2.pwe_stage1 <- optim(init.values[-(num_pieces+1)],loglik.1stage_GH,marginal=marginal,stage2part=TRUE,
				cutpoints=cutpoints,num_pieces=num_pieces,data=data,time=time,status=status,covariates=covariates,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				hessian=TRUE,method="BFGS")
		
		# (estimates)
		lambdas2_pwe <- exp(res2.pwe_stage1$par[1:num_pieces])
		betas2_pwe <- res2.pwe_stage1$par[(num_pieces+1):(length(res2.pwe_stage1$par))]
		
		#plug in into Clayton loglikelihood
		haz <- approx(cutpoints[1:num_pieces],lambdas2_pwe,xout=data[,time],method='constant',rule=2)$y
		DiffInter <- cutpoints[-1]-cutpoints[1:num_pieces]
		Inter <- c(0,cumsum(lambdas2_pwe*DiffInter))
		cumhaz <- approx(cutpoints,Inter,xout=data[,time],method='linear',rule=2)$y
		
		if(length(covariates)==1){
			cov.lincomb <- (betas2_pwe * data[,covariates])
		}else{
			
			# ALTERNATIVE 1 (use vectorised solution)
#			cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas2_pwe * row)})
			
			# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
			data.list <- split(as.matrix(data[,covariates,drop=FALSE]), seq(nrow(data)))
			cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas2_pwe*row)}))
			names(cov.lincomb) <- NULL
			rm(list="data.list")
		}
		
		s <- exp(-cumhaz*exp(cov.lincomb)) #u_ij 
		f <- haz*exp(cov.lincomb)*s #-du_ij/dy_ij
		
		for(i in 1:length(ClusterDataList)){
			ClusterDataList[[i]]$S <- s[data[,clusters]==i] #contains estimated survival probabilities
			ClusterDataList[[i]]$F <- f[data[,clusters]==i] #contains estimated distribution
		}
		
		#plug marginal estimates into loglikelihood
		if(verbose){cat("Stage 1 finalized. \n \n")}
		if(verbose){cat("Stage 2 initialized. \n")}
		
		if(is.null(optim.method)){
			optim.method <- "Brent"
			lower <- min(c(init.values[num_pieces+1]-1,0))
			upper <- init.values[num_pieces+1]+1
		}
		
		
		# theta input

		#res2_pweGH <- optim(init.values[num_pieces+1],loglik.2stage_GH,
		#		marginal=marginal,status=status,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
		#		hessian=TRUE,method="BFGS")
		res2_pweGH <- optim(init.values[num_pieces+1],loglik.2stage_GH,
				marginal=marginal,status=status,ClusterData=ClusterData,ClusterDataList=ClusterDataList,
				hessian=TRUE,control=list(maxit=3000),method=optim.method,lower=lower,upper=upper)

		theta2_pweGH <- exp(res2_pweGH$par)/(1+exp(res2_pweGH$par)) 

		
		if(verbose){
			cat("Estimates: \n")
			cat(c(lambdas2_pwe,betas2_pwe,theta2_pweGH),fill=1,labels=c(paste0("Lambda",c(1:n.piecewise)," ="),paste0("beta_",covariates," ="),"Theta ="))
		}
		
		init.values.jack <-  c(log(lambdas2_pwe),log(theta2_pweGH/(1-theta2_pweGH)),betas2_pwe)# lambda's,theta, beta's
		jack_out <- jack_2stage_pweGH(data=data,covariates=covariates,status=status,time=time,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,num_pieces=num_pieces,init.values=init.values.jack,cutpoints=cutpoints,verbose=verbose,optim.method=optim.method,lower=lower,upper=upper)
		
		if(verbose){cat("Stage 2 finalized. \n \n")}
		betas2jack_pwe <- jack_out$betas
		theta2jack_pweGH <- jack_out$theta
		lambdas2jack_pwe <- jack_out$lambdas
		
		stderrbetas2_pwe <- c()
		for(i in 1:dim(betas2jack_pwe)[2]){
			stderrbetas2_pwe <- c(stderrbetas2_pwe,sqrt((length(ClusterDataList)-(num_pieces+length(covariates)))/length(ClusterDataList)*sum((betas2jack_pwe[,i]-betas2_pwe[i])^2)))
		}
		stderrtheta2_pweGH <- sqrt((length(ClusterDataList)-(num_pieces+length(covariates)+1))/length(ClusterDataList)*sum((theta2jack_pweGH-theta2_pweGH)^2)) 
		
		stderrlambdas2_pwe <- c()
		for(i in 1:dim(lambdas2jack_pwe)[2]){
			stderrlambdas2_pwe <- c(stderrlambdas2_pwe,sqrt((length(ClusterDataList)-(num_pieces+length(covariates)))/length(ClusterDataList)*sum((lambdas2jack_pwe[,i]-lambdas2_pwe[i])^2)))
		}
		
		stderr_all <- c(stderrlambdas2_pwe,stderrtheta2_pweGH,stderrbetas2_pwe)
		
		#kendall's tau
		tau2_pweGH <- 1-theta2_pweGH 
		stderrtau2_pweGH <- stderrtheta2_pweGH 
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(lambdas2_pwe,theta2_pweGH,betas2_pwe),StandardErrors=c(stderr_all))
		rownames(out1) <- c(paste0("lambda",c(1:n.piecewise)),"theta",paste0("beta_",covariates))
		
		out2 <- c(tau2_pweGH,stderrtau2_pweGH)
		names(out2) <- c("Estimates","StandardErrors")
		
		# SHOULD COVARIANCE BE HERE?
		out3 <- NA
		##		rownames(out3) <- colnames(out3) <- rownames(out1)
#		
		
		# CORRECT LIKELIHOOD???
		out4 <- res2_pweGH$value
		
#		Give init.values correct naming 
		names(init.values) <- rownames(out1)
		
		
		
	}
	
	#############################
	### GUMBEL-HOUGAARD & COX ###
	#############################
	
	if(marginal=="Cox" & copula=="GH"){
		
		if(verbose){cat("Stage 1 initiated. \n")}
		
		
		eval(parse(text=paste0("temp_formula <- Surv(",time,",",status,") ~ cluster(",clusters,")")))
		
		for(i.covariate in covariates){
		  eval(parse(text=paste0("temp_formula <- update(temp_formula, ~ . + ",i.covariate,")")))
		}

		PH_COX <- coxph(temp_formula,data=data)
		betas2_semipar <- as.numeric(PH_COX$coeff) 
		stderrbetas2_semipar <- sqrt(PH_COX$var) 
		
		newdata.zero <- as.data.frame(matrix(0,nrow=1,ncol=length(covariates)))
		colnames(newdata.zero) <- names(PH_COX$coeff)
		
		fit0 <- survfit(PH_COX,newdata=newdata.zero) # put all covariates at 0
		S0 <- fit0$surv
		
		for(i in 1:length(ClusterDataList)){
			
	
			if(length(covariates)==1){
				cov.lincomb <- (betas2_semipar * ClusterDataList[[i]][,covariates])
			}else{
				
				# ALTERNATIVE 1 (use vectorised solution)
#				cov.lincomb <- apply(data[,covariates,drop=FALSE],MARGIN=1,FUN=function(row){sum(betas2_semipar * row)})
				
				# ALTERNATIVE 2 (Use List -> faster but probably more RAM)
				data.list <- split(as.matrix(ClusterDataList[[i]][,covariates,drop=FALSE]), seq(nrow(ClusterDataList[[i]])))
				cov.lincomb <- unlist(lapply(data.list,FUN=function(row){sum(betas2_semipar*row)}))
				names(cov.lincomb) <- NULL
				rm(list="data.list")
			}
			
			temp.S <- c(1:length(ClusterDataList[[i]][,time]))
			for(j in 1:length(ClusterDataList[[i]][,time])){
				temp.S[j] <- (S0[fit0$time==ClusterDataList[[i]][,time][j]])^exp(cov.lincomb[j])
			}
			
			ClusterDataList[[i]]$S <- temp.S
		}
		
		if(verbose){cat("Stage 1 finalized. \n \n")}
		if(verbose){cat("Stage 2 initiated. \n")}
		
		if(is.null(optim.method)){
			optim.method <- "Brent"
			lower <- min(c(init.values-1,0))
			upper <- init.values+1
		}
		
		#plug marginal estimates into loglikelihood 
		
		# init.values = theta
		res2_semiparGH <- optim(init.values,loglik.2stage_GH,marginal=marginal,status=status,ClusterData=ClusterData,ClusterDataList=ClusterDataList,control=list(maxit=3000),method=optim.method,lower=lower,upper=upper)
		theta2_semiparGH <- exp(res2_semiparGH$par)/(1+exp(res2_semiparGH$par))
		
		if(verbose){
			cat("Estimates: \n")
			cat(c(betas2_semipar,theta2_semiparGH),fill=1,labels=c(paste0("beta_",covariates," ="),"Theta ="))
		}
		
		init.values.jack <-  c(log(theta2_semiparGH/(1-theta2_semiparGH)))# theta
		jack_out <- jack_2stage_coxGH(data=data,covariates=covariates,status=status,time=time,clusters=clusters,ClusterData=ClusterData,ClusterDataList=ClusterDataList,init.values=init.values.jack,verbose=verbose,optim.method=optim.method,lower=lower,upper=upper)
		
		stderrtheta2_semiparGH <- sqrt((length(ClusterDataList)-2)/length(ClusterDataList)*sum((jack_out$theta-theta2_semiparGH)^2))
		
		if(verbose){cat("Stage 2 finalized. \n")}
		
		#kendall's tau
		tau2_semiparGH <- 1-theta2_semiparGH
		stderrtau2_semiparGH <- stderrtheta2_semiparGH
		
		# Standard Errors
		
		stderr_all <- c(stderrtheta2_semiparGH,stderrbetas2_semipar)
		
		# Put estimates and standard errors in dataframe
		out1 <- data.frame(Estimates=c(theta2_semiparGH,betas2_semipar),StandardErrors=c(stderr_all))
		rownames(out1) <- c("theta",paste0("beta_",covariates))
		
		out2 <- c(tau2_semiparGH,stderrtau2_semiparGH)
		names(out2) <- c("Estimates","StandardErrors")
#		
		# SHOULD COVARIANCE BE HERE?
		out3 <- NA
		##		rownames(out3) <- colnames(out3) <- rownames(out1)
		
		
		# CORRECT LIKELIHOOD???
		out4 <- res2_semiparGH$value
		
#		Give init.values correct naming 
		names(init.values) <- "theta"
		
	}
	
	
	# save all parameters
	parameters <- list(time=time,status=status,clusters=clusters,covariates=covariates,init.values=init.values,
			marginal=marginal,copula=copula,n.piecewise=n.piecewise)
	
	# Add factor base levels (if necessary)
	if(!is.null(factorbasenames)){
	  out1 <- insert_factorbase(df=out1,factorbasenames=factorbasenames)
	}
	
	out <- list(Parameters=out1,Kendall_Tau=out2,ParametersCov=out3,logllh=out4,parameter.call=parameters)
	class(out) <- "Sunclarco"
	return(out)
}


##################################################
##################################################



#' Sunclarco Model
#' 
#' @description Model for Survival Analysis of Unbalanced Clusters using Archimedes Copula's.
#' @references Prenen L, Braekers R, Duchateau L (2017). Extending the Archimedean copula methodology to model multivariate survival data grouped in clusters of variable size. Journal of the Royal Statistical Society, 6, 1-24. 
#' 
#' @export
#' @param data Input dataframe containing all variables.
#' @param time Which variable name is the time covariate?
#' @param status The status indicator, 0=alive, 1=dead. 
#' @param clusters The variable name describing the clusters.
#' @param covariates A vector of one or more covariates to be included in the model. Categorical covariates should be a factor in this data frame.
#' @param stage Denotes whether the one-stage (\code{stage=1}, default) or the two-stage (\code{stage=2}) approach should be used. See Details for more information.
#' @param copula Denotes which copula to use. Can be \code{"Clayton"} (default) or \code{"GH"} for Gumbel-Hougaard.
#' @param marginal Denotes which marginal survival function to use. Can be \code{"Weibull"} (default), \code{"PiecewiseExp"} for Piecewise Exponential or \code{"Cox"} for non-parametric.
#' @param n.piecewise For \code{marginal="PiecewiseExp"}, denotes how many pieces the Piecewise Exponential should have (Default = 20).
#' @param init.values A List object which contains the initial values for the parameters. This depends on the choice of the parameters \code{stage}, \code{copula} and \code{marginal}. See the \emph{Initial Values} Section for more information. If no initial parameters are given, they will be chosen automatically (See Details for more information). 
#' @param baselevels Denotes the level of a categorical covariate in the \code{covariates} vector to be used as baseline. If not set, the first appearing level will be used as the baseline level. The specification should be done as a character vector and the names of this vector should coincide with the chosen factor variable (e.g. \code{c(disease='Other',region='Region1')} in which \code{disease} and \code{region} are factor covariates).
#' @param verbose Print some in-between results as well as computation progress.
#' @param summary.print Logical value to print a short summary at the end of the computation.
#' @param optim.method Method used for optimization in one-stage estimation or in second stage of two-stage estimation. Can either be \code{"Nelder-Mead"}, \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"}, \code{"SANN"} or \code{"Brent"}. Default in one-stage estimation is \code{"Nelder-Mead"} with Weibull margins and \code{"BFGS"} with piecewise exponential margins. Default in two-stage estimation is \code{"Brent"}, except for the combination of Gumbel copula with Weibull margins, where the default is \code{"BFGS"}. 
#' @param optim.bounds Lower and upper bounds on the variables for the \code{"L-BFGS-B"} method, or bounds in which to search for method \code{"Brent"}. Should be a vector of length 2 in which the first element is the lower and the second the upper bound (e.g. \code{c(-Inf,Inf)}). 
#' If \code{optim.method = NULL} and \code{"Brent"} is used, then default bounds will be chosen. Otherwise, \code{optim.bounds} is defaulted to \code{c(-Inf,Inf)}.
#' @return S3 List object
#' \itemize{
#' \item \code{Parameters}: Data frame containing estimates and standard errors of parameters.
#' \item \code{Kendall_Tau}: Vector containing estimate and standard error of Kendall's Tau.
#' \item \code{ParametersCov}: If available, covariance matrix of the parameters. For 2-stage approaches this is only available for the Weibull marginal.
#' \item \code{logllh}: The log-likelihood value.
#' \item \code{parameter.call}: A list containing all arguments given to the function, as well as the initial parameter values and the elapsed time.
#' }
#' @details 
#' All copula models, regardless the choice of the marginal survival function, can be fitted with the two-stage approach.
#' The one-stage approach, however, is only available for the \code{"Weibull"} and \code{"PiecewiseExp"} marginal survival functions choice.
#' The one-stage approach is preferred as it leads to less biased estimates in the case of small sample sizes.  
#' When no initial values for the parameters are given, initial values for the optimisation procedure will be derived in the following way.
#' Initial values for the marginal survival functions are obtained by estimating the parameters marginally, i.e., without taking into consideration the copula function.
#' In the two-stage approach, these estimates are fixed, whereas in the one-stage approach, they are parameters in the optimisation.
#' The association parameter is set arbitrarily to 0.5 for \code{"Clayton"} 0.55 for \code{"GH"}.
#' An initial value for the association parameter can be supplied as \code{c(theta=value)} 
#' Initial values for the marginal survival function parameters can only be supplied for the \code{"Weibull"} choice as \code{c(lambda=value,rho=value)}
#' Initial values for the beta parameters from continuous covariates can be supplied as \code{c(beta_variablename=value)}
#' Initial values for the beta parameters from categorical covariates can be supplied as \code{c(beta_variablename_level=value)}
#' 
#' @section Initial Values:
#' Initial values are provided in a \code{list()} object as following:
#' 
#' \code{list(lambda=c(0.5),
#'            rho=0.5,
#'            theta=0.5
#'            beta=c(0.5))}
#'          
#' \strong{Not all initial values need to be provided!} If only some of the initial values are provided, all initial parameters will be estimated (see Details),
#' but the provided initial values will overwrite the generated ones.
#' 
#' Depending on the \code{stage} and \code{marginal} parameter,
#' different initial values can be provided:
#' \itemize{
#' \item \emph{One-Stage}:
#' \itemize{
#' \item Weibull Marginal
#' \itemize{
#' \item \code{lambda}: Single initial value for marginal survival function.
#' \item \code{rho}: Single initial value for marginal survival function.
#' \item \code{theta}: Single initial value for the association parameter.
#' \item \code{beta:}  Vector of multiple initial values for the continuous/categorical covariates.
#' }
#' \item Piecewise Exponential Marginal
#' \itemize{
#' \item \code{lambda}: Vector of multiple initial value for marginal survival function. The length of this vector should be the number of \code{n.piecewise} (see note down below).
#' \item \code{theta}: Single initial value for the association parameter.
#' \item \code{beta:}  Vector of multiple initial values for the continuous/categorical covariates.
#' } 
#' }
#' \item \emph{Two-Stage}:
#' \itemize{
#' \item Weibull or Cox Marginal
#' \itemize{
#' \item \code{theta}: Single initial value for the association parameter.
#' }
#' \item Piecewise Exponential Marginal
#' \itemize{
#' \item \code{lambda}: Vector of multiple initial value for marginal survival function. The length of this vector should be the number of \code{n.piecewise} (see note down below).
#' \item \code{theta}: Single initial value for the association parameter.
#' \item \code{beta:}  Vector of multiple initial values for the continuous/categorical covariates.
#' }
#' } 
#' }
#' 
#' \strong{Initial Values Boundaries}
#' \itemize{
#' \item \eqn{\lambda>0}
#' \item \eqn{\rho>0}
#' \item \eqn{\theta}:
#' \itemize{
#' \item GH Copula: \eqn{\theta>0} & \eqn{\theta<1} 
#' \item Clayton Copula: \eqn{\theta>0}
#' }
#' }
#' 
#' \strong{Note on \code{lambda} and \code{beta}}
#' 
#' For the Piecewise Exponential marginal, multiple \eqn{\lambda}'s should be provided
#' in the \code{lambda} slot as a vector. This vector can have a maximum length of the
#' number of pieces there were chosen (\code{n.piecewise}). In the scenario not all \eqn{\lambda}'s
#' are provided, only the first few \eqn{\lambda}'s are overwritten.
#' 
#' In the \code{beta} slot, as many \eqn{\beta}'s should be provided as there are 
#' \code{covariates} (as well as in the same order of the \code{covariates} parameter). If one of the covariates is a
#' categorical variable (\code{factor}), multiple \eqn{\beta}'s should be provided for
#' a single covariate (namely the number of levels minus 1). In the scenario not all \eqn{\beta}'s
#' are provided, only the first few \eqn{\beta}'s are overwritten.
#' 
#' @examples
#' \dontrun{
#' data("insem",package="Sunclarco")
#' result1 <- SunclarcoModel(data=insem,time="Time",status="Status",
#'                           clusters="Herd",covariates="Heifer",
#'                           stage=1,copula="Clayton",marginal="Weibull")
#' 
#' summary(result1)
#' 
#' result2 <- SunclarcoModel(data=insem,time="Time",status="Status",
#'                           clusters="Herd",covariates="Heifer",
#'                           stage=1,copula="GH",marginal="PiecewiseExp")
#' summary(result2)
#' 
#' 
#' result3 <- SunclarcoModel(data=kidney,time="time",status="status",
#'                           clusters="id",covariates="sex",
#'                           stage=2,copula="Clayton",marginal="Weibull")
#' 
#' summary(result3)
#' 
#' result4 <- SunclarcoModel(data=kidney,time="time",status="status",
#'                           clusters="id",covariates="sex",
#'                           stage=2,copula="Clayton",marginal="Cox")
#' 
#' summary(result4)
#' }
SunclarcoModel <- function(data,time,status,clusters,covariates,stage=1,copula="Clayton",marginal="Weibull",n.piecewise=20,init.values=NULL,baselevels=NULL,verbose=TRUE,summary.print=TRUE,optim.method=NULL,optim.bounds=NULL){
	
	## Check optimization method
	if(!(is.null(optim.method))){
		# Checking correct name of optimization method
		if(!(optim.method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))){
			stop(paste0("Optimization method can not be ",optim.method,". It should be either \"Nelder-Mead\", \"BFGS\", \"CG\", \"L-BFGS-B\", \"SANN\" or \"Brent\""),call.=FALSE)
		}}
  
  if(is.null(optim.bounds)){optim.bounds <- c(-Inf,Inf)}
  if(length(optim.bounds)!=2){stop("optim.bounds should be a vector of length 2",call.=FALSE)}
  if(class(optim.bounds)!="numeric" & class(optim.bounds)!="integer"){stop("optim.bounds should be a numeric vector",call. = FALSE)}
  lower <- optim.bounds[1]
  upper <- optim.bounds[2]
	
	### Check Dataframe
  if(class(data)!="data.frame"){stop("data should be a data.frame")}
	for(i in c(time,status,clusters,covariates)){
		if(!(i %in% colnames(data))){stop(paste0("\"",i,"\" is not a variable of the data frame."),call.=FALSE)}
	}
  
  

  ## Check for covariates and transform data if necessary
  factor_index <- which(unlist(lapply(data,FUN=class)[covariates])=="factor")
  if(length(factor_index)==0){factor_index <- NULL}
  if(!is.null(factor_index)){
    factorinfo <- factor2cont(data=data,factorcovariates=names(factor_index),baselevels=baselevels)
    data <- cbind(data[,!(colnames(data)%in%names(factor_index))],factorinfo$extra_columns)
    covariates <- setdiff(covariates,names(factor_index))
    covariates <- c(covariates,factorinfo$newcovariates)
    factorbasenames <- factorinfo$factorbasenames
  }else{
    factorbasenames <- NULL
  }
  
  
  
  # Check if initial values of correct format + check if within correct bounds
  if(!is.null(init.values)){
    if(class(init.values)!="list"){stop("init.values should be a list object")}
    
    if(marginal=="PiecewiseExp"){
      if(!all(names(init.values)%in%c("lambda","theta","beta"))){
        stop("init.values should only contain one or more of the following elements: 'lambda', 'theta' and 'beta'")
      }
      if("lambda"%in%names(init.values)){
        if(length(init.values$lambda)>n.piecewise){stop("lambda contains too many values in the init.values parameter")}
        if(length(init.values$lambda)<n.piecewise){warning("lambda does not contain enough values in the init.values parameter. Only part of the automatically generated initial lambda values will be overwritten.")}
      }
      if("beta"%in%names(init.values)){
        if(length(init.values$beta)>length(covariates)){stop("beta contains too many values in the init.values parameter")}
        if(length(init.values$beta)<length(covariates)){warning("beta does not contain enough values in the init.values parameter. Only part of the automatically generated initial beta values will be overwritten.")}
      }
      
    }else{
      if(stage==1){
        if(!all(names(init.values)%in%c("lambda","rho","theta","beta"))){
          stop("init.values should only contain one or more of the following elements: 'lambda', 'rho', 'theta' and 'beta'")
        }
        if("beta"%in%names(init.values)){
          if(length(init.values$beta)>length(covariates)){stop("beta contains too many values in the init.values parameter")}
          if(length(init.values$beta)<length(covariates)){warning("beta does not contain enough values in the init.values parameter. Only part of the automatically generated initial beta values will be overwritten.")}
        }
      }else if(stage==2){
        if(!all(names(init.values)%in%c("theta"))){
          stop("init.values should only contain one or more of the following elements: 'theta'")
        }
      }
    }
    
    
    # Check if parameters in init.values are within correct bounds
    for(i in c("lambda","rho")){
      if(i %in% names(init.values)){
        if(any(init.values[[i]]<=0)){stop(paste0(i," parameter should be strictly larger than 0."))}
      }
    }
    if("theta" %in% names(init.values)){
      if(copula=="GH"){
        if(init.values$theta[1]<=0 | init.values$theta[1]>=1){stop("theta parameter should be between 0 and 1 for copula=\"GH\".")}
      }else{
        if(init.values$theta[1]<=0){stop("theta parameter should be strictly larger than 0.")}
      }
    }
    
  }

  
  
	proc_start <- proc.time()
	
	# Check if we can add more stuff here (taken from the other functions)
	
	if(stage==1){
		out <- CopulaModel_1stage(data=data,time=time,status=status,clusters=clusters,covariates=covariates,init.values=init.values,marginal=marginal,copula=copula,n.piecewise=n.piecewise,verbose=verbose,factorbasenames=factorbasenames,optim.method=optim.method,lower=lower,upper=upper)
			
	}
	
	if(stage==2){
		out <- CopulaModel_2stage(data=data,time=time,status=status,clusters=clusters,covariates=covariates,init.values=init.values,marginal=marginal,copula=copula,n.piecewise=n.piecewise,verbose=verbose,factorbasenames=factorbasenames,optim.method=optim.method,lower=lower,upper=upper)
		
	}
	
	proc_end <- proc.time() - proc_start
	runtime <- proc_end['elapsed']/60  #in min
	names(runtime) <- NULL
	
	# out$parameter.call$elapsedtime <- runtime
	out$info <- list(parameters=match.call(),init.values=out$parameter.call$init.values,runtime=runtime)
	out$parameter.call <- NULL
	
	# Summary Print
	if(summary.print){summary(out)}
	
	return(out)
}



