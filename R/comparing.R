# Project: Copula Package
# 
# Author: lucp8394
###############################################################################


## 2 - STAGE CLAYTON AND WEIBULL (compare with 1-stage Clayton and Weibull)

##compare with
#I.inv <- VarCov1_weibCL
#I <- solve(I.inv)
#I11 <- I[1:3,1:3]
#I21 <- I[4,1:3]
#I12 <- I[1:3,4]
#I22 <- I[4,4]
#stderrtheta1check1 <- sqrt(solve(I22-I21%*%solve(I11)%*%I12))
#stderrtheta1check2 <- sqrt(1/I22+(I21%*%I.inv[1:3,1:3]%*%I12)/(I22^2))
#stderrtheta1_weibCL
#stderrtheta1check1
#stderrtheta1check2
#vartheta2_weibCLa <- 1/I22+I21%*%IVI%*%I12/(I22^2) #I12 en I22 komen nog uit one-stage 
#stderrtheta2_weibCLa <- sqrt(vartheta2_weibCLa) #0.02129424
##large difference because
#I
#II_weibCL
##are different
