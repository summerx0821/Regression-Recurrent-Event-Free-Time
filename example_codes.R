
#An example codes for univariate analysis of treatment

library(tidyr)
library(survival)
library(geeM)
library(gdata)

setwd("~/Desktop/")
source("PO_MI functions.R")
data_1 <- read.csv("example_data.csv")

#Get some outcomes
n <- nrow(data_1) #number of subjects
tr <- data_1$treatment #treatment
X <- data_1$follow_up_days #time to death
delta_X <- data_1$death #death indicator
R_star <- as.matrix(data_1[,paste0("recurrent",1:11)]) #recurrent events

#Specify paramters of window framework
space=60 #space between windows
tau=182 #length of window
t=seq(from=0, to=max(X)-tau, by=space) #window starting times
b=length(t) #number of windows


#Transform data frame
#get T(t) for each window
get_X_t <- function(n,b,t,R_star,X,delta_X){
  X_t=array(NA,c(n,b))
  delta=array(NA,c(n,b))
  for(j in 1:b)
  {
    R_star_j=R_star-t[j]
    R_star_j[R_star_j<=0]=NA #recurrent events observed before time t[j]
    X_j=X-t[j]
    X_j[X_j<=0]=NA #death observed before time t[j]
    R_star_j_prime=ifelse(is.na(X_j), NA, 
                          apply(cbind(X_j,R_star_j), 1, 
                                function(x) {ifelse(sum(!is.na(x))>0, min(x, na.rm = T), NA)})) #time to next recurrent event
    X_t[,j]=R_star_j_prime
    delta[,j]=ifelse(R_star_j_prime==X_j, delta_X, 1)
  }
  return(list(X_t=X_t, delta=delta))
}
X_t <- get_X_t(n,b,t,R_star,X,delta_X)$X_t
delta <- get_X_t(n,b,t,R_star,X,delta_X)$delta


############################################################################################
#PO
############################################################################################

#Format to long data - get PO within each window
min_T_tau <- apply(X_t, c(1,2), function(x) min(x, tau))
data_wide <- data.frame(cbind(1:n, tr, min_T_tau))
colnames(data_wide) <- c("ID", "tr", paste0("min_T_tau_",1:b))
#Transform from wide to long
data_long <- gather(data_wide, tj, min_T_tau, paste0("min_T_tau_",1:b), factor_key=TRUE)
data_long <- data_long[!is.na(data_long$min_T_tau),] #remove missing
#Get delta
delta_long=array(delta)
data_long$delta=delta_long[!is.na(delta_long)]
data_long$delta=ifelse(data_long$min_T_tau==tau,0,data_long$delta)


#Get the pseudo observation for each window
pseudoval=NULL
for (j in 1:b){
  subdata=data_long[data_long$tj==paste0("min_T_tau_",j),]
  pseudoval_each=get_log_pseudo_obs(X_po=subdata$min_T_tau, delta_po=subdata$delta, tau_po=tau)
  pseudoval=c(pseudoval, pseudoval_each)
}
data_long$pseudoval <- pseudoval

#Get log value and remove missing
data_long$log_min_T_tau <- data_long$pseudoval
data_long <- data_long[!is.na(data_long$log_min_T_tau),] #remove those negative pseudoval


#toeplitz
toep <- matrix(4,b,b)
diag(toep)<- 1
for(p in 1:b){for (q in 1:b){if (abs(p-q)==1){toep[p,q]=2}}}
for(p in 1:b){for (q in 1:b){if (abs(p-q)==2){toep[p,q]=3}}}

g <- summary(geem(log_min_T_tau~tr,
                  data=data_long, id=ID,
                  family=gaussian(link = "identity"), 
                  corstr = "userdefined", corr.mat = toep))
results <- cbind(Est=g$beta,
                 fold=exp(g$beta),
                 robust.se=g$se.robust,
                 robust.LB=exp(g$beta-1.96*g$se.robust),
                 robust.UB=exp(g$beta+1.96*g$se.robust),
                 p=g$p)
results
#            Est      fold  robust.se robust.LB robust.UB          p
# [1,] 4.5184546 91.693788 0.03141124 86.218830 97.516411 0.00000000
# [2,] 0.1216923  1.129407 0.04468999  1.034688  1.232796 0.00646856



############################################################################################
#MI
############################################################################################
Z=cbind(1, tr)
po_beta=results[,1] #get beta from PO results
ID=1:n
M=10
#Get imputed data from function
set.seed(821)
imputed_data_list=getMI(n=n,M=M,ID=ID,Z=Z,X=X,delta=delta,T_t=X_t, po_beta=po_beta)

#Analysis of the M multiply imputed datasets
beta_unstr=beta_defin=array(NA,c(M,ncol(Z)))
var_beta_robust_unstr=var_beta_robust_defin=array(NA,c(M,ncol(Z)))
for(m in 1:M)
{
  #Transformed variables
  imputed_T_t=imputed_data_list[[m]]
  imputed_min_T_tau <- apply(imputed_T_t, c(1,2), function(x) min(x, tau))
  imputed_data_wide <- data.frame(cbind(ID, Z[,-1], imputed_min_T_tau))
  colnames(imputed_data_wide) <- c("ID","tr", paste0("min_T_tau_",1:b))
  #Transform from wide to long
  imputed_data_long <- gather(imputed_data_wide, tj, min_T_tau, paste0("min_T_tau_",1:b), factor_key=TRUE)
  imputed_data_long <- imputed_data_long[!is.na(imputed_data_long$min_T_tau),] #remove missing
  imputed_data_long$log_min_T_tau <- log(imputed_data_long$min_T_tau)
  imputed_data_long <- imputed_data_long[!is.na(imputed_data_long$log_min_T_tau),] #remove those negative pseudoval
  
  #defined
  g <- summary(geem(log_min_T_tau~tr,
                    data=imputed_data_long, id=ID,
                    family=gaussian(link = "identity"), 
                    corstr = "userdefined", corr.mat = toep))
  beta_defin[m,] <- g$beta
  var_beta_robust_defin[m,] <- g$se.robust^2
  print(paste("Get beta for impute",m))
}

#Get MI beta and variance
beta_MI=apply(beta_defin,2,mean)
W_robust=apply(var_beta_robust_defin,2,mean)
B=diag(cov(beta_defin))
se_beta_MI_robust=sqrt(W_robust + (1+1/M)*B)
v=(M-1)*(1+W_robust/(B*(M+1)))^2

results <- cbind(Est=exp(beta_MI),
                 robust.se=se_beta_MI_robust,
                 robust.LB=exp(beta_MI-1.96*se_beta_MI_robust),
                 robust.UB=exp(beta_MI+1.96*se_beta_MI_robust),
                 p=ifelse(beta_MI/se_beta_MI_robust>0,
                          (1-pt(beta_MI/se_beta_MI_robust,df=v))*2,
                          pt(beta_MI/se_beta_MI_robust,df=v)*2))
results
#            Est  robust.se robust.LB robust.UB          p
# [1,] 92.940349 0.03055103 87.538426 98.675620 0.00000000
# [2,]  1.128984 0.04339621  1.036927  1.229214 0.00519133



