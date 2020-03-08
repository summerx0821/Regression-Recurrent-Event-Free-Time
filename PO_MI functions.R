
#############################################################
#Functions to get pseudo observations or multiple imputations
#############################################################


#Get PO (log version)
get_E_log=function(X_po,delta_po,tau_po)
{
  time=survfit(Surv(X_po, delta_po)~1)$time
  surv=survfit(Surv(X_po, delta_po)~1)$surv
  
  a=c(0,time)
  K=length(a)-1
  d=a[2:(K+1)]-a[1:K]
  survival=c(1,surv)
  dsurvival=survival[2:(K+1)]-survival[1:K]
  #Get Pseudo obs
  delta_hat=-sum(log(a[1:K]+d/2)*dsurvival) + log(tau_po)*min(surv[time<=tau_po])
  return(delta_hat)
}

get_log_pseudo_obs=function(X_po,delta_po,tau_po){
  n_po=length(X_po)
  all_obs_results=get_E_log(X_po,delta_po,tau_po)
  pseudo_obs=array()
  for(k in 1:n_po)
  {
    remove_k_results=get_E_log(X_po[-k],delta_po[-k],tau_po)
    pseudo_obs[k]=n_po*all_obs_results-(n_po-1)*remove_k_results
  }
  return(pseudo_obs)
}



#Get MI
sum_function=function(X)
{
  temp=array(X,c(length(X),length(X)))
  lowerTriangle(temp,diag=F)=0
  apply(temp,2,sum,na.rm=TRUE)
}

getMI <- function(n,M,ID,Z,X,T_t,delta,po_beta){
  #Identify patients needing imputation
  impute_subjects=ID[round(X,digits = 1)<(tau+t[length(t)]) & apply(delta,1,function(x){x[sum(!is.na(x))]})==0] #Only impute those who loss of follow up, and the last window at risk captures the censoring
  n_impute_subjects=length(impute_subjects)
  imputed_data=T_t
  imputed_data[impute_subjects,]=imputed_data[impute_subjects,]*delta[impute_subjects,] #anywhere 0 needs to be imputed
  imputed_data_list=rep(list(imputed_data),M)
  
  #Identify times when dN>0
  for(k in 1:n_impute_subjects)
  {
    #Identify risk set for patient k
    patient_id=impute_subjects[k]
    largest_possible_t=max(t[t<X[patient_id]]) #time of last window for patient k
    risk_set=rep(0,n)
    j=0
    epsilon=0.01
    while(sum(risk_set)<5)
    {
      epsilon=epsilon+j*0.001
      risk_set[X>X[patient_id] & T_t[,which(t==largest_possible_t)]>T_t[patient_id,which(t==largest_possible_t)] 
               & abs(t(t(po_beta)%*%t(Z)) - sum(po_beta*Z[patient_id,]))<epsilon]=1
      j=j+1
    }
    
    X_riskset=T_t[risk_set==1,which(t==largest_possible_t)]
    delta_riskset=delta[risk_set==1,which(t==largest_possible_t)]
    
    n_riskset=length(X_riskset)
    time_riskset=sort(unique(X_riskset))
    n_time_riskset=length(time_riskset)
    
    #Estimate survival - NA estimator
    temp3=array(X_riskset,c(n_riskset,n_time_riskset))
    temp4=t(array(time_riskset,c(n_time_riskset,n_riskset)))
    delta_array_riskset=array(delta_riskset,c(n_riskset,n_time_riskset))
    dN_T_riskset=array(as.numeric(temp3==temp4 & delta_array_riskset==1),c(n_riskset,n_time_riskset))
    Y_riskset=array(as.numeric(temp3>=temp4),c(n_riskset,n_time_riskset))
    
    denominator=apply(Y_riskset,2,sum,na.rm=TRUE)
    temp5=apply(t(dN_T_riskset)/denominator,1,sum,na.rm=TRUE)
    CH_Rk=sum_function(temp5)
    surv_Rk=exp(-CH_Rk)
    time_Rk=time_riskset
    #plot(time_Rk,surv_Rk,type="s",col="red")
    
    
    #Sampling a valid impute
    l=1
    subject_imputes=NULL
    while(length(subject_imputes)<M & l<=100)
    {
      u=runif(1)
      if(u<=min(surv_Rk)){subject_imputes_one=largest_possible_t+tau}
      if(u>min(surv_Rk)){
        impute_index=(length(time_Rk)-apply(array(surv_Rk,c(length(surv_Rk),1))<=t(array(u,c(1,length(surv_Rk)))),2,sum)) + 1
        imputation_times=time_Rk[impute_index]
        
        temp=array(X_riskset,c(n_riskset,1))==t(array(time_Rk[impute_index],c(1,n_riskset)))
        risk_set_imputes_subjects=array(ID[risk_set==1],c(length(ID[risk_set==1]),1))[temp]
        residuals=log(apply(cbind(T_t[risk_set_imputes_subjects,which(t==largest_possible_t)],rep(tau,1)),1,min))-sum(po_beta*Z[risk_set_imputes_subjects,]) 
        subject_imputes_one=exp(sum(po_beta*Z[patient_id,]) + residuals) + largest_possible_t
      }
      if(subject_imputes_one[1]>X[patient_id])
      {
        subject_imputes=c(subject_imputes,subject_imputes_one[1])
      }
      print(paste("Imputing",l))
      l=l+1
    }
    #If don't have enough imputes after 100 sampling use the largest censoring time
    # subject_imputes <- subject_imputes[1:M]
    # subject_imputes <- c(subject_imputes, rep(largest_possible_t+tau,M-length(subject_imputes)))
    
    #Put imputes in each data
    for (m in 1:M){
      for (ti in 1:which(t==largest_possible_t)){
        if (imputed_data_list[[m]][patient_id,ti]==0){
          imputed_data_list[[m]][patient_id,ti]=subject_imputes[m]-(ti-1)*space
        }
      }
    }
    print(paste("Patient",k))
  }
  return(imputed_data_list)
}



