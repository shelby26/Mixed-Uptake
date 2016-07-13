##Resampling Stats for Frequency of Zeros in Mixed infections
##If I took a random sample from the data (how much to pool?)
##how many zeros would I get by chance? Did I observe more?

##The most cumbersome way to separate and pool data ever.
AB_onlyT1_Light=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==1&Env=="Light")
AD_onlyT1_Light=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==1&Env=="Light")
BD_onlyT1_Light=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==1&Env=="Light")
AB_onlyT1_Dark=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==1&Env=="Dark")
AD_onlyT1_Dark=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==1&Env=="Dark")
BD_onlyT1_Dark=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==1&Env=="Dark")
AB_onlyT2_Light=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==2&Env=="Light")
AD_onlyT2_Light=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==2&Env=="Light")
BD_onlyT2_Light=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==2&Env=="Light")
AB_onlyT2_Dark=subset(qPCR,qPCR$Sym_ID=="AB"&qPCR$D_Count==0&qPCR$Time==2&Env=="Dark")
AD_onlyT2_Dark=subset(qPCR,qPCR$Sym_ID=="AD"&qPCR$B_Count==0&qPCR$Time==2&Env=="Dark")
BD_onlyT2_Dark=subset(qPCR,qPCR$Sym_ID=="BD"&qPCR$A_Count==0&qPCR$Time==2&Env=="Dark")

Pool_1_Light=c(AB_onlyT1_Light$A_Count,AD_onlyT1_Light$A_Count,BD_onlyT1_Light$B_Count,AB_onlyT1_Light$B_Count,AD_onlyT1_Light$D_Count,BD_onlyT1_Light$D_Count)
Pool_1_Dark=c(AB_onlyT1_Dark$A_Count,AD_onlyT1_Dark$A_Count,BD_onlyT1_Dark$B_Count,AB_onlyT1_Dark$B_Count,AD_onlyT1_Dark$D_Count,BD_onlyT1_Dark$D_Count)
Pool_2_Light=c(AB_onlyT2_Light$A_Count,AD_onlyT2_Light$A_Count,BD_onlyT2_Light$B_Count,AB_onlyT2_Light$B_Count,AD_onlyT2_Light$D_Count,BD_onlyT2_Light$D_Count)
Pool_2_Dark=c(AB_onlyT2_Dark$A_Count,AD_onlyT2_Dark$A_Count,BD_onlyT2_Dark$B_Count,AB_onlyT2_Dark$B_Count,AD_onlyT2_Dark$D_Count,BD_onlyT2_Dark$D_Count)

##Just need this script for the pooled mean and CI, then count zeros in individual data and compare within light/time. 
#overly complex to include the sym1,sym2 since they are not used for the pool, but are for the sample size...
pZero_pool<-function(n_obs,pool,nsim)
{

  a=rep(0,n_obs)
  a_sumzero=rep(0,nsim)
  
  for(i in 1:nsim)
  {
    a=sample(pool,n_obs,replace=T)
    a_sumzero[i]=sum(a==0)
    ##cat("\n","\n",a_sumzero)
  }
  ##cat(a_sumzero)
  perm_mean=mean(a_sumzero)
  perm_SD=sd(a_sumzero)
  up_qt=quantile(a_sumzero,0.95)
  error<-qnorm(0.95)*perm_SD/sqrt(n_obs)
  #answer1=obs1_zero>up_qt
  #answer2=obs2_zero>up_qt
  cat("\n Average zeros by chance =",perm_mean,"\n Upper CI =",up_qt)
  hist(a_sumzero)
}