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
pZero<-function(sym1,sym2,pool,nsim)
  {
  sim_cct_all=pool #pool data
  samp_size=length(sym1) #create sample size
  a=rep(0,samp_size)
  a_sumzero=rep(0,nsim)
  
  for(i in 1:nsim)
    {
    a=sample(sim_cct_all,samp_size,replace=T)
    a_sumzero[i]=sum(a==0)
    ##cat("\n","\n",a_sumzero)
    }
  ##cat(a_sumzero)
  perm_zero=mean(a_sumzero)
  obs1_zero=sum(sym1==0)
  obs2_zero=sum(sym2==0)
  up_qt=quantile(a_sumzero,0.95)
  answer1=obs1_zero>up_qt
  answer2=obs2_zero>up_qt
  cat("\n Average zeros by chance =",perm_zero,"\n",up_qt)
  cat("\n Observed zeros in 1=",obs1_zero,"significant?",answer1,"\n")
  cat("\n Observed zeros in 2=",obs2_zero,"significant?",answer2,"\n")
  }