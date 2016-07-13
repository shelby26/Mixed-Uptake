##Resampling test for qPCR data to look for dominance i.e., avoidance of 50:50 strain

##simulated samples should come from a common pool

qPCR_percentdom<-function(cct1,cct2,CI,nsim) {   ##input is raw cell number (cell count) data for two mixed strains, CI (as decimal), # simulations

##upper=CI+((1-CI)/2)					##calcualte fraction for upper CI limit
upper=CI
##lower=((1-CI)/2)
lower=0					##calculate fraction for lower CI limit
sim_cct_all=c(cct1,cct2)
sim_cct1=cct1
sim_cct2=cct2
num_pairs=length(cct1)
samp_num1=rep(0,num_pairs)			##establish vector to store ramdomization numbers
samp_num2=rep(0,num_pairs)			##establish vector to store ramdomization numbers
sim_diff=rep(0,nsim)			##establish vector to store simulated paired difference values
sim_diff_mean=rep(0,nsim)
emp_diff=rep(0,num_pairs)
cnt=length(cct1)

	emp_diff=abs((cct1/(cct1+cct2))-(cct2/(cct1+cct2)))
	##cat(emp_diff)
	mean_emp_diff=mean(emp_diff)
	##cat(emp_diff)
	cat(" data mean=",mean_emp_diff)

for(i in 1:nsim)					##make loop to sample from vectors samp_num times
	{
	sim_num1=sample(sim_cct_all,num_pairs,replace=T)	##simulate random position to draw samples from data with replacement
##	sim_num1=sample(sim_cct1,num_pairs,replace=T)	##simulate random position to draw samples from data with replacement

	sim_num2=sample(sim_cct_all,num_pairs,replace=T)
	##sim_num2=sample(sim_cct2,num_pairs,replace=T)

	sim_diff=abs((sim_num1/(sim_num1+sim_num2))-(sim_num2/(sim_num1+sim_num2)))
	sim_diff_mean[i]=mean(sim_diff)
	##cat("\n","\n",sim_diff)
	}
perm_mean=mean(sim_diff_mean,na.rm=T)
##cat("\n","\n",perm_mean)
up_qt=quantile(sim_diff_mean,0.95,na.rm=T)
low_qt=quantile(sim_diff_mean,0.0,na.rm=T)
#cat("\n","\n",perm_mean)
cat("\n"," perm mean ",perm_mean,"\n",upper," CI = ",up_qt,"\n",lower," CI = ",low_qt,"\n")
answer=mean_emp_diff>up_qt|mean_emp_diff<low_qt ##determine if observed mean is outside of CI
cat("\n n=",num_pairs,"observed mean outside CI?",answer,"\n")
}
